/* Assumes parameter list is opened and reads the things it needs. 
   Variables are initialized etc, default values are set */


#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h> /* for "system" command */
#include <strings.h>

void set_convection_defaults(E)
     struct All_variables *E;
{
    void PG_timestep_with_melting();
    void PG_timestep();
    void read_convection_settings();
    void convection_derived_values();
    void convection_allocate_memory();
    void convection_boundary_conditions();
    void node_locations();
    void convection_initial_fields();
    void twiddle_thumbs();
 
    E->next_buoyancy_field = PG_timestep;
    E->special_process_new_buoyancy = twiddle_thumbs; 
    E->problem_settings = read_convection_settings;
    E->problem_derived_values = convection_derived_values;
    E->problem_allocate_vars = convection_allocate_memory;
    E->problem_boundary_conds = convection_boundary_conditions;
    E->problem_initial_fields = convection_initial_fields;
    E->problem_node_positions = node_locations;
    E->problem_update_node_positions = twiddle_thumbs;
    E->problem_update_bcs = twiddle_thumbs;

    sprintf(E->control.which_data_files,"Temp,Strf,Pres");
    sprintf(E->control.which_horiz_averages,"Temp,Visc,Vrms");
    sprintf(E->control.which_running_data,"Step,Time,");
    sprintf(E->control.which_observable_data,"Shfl");
 
return;
}

void read_convection_settings(E)
     struct All_variables *E;
    
{ 
    void advection_diffusion_parameters();
    
/* parameters */

    input_float("rayleigh",&(E->control.Ra_temp),"essential");

    input_float("density_jump",&(E->control.Ra_comp),"essential");
    
    input_boolean("halfspace",&(E->convection.half_space_cooling),"off");
    input_float("halfspage",&(E->convection.half_space_age),"nodefault");
    
    input_int("temperature_blobs",&(E->convection.temp_blobs),"0");
    input_float_vector("temperature_blobx",E->convection.temp_blobs,E->convection.temp_blob_x);
    input_float_vector("temperature_bloby",E->convection.temp_blobs,E->convection.temp_blob_y);
    input_float_vector("temperature_blobz",E->convection.temp_blobs,E->convection.temp_blob_z);
    input_float_vector("temperature_blobsize",E->convection.temp_blobs,E->convection.temp_blob_radius);
    input_float_vector("temperature_blobDT",E->convection.temp_blobs,E->convection.temp_blob_T);
    input_float_vector("temperature_blobbg",E->convection.temp_blobs,E->convection.temp_blob_bg);
    input_int_vector("temperature_blobsticky",E->convection.temp_blobs,E->convection.temp_blob_sticky);
    
    input_int("temperature_zones",&(E->convection.temp_zones),"0");
    input_float_vector("temperature_zonex1",E->convection.temp_zones,E->convection.temp_zonex1);
    input_float_vector("temperature_zonex2",E->convection.temp_zones,E->convection.temp_zonex2);
    input_float_vector("temperature_zonez1",E->convection.temp_zones,E->convection.temp_zonez1);
    input_float_vector("temperature_zonez2",E->convection.temp_zones,E->convection.temp_zonez2);
    input_float_vector("temperature_zoney1",E->convection.temp_zones,E->convection.temp_zoney1);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zoney2",E->convection.temp_zones,E->convection.temp_zoney2);
    input_float_vector("temperature_zonehw",E->convection.temp_zones,E->convection.temp_zonehw);
    input_float_vector("temperature_zonemag",E->convection.temp_zones,E->convection.temp_zonemag);
    input_int_vector("temperature_zonesticky",E->convection.temp_zones,E->convection.temp_zone_sticky);
    
    input_int("num_perturbations",&(E->convection.number_of_perturbations),"0,0,32");
    input_float_vector("perturbmag",E->convection.number_of_perturbations,E->convection.perturb_mag);
    input_float_vector("ll",E->convection.number_of_perturbations,E->convection.perturb_ll);
    input_float_vector("mm",E->convection.number_of_perturbations,E->convection.perturb_mm);
    
    input_string("prevT",E->convection.old_T_file,"initialize");
    
    advection_diffusion_parameters(E);
    
    if (E->control.restart)    {
       input_int("restart_frame",&(E->control.restart_frame),"0");
       input_int("restart_timesteps",&(E->monitor.solution_cycles),"0");
       input_float("restart_time",&(E->monitor.elapsed_time),"0.0");
       }

    return;
}

/* =================================================================
   Any setup which relates only to the convection stuff goes in here
   ================================================================= */

void convection_derived_values(E)  
     struct All_variables *E;
 
{ 

return;
}

void convection_allocate_memory(E)
     struct All_variables *E;

{ void advection_diffusion_allocate_memory();

  advection_diffusion_allocate_memory(E);

return;
}

/* ============================================ */
    
void convection_initial_fields(E)
     struct All_variables *E;

{ 
    void convection_initial_temperature();

    report(E,"convection, initial temperature");
    convection_initial_temperature(E);

  return; }

/* =========================================== */

void convection_boundary_conditions(E)
     struct All_variables *E;

{
    void velocity_boundary_conditions();
    void temperature_boundary_conditions();
    void temperatures_conform_bcs();
    void composition_boundary_conditions();
 
    velocity_boundary_conditions(E);      /* universal */
    temperature_boundary_conditions(E);

    temperatures_conform_bcs(E);

    composition_boundary_conditions(E);

    return;
}

/* ===============================
   Initialization of fields .....
   =============================== */

void convection_initial_temperature(E)
     struct All_variables *E;
{
    int i,j,k,p,node,ii,jj;
    double temp,base,radius,radius2;
    double modified_plgndr_a(),drand48();
    FILE *fp;
    void remove_horiz_ave();
    void temperatures_conform_bcs();
    void thermal_buoyancy();
    void process_restart();
    
    int in1,in2,in3,instance,nox,noy,noz,nfz,ok,noz2,ll,mm;
    char output_file[255];
    double ts,tm,t1,r1,weight,const1,para1,plate_velocity,delta_temp,age;
    double x00,x01,x02,slope,con;

    int read_previous_field();

    const int dims=E->mesh.nsd;
    const float e_5=1.0e-5;

    noy=E->mesh.noy;  
    noz=E->mesh.noz;  
    nox=E->mesh.nox;  

  if(read_previous_field(E,E->T,"temperature","Temp")==0)      {
    if(strstr(E->convection.old_T_file,"initialize") != NULL) {

        mm = E->convection.perturb_mm[0];
        ll = E->convection.perturb_ll[0];

        con = E->convection.perturb_mag[0];
        noz2 = (noz-1)/2+1;
        noz2 = (noz-1);
        con = (noz-1);

        age  = E->control.plate_age*1e6 * 3600*24*365;
        const1 = 1.0/(2.0*sqrt(E->data.therm_diff*age));

          for(i=1;i<=noy;i++)
            for(j=1;j<=nox;j++)
              for(k=1;k<=noz;k++)  {
                node=k+(j-1)*noz+(i-1)*nox*noz;
                t1=E->X[1][node];
                r1=(1.0-E->X[2][node])*E->sphere.radius;
                x00 = r1*const1;

                E->T[node] = erf(x00);

                E->C[node] = 0.0;
              E->node[node] = E->node[node] | (INTX | INTZ | INTY);
			
              }    /* close the loop for node */

       }
	
	else
	    if ( (fp = fopen(E->convection.old_T_file,"r")) == NULL) {
		fprintf(E->fp,"Error in previous temperature field given\n");fflush(E->fp);
		exit(1);  
                }
	    else {
		fprintf(E->fp,"Reading an old-format temperature file: %s\n",E->convection.old_T_file);
		fscanf(fp,"%*d %*d %*d");
		for(i=1;i<=E->mesh.nno;i++) {
		    fscanf(fp,"%f",&(E->T[i]));
		    E->node[i] = E->node[i] | (INTX | INTZ);
                    } 
		fclose(fp);
	        }
	}

   if (E->control.restart==1)
         process_restart(E);

   temperatures_conform_bcs(E);

/*
   sprintf(output_file,"%s.XandT",E->control.data_file);
   if ( (fp = fopen(output_file,"w")) != NULL) {
      for (j=1;j<=E->mesh.nno;j++)
         fprintf(fp,"X[%05d] = %.6e Z[%05d] = %.6e T[%05d] = %.6e C[%05d] = %.6e\n",j,E->X[1][j],j,E->X[2][j],j,E->T[j],j,E->C[j]);
      }        

   fclose(fp);     
*/

    thermal_buoyancy(E);

    return; 
    }
