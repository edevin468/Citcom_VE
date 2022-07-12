/*****************************************
 *   CC  III  TTTTT   CC   OO   MM MM    *
 *  C     I     T    C    O  O  M M M    *
 *  C     I     T    C    O  O  M   M    *
 *   CC  III    T     CC   OO   M   M    *
 *                                       *  
 * Developed at CIT for COnvection in    *
 * the Mantle by Louis Moresi 1992-today *
 *                                       *
 * You are free to use this code but it  * 
 * is distrubuted as BeWare i.e. it does *
 * not carry any guarantees or warranties *
 * of reliability.                       *
 *                                       *
 * Please respect all the time and work  *
 * that went into the development of the *
 * code.                                 *  
 *                                       *
 * LM                                    *
 *****************************************/
/* Functions relating to the determination of viscosity field either
   as a function of the run, as an initial condition or as specified from
   a previous file */
/* 11/18/02 Changed to use new derivation for stress dep viscosity. 
 * Old code remains as Viscosity_structures.c.bak
 */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void get_viscosity_option(E)
     struct All_variables *E;
{
    void viscosity_for_system();
 
    /* general, essential default */
  
    E->viscosity.update_allowed = 0; 
    E->viscosity.SDEPV = E->viscosity.TDEPV =  E->viscosity.YIELD_STRESS =0;
    E->viscosity.EXPX=0;
  
    input_string("Viscosity",E->viscosity.STRUCTURE,NULL);   /* Which form of viscosity */
    
    input_boolean("VISC_EQUIVDD",&(E->viscosity.EQUIVDD),"off");    /* Whether to average it */
    input_int("equivdd_opt",&(E->viscosity.equivddopt),"1");
    input_int("equivdd_x",&(E->viscosity.proflocx),"1");
    input_int("equivdd_y",&(E->viscosity.proflocy),"1");
  
    input_boolean("VISC_SMOOTH",&(E->viscosity.SMOOTH),"off");
    input_int ("visc_smooth_cycles",&(E->viscosity.smooth_cycles),"0");
    
    if ( strcmp(E->viscosity.STRUCTURE,"system") == 0) /* Interpret */ {
      fprintf(E->fp,"Viscosity derived from system state\n");
      E->viscosity.FROM_SYSTEM = 1;
      viscosity_for_system(E);
    }

  return;

}

/* ============================================ */


void viscosity_for_system(E)
     struct All_variables *E;
{
    void get_system_viscosity();
    void twiddle_thumbs();
    int i;

  /* default values .... */

   for(i=0;i<40;i++) {
       E->viscosity.N0[i]=1.0;
       E->viscosity.T[i] = 0.0;
       E->viscosity.Z[i] = 0.0;
       E->viscosity.E[i] = 0.0;
       E->viscosity.T0[i] = 0.0;
   }

  /* read in information */
    input_int("rheol",&(E->viscosity.RHEOL),"essential");
    input_int("num_mat",&(E->viscosity.num_mat),"1");
 
    input_float_vector("viscT",E->viscosity.num_mat,(E->viscosity.T));  /* redundant */
    input_float_vector("viscZ",E->viscosity.num_mat,(E->viscosity.Z));
    input_float_vector("viscE",E->viscosity.num_mat,(E->viscosity.E));
    input_float_vector("visc0",E->viscosity.num_mat,(E->viscosity.N0)); /* redundant */
    input_float_vector("shearModulus",E->viscosity.num_mat,(E->viscosity.G));
  
    input_boolean("TDEPV",&(E->viscosity.TDEPV),"on");
    input_boolean("SDEPV",&(E->viscosity.SDEPV),"off");
    input_boolean("YIELD_STRESS",&(E->viscosity.YIELD_STRESS),"off");
    input_float("cohesion",&(E->viscosity.cohesion),"nodefault");
    input_float("friction_coef",&(E->viscosity.friction_coef),"nodefault");

    E->viscosity.friction_coef = E->viscosity.friction_coef*E->sphere.radius*E->data.density*E->data.grav_acc/E->data.shear_mod;
    E->viscosity.cohesion = E->viscosity.cohesion/E->data.shear_mod;
    
    input_float("sdepv_misfit",&(E->viscosity.sdepv_misfit),"0.001");
    input_float_vector("sdepv_expt",E->viscosity.num_mat,(E->viscosity.sdepv_expt));
    input_float_vector("ref_strain_rate",E->viscosity.num_mat,(E->viscosity.sdepv_trns));
    input_float_vector("background_stress",E->viscosity.num_mat,(E->viscosity.sdepv_bg)); /* added 11/22/02 */

    for (i=0;i<E->viscosity.num_mat;i++)   {
       E->viscosity.T[i] = E->viscosity.T[i]/E->data.ref_temperature;
       E->viscosity.E[i] = E->viscosity.E[i]/(E->data.gas_const*E->data.ref_temperature);
       E->viscosity.Z[i] = E->viscosity.Z[i]*E->data.density*E->data.grav_acc*E->sphere.radius/(E->data.gas_const*E->data.ref_temperature);
       /* mult by Maxwell time for some reason */
       E->viscosity.sdepv_trns[i] = E->viscosity.sdepv_trns[i]/E->data.shear_mod;
       E->viscosity.sdepv_bg[i] = E->viscosity.sdepv_bg[i]/E->data.shear_mod;

       /*note, 2.0* may be unnecessary, is mult by 0.5 in visc_from_s
	*also, this is exactly the formula for n_ref based on refrence stress
	*even though it is called sdepv_trns (transition) and is actually
	*ref_strain_rate in input file (weird, huh?)
	*/

       if (E->viscosity.SDEPV ) 
         E->viscosity.N0[i] = E->viscosity.N0[i]
	   *exp(-E->viscosity.E[i]/(1.0+E->viscosity.T[i]));
        /*
	 removed and now included in visc_from_S()
	   *pow(E->viscosity.sdepv_trns[i],(double)(E->viscosity.sdepv_expt[i]-1.0));
	 */
       else if (E->viscosity.TDEPV)
	 /* changed for half-plate cooling model 10/28/02
         E->viscosity.N0[i] = E->viscosity.N0[i]/exp((E->viscosity.E[i]+2.0e5/E->sphere.radius*E->viscosity.Z[i])/(1.0+E->viscosity.T[i]));
	 */
	 E->viscosity.N0[i] = E->viscosity.N0[i]*exp(-E->viscosity.E[i]/(1.0+E->viscosity.T[i]));
       }
       
     
    input_boolean("TDEPV_AVE",&(E->viscosity.TDEPV_AVE),"off");
    input_boolean("VFREEZE",&(E->viscosity.FREEZE),"off");
    input_boolean("VMAX",&(E->viscosity.MAX),"off");
    input_boolean("VMIN",&(E->viscosity.MIN),"off");
    input_boolean("VISC_UPDATE",&(E->viscosity.update_allowed),"on");

     input_float("alpha",&(E->viscosity.alpha),"0.0");
     input_float("beta",&(E->viscosity.beta),"0.0");
  
    input_float("freeze_thresh",&(E->viscosity.freeze_thresh),"0.0");
    input_float("freeze_value",&(E->viscosity.freeze_value),"1.0");
    input_float("visc_max",&(E->viscosity.max_value),"nodefault");
    input_float("visc_min",&(E->viscosity.min_value),"nodefault");

    input_boolean("VISC_GUESS",&(E->viscosity.guess),"off");
    input_string("visc_old_file",E->viscosity.old_file," ");

    if(!E->viscosity.update_allowed)  {
      get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);
      }
 
return;
}


void get_system_viscosity(E,propogate,evisc,visc)
     struct All_variables *E;
     int propogate;
     float *visc,*evisc;     
{ 
    void visc_from_mat();
    void visc_from_T();
    void visc_from_S();
    void apply_viscosity_smoother();
    void add_viscoelasticity();
    void p_to_nodes();
    void v_to_nodes();

    int i,j,e;
    float *visc_old,*evisc_old;
    FILE *fp;
    char output_file[255];
    static int been=0;

    const int vpts = vpoints[E->mesh.nsd];

       for(e=1;e<=E->mesh.nel;e++)   {
          E->EVolder[e] = E->EVold[e];
          E->MaxwelltimeO[e] = E->Maxwelltime[e];
          E->EVold[e] = 0.25*(evisc[(e-1)*vpts+1]+evisc[(e-1)*vpts+2]
                             +evisc[(e-1)*vpts+3]+evisc[(e-1)*vpts+4]);
          } 
    if(E->viscosity.SDEPV || E->viscosity.YIELD_STRESS)   {
      if (been==0)  
         for(e=1;e<=E->mesh.nel;e++)   {
            E->EVolder[e] = E->EVold[e];
            E->EViO[(e-1)*vpts+1]=1.0;
            E->EViO[(e-1)*vpts+2]=1.0;
            E->EViO[(e-1)*vpts+3]=1.0;
            E->EViO[(e-1)*vpts+4]=1.0;
            }
     } 

    if(E->viscosity.TDEPV)
       visc_from_T(E,visc,evisc,propogate); 
    else 
       visc_from_mat(E,visc,evisc);   

    if(E->viscosity.SDEPV)
       visc_from_S(E,visc,evisc,propogate);

   if (been==0)  {
      sprintf(output_file,"%s/visc",E->control.data_file);
      fp=fopen(output_file,"w");
      fprintf(fp,"strain rate %g %g %g %g\n",E->viscosity.sdepv_trns[1],E->viscosity.N0[1],E->viscosity.E[1],E->viscosity.T[1]);
      for(i=1;i<=E->mesh.nel;i++)
         fprintf(fp,"%d %d %.4e\n",i,E->mat[i],evisc[ (i-1)*vpoints[E->mesh.nsd]+1 ]);
      fclose(fp);
      }
      been++;

  
/*    if (E->viscosity.SMOOTH) 
       apply_viscosity_smoother(E,visc,evisc);
*/
    if(E->viscosity.MAX) {
      for(i=1;i<=E->mesh.nel;i++)
          for(j=1;j<=vpts;j++) {
            if(evisc[(i-1)*vpts + j] > E->viscosity.max_value)
               evisc[(i-1)*vpts + j] = E->viscosity.max_value;
            }
      }

    if(E->viscosity.MIN) {
      for(i=1;i<=E->mesh.nel;i++)
          for(j=1;j<=vpts;j++)
            if(evisc[(i-1)*vpts + j] < E->viscosity.min_value)
               evisc[(i-1)*vpts + j] = E->viscosity.min_value;
      }


  if (E->monitor.solution_cycles % (E->control.record_every)  == 0)  {

//     sprintf(output_file,"%s/visc1.%d.%d",E->control.data_file,E->monitor.solution_cycles,been);
     sprintf(output_file,"%s/visc1.%d",E->control.data_file,E->monitor.solution_cycles);
      fp=fopen(output_file,"w");
      fprintf(fp,"strain rate %g %g %g %g\n",E->viscosity.sdepv_trns[1],E->viscosity.N0[1],E->viscosity.E[1],E->viscosity.T[1]);
      for(i=1;i<=E->mesh.nel;i++)
         fprintf(fp,"%d %d %.5e %5e %.4e %.4e %.4e %.4e\n",i,E->mat[i],E->Maxwelltime[i],E->Ezz[i],E->Exx[i],E->Eyy[i],E->Exz[i],evisc[ (i-1)*vpoints[E->mesh.nsd]+1 ]);
//Ezz: stress, Exx: strength from byerlee's law, Eyy: strain rate, Exz: total strain.
      fclose(fp);
     } 

   for(e=1;e<=E->mesh.nel;e++)   {
      E->viscosity.element1[e] = 0.25*(evisc[(e-1)*vpts+1]+evisc[(e-1)*vpts+2]
                       +evisc[(e-1)*vpts+3]+evisc[(e-1)*vpts+4]);
      } 

  p_to_nodes(E,E->viscosity.element1,visc,E->mesh.levmax);

      for(e=1;e<=E->mesh.nel;e++)
        for(j=1;j<=vpts;j++)
            E->viscosity.element[(e-1)*vpts+j] = evisc[(e-1)*vpts+j];


  if (E->control.pure_visc==0)
    add_viscoelasticity(E,evisc);

/*    v_to_nodes(E,evisc,visc,E->mesh.levmax);
*/
 return;
}


void apply_viscosity_smoother(E,visc,evisc)
     struct All_variables *E;
     float *visc,*evisc;

{
    void p_to_centres();
    void p_to_nodes();
    
    float  *ViscCentre;
    int i;

    ViscCentre = (float *)malloc((E->mesh.nno+10)*sizeof(float));
   
  
    for(i=1;i<=E->viscosity.smooth_cycles;i++)  {
	p_to_centres(E,visc,ViscCentre,E->mesh.levmax);
	p_to_nodes(E,ViscCentre,visc,E->mesh.levmax);
    }

 
    free ((void *)ViscCentre);
 
  return;
}

void visc_from_mat(E,Eta,EEta)
     struct All_variables *E;
     float *Eta,*EEta;
{

    int i,j,k,l,z,jj,kk,lv,ll,mm,basin,ixbasin;
    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int ends=enodes[dims];
    float x1,xn;
    double modified_plgndr_a(),con,t1;
    static int been_here=0;

    lv = E->mesh.levmax;

    for(j=1;j<=E->mesh.elz;j++)     
    for(k=1;k<=E->mesh.elx;k++)       {
      i = j + (k-1)*E->mesh.elz;
      E->Maxwelltime[i] = E->viscosity.G[E->mat[i]-1];
      }


  if (been_here==0)    {
/*    been_here ++;
*/
    for(j=1;j<=E->mesh.elz;j++)     
    for(k=1;k<=E->mesh.elx;k++)       {
      i = j + (k-1)*E->mesh.elz;
      for(jj=1;jj<=vpoints[E->mesh.nsd];jj++)
          EEta[ (i-1)*vpoints[E->mesh.nsd]+jj ]=E->viscosity.N0[E->mat[i]-1];
/*    fprintf(E->fp,"%d %.4e\n",i,EEta[ (i-1)*vpoints[E->mesh.nsd]+1 ]);
 */     }

    }

    return;
}

void visc_from_T(E,Eta,EEta,propogate)
     struct All_variables *E;
     float *Eta,*EEta;
     int propogate;

{
    void remove_horiz_ave();
    int e,i,j,k,l,z,jj,kk,imark;
    float c1,c2,c3,zero,e_6,one,eta0,Tave,depth,temp,tempa,TT[9];
    int node1,node2;
    double temp1;
    static int visits=0;
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int nel = E->mesh.nel;

    one = 1.0;
    zero = 0.0;

    switch (E->viscosity.RHEOL)   {
    case 1:
        fprintf(stderr,"not supported for rheol=1\n");
    	break;

    case 2:
        fprintf(stderr,"not supported for rheol=2\n");
    	break;

    case 3:
        if(propogate && visits==0) {
            fprintf(E->fp,"\tRheological option 3:\n");

            for(l=1;l<=E->viscosity.num_mat;l++) {
              fprintf(E->fp,"\tlayer %d/%d: E=%g T1=%g \n",
                      l,E->viscosity.num_mat,
                      E->viscosity.E[l-1],E->viscosity.T[l-1]);
            }
            fflush(E->fp);
        }


     temp =30;
     temp1 =1e10;
     for(e=1;e<=E->mesh.elx;e++)   {
       node1 = (1+e)*E->mesh.noz;
       node2 = (e)*E->mesh.noz;
       E->slice.surf[5][e] = fabs((E->U[E->id[node1].doff[2]]-E->U[E->id[node2].doff[2]])/(E->X[1][node1]-E->X[1][node2])*(E->advection.timestep+1e-8));
       }
     for(e=1;e<=E->mesh.elx;e++)   {
       E->slice.surf[5][e] = E->slice.surf[5][e]/E->viscosity.sdepv_trns[1];
       }

     for(e=1;e<=nel;e++)   {
        k = 1 + (e-1)/E->mesh.elz;
        E->Maxwelltime[e] = E->viscosity.G[E->mat[e]-1];

        for(jj=1;jj<=vpts;jj++) {
            temp=1.0e-32;
            depth=1.0e-32;
            for(kk=1;kk<=ends;kk++)   {
               temp += E->T[E->ien[e].node[kk]] * E->N.vpt[GNVINDEX(kk,jj)];
               depth += E->X[2][E->ien[e].node[kk]] * E->N.vpt[GNVINDEX(kk,jj)];
               }

            depth = 1.0-depth;

	    /*do not want depth dependence 10/29/02
            temp1 = (E->viscosity.E[E->mat[e]-1]+depth*E->viscosity.Z[E->mat[e]-1])/(temp+E->viscosity.T[E->mat[e]-1]);
	    */
	    temp1 = E->viscosity.E[E->mat[e]-1]*(1/(temp+E->viscosity.T[E->mat[e]-1]));
	    /*just use the same viscosity profile*/
	    EEta[(e-1)*vpts + jj] = E->viscosity.N0[E->mat[e]-1]*exp(temp1);

	    /* not this stuff 10/29/02
	       if (E->mat[e]==2)
	       EEta[(e-1)*vpts + jj] = E->viscosity.N0[E->mat[e]-1]*exp(temp1);
	       else
	       EEta[(e-1)*vpts + jj] = E->viscosity.N0[E->mat[e]-1]/pow((double)(1.0+E->slice.surf[5][k]),(double)(0.66667))*exp(temp1);
	    */
	}
	
     }
     
     break;
     
    }

    visits++;

  return;  }


void visc_from_S(E,Eta,EEta,propogate)
     struct All_variables *E;
     float *Eta,*EEta;
     int propogate;
{
    static int visits = 0;
    static int steps = 0;
    double temp,temp1,one,two,scale,stress_magnitude,depth,exponent1,exponent;
    
    void stress_strain_rate_2_inv();
    void strain_rate_2_inv();
    void stress_2_inv();
    int e,l,z,jj,kk;

    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int nel = E->mesh.nel;

    one = 1.0;
    two = 2.0;

//     stress_2_inv(E,E->Ezz,1);

   stress_strain_rate_2_inv(E,E->Ezz,E->Eyy,E->Exz,1,E->viscosity.iterate);


     for(e=1;e<=nel;e++)   {

	/* added sdepv_bg 11/22/02 */
        temp1 = (E->Ezz[e]+E->viscosity.sdepv_bg[E->mat[e]-1])/E->viscosity.sdepv_trns[E->mat[e]-1];
        exponent=E->viscosity.sdepv_expt[E->mat[e]-1] - one;
        scale=1.0/(1.0 + pow(temp1,exponent));

        temp=1.0e-32;
        for(kk=1;kk<=ends;kk++)
               temp += E->X[2][E->ien[e].node[kk]];
        temp = temp/ends;
        temp = 1.0-temp;

        temp1 = E->viscosity.cohesion + E->viscosity.friction_coef*temp;

        E->Exx[e]=temp1;

        if (E->Ezz[e]>temp1 && E->viscosity.YIELD_STRESS) {
          for(jj=1;jj<=vpts;jj++) {
//            E->Maxwelltime[e] = (temp1/(E->Exz[e]+1e-32))*(1-E->viscosity.beta) +E->viscosity.beta*E->ShearMO[e];
//            E->ShearMO[e] = E->Maxwelltime[e];
            EEta[(e-1)*vpts + jj] = (temp1/(E->Eyy[e]+1e-32))*(1-E->viscosity.beta)+E->viscosity.beta*E->EViO[(e-1)*vpts+ jj];
            E->EViO[(e-1)*vpts + jj] = EEta[(e-1)*vpts+ jj];
            }
        }
        else {
          for(jj=1;jj<=vpts;jj++) {
            EEta[(e-1)*vpts + jj] = scale*EEta[(e-1)*vpts + jj];
            }
          }

      }



       visits ++;

    return;  
}

void stress_strain_rate_2_inv(E,SSDOT,EEDOT,totalE,SQRT,iterate)
     struct All_variables *E;
     float *SSDOT,*EEDOT,*totalE;
     int SQRT,iterate;
{
    void get_global_shape_fn();
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;

    double xk[3][5],edot[4][4],sdot[4][4],dudx[4][4],dudxt[4][4],etdot[4][4];
    float VV[4][9],VVT[4][9];

    int e,i,p,q,n,nel,k;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->mesh.nno;
    const int vpts = vpoints[dims];

    nel = E->mesh.nel;
    if (iterate==0)
          return;

   for(e=1;e<=nel;e++) {
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,xk,2,lev);

      for(i=1;i<=ends;i++)   {
        n=E->ien[e].node[i];
        VV[1][i] = E->V[1][n];
        VV[2][i] = E->V[2][n];
        if(dims==3) VV[3][i] = E->V[3][n];
        VVT[1][i] = E->X0[1][n] + E->V[1][n];
        VVT[2][i] = E->X0[2][n] + E->V[2][n];
        if(dims==3) VVT[3][i] = E->X0[3][n] + E->V[3][n];
        }

      for(p=1;p<=dims;p++)
        for(q=1;q<=dims;q++)  {
           dudx[p][q] = 0.0;
           dudxt[p][q] = 0.0;
           }

      for(i=1;i<=ends;i++)
        for(p=1;p<=dims;p++)
           for(q=1;q<=dims;q++)  {
              dudx[p][q] += VV[p][i] * GNx.ppt[GNPXINDEX(q-1,i,1)];
              dudxt[p][q] += VVT[p][i] * GNx.ppt[GNPXINDEX(q-1,i,1)];
              }

      if (E->control.AXI)   {
         edot[3][3] = 0.0;
         etdot[3][3] = 0.0;
         for(i=1;i<=ends;i++)  {
            edot[3][3] += 2.0*VV[1][i] * xk[1][1] * E->N.ppt[GNPINDEX(i,1)];
            etdot[3][3] += 2.0*VVT[1][i] * xk[1][1] * E->N.ppt[GNPINDEX(i,1)];
            }

         sdot[3][3] = ((1-E->viscosity.alpha)*E->EVolder[e] + E->viscosity.alpha*E->EVold[e])*edot[3][3] + E->S2yy[e]*E->MaxwelltimeO[e];
         }

      for(p=1;p<=dims;p++)
        for(q=1;q<=dims;q++)  {
           edot[p][q] = dudx[p][q] + dudx[q][p];
           etdot[p][q] = dudxt[p][q] + dudxt[q][p];
         }

      sdot[1][1] = ((1-E->viscosity.alpha)*E->EVolder[e] + E->viscosity.alpha*E->EVold[e])*edot[1][1] + E->S2xx[e]*E->MaxwelltimeO[e];
      sdot[2][2] = ((1-E->viscosity.alpha)*E->EVolder[e] + E->viscosity.alpha*E->EVold[e])*edot[2][2] + E->S2zz[e]*E->MaxwelltimeO[e];
      sdot[1][2] = ((1-E->viscosity.alpha)*E->EVolder[e] + E->viscosity.alpha*E->EVold[e])*edot[1][2] + E->S2xz[e]*E->MaxwelltimeO[e];

      if (dims==2)   {
         SSDOT[e] = sdot[1][1]*sdot[1][1] + sdot[2][2]*sdot[2][2]
                  + sdot[1][2]*sdot[1][2]*2.0;
         EEDOT[e] = edot[1][1]*edot[1][1] + edot[2][2]*edot[2][2]
                  + edot[1][2]*edot[1][2]*2.0;
         totalE[e] = etdot[1][1]*etdot[1][1] + etdot[2][2]*etdot[2][2]
                  + etdot[1][2]*etdot[1][2]*2.0;
         if (E->control.AXI)  {
            SSDOT[e] += sdot[3][3]*sdot[3][3];
            EEDOT[e] += edot[3][3]*edot[3][3];
            totalE[e] += etdot[3][3]*etdot[3][3];
            }
         }

      }

    if(SQRT)
        for(e=1;e<=nel;e++)  {
            SSDOT[e] =  sqrt(0.5 *SSDOT[e]);
            EEDOT[e] =  sqrt(0.5 *EEDOT[e]);
            totalE[e] =  sqrt(0.5 *totalE[e]);
            }
    else
        for(e=1;e<=nel;e++)  {
            SSDOT[e] *=  0.5;
            EEDOT[e] *=  0.5;
            totalE[e] *=  0.5;
            }

    for(e=1;e<=nel;e++)
            EEDOT[e] /= E->advection.timestep;

    return;
  }


void stress_2_inv(E,EEDOT,SQRT)
     struct All_variables *E;
     float *EEDOT;
     int SQRT;
{
    void get_global_shape_fn();
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    
    double xk[3][5],edot[4][4],dudx[4][4];
    float VV[4][9];

    int e,i,p,q,n,nel,k;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->mesh.nno;
    const int vpts = vpoints[dims];

    const float alfa=0.5;
 
    nel = E->mesh.nel;

    for(e=1;e<=nel;e++) {
  
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,xk,2,lev);

      for(i=1;i<=ends;i++)   {
        n=E->ien[e].node[i];
        VV[1][i] = E->V[1][n];
        VV[2][i] = E->V[2][n];
        if(dims==3) VV[3][i] = E->V[3][n];
        }
      
      for(p=1;p<=dims;p++)
        for(q=1;q<=dims;q++)
           dudx[p][q] = 0.0;  
      
      for(i=1;i<=ends;i++)
        for(p=1;p<=dims;p++)
           for(q=1;q<=dims;q++)
              dudx[p][q] += VV[p][i] * GNx.ppt[GNPXINDEX(q-1,i,1)];  

      if (E->control.AXI)   {
         edot[3][3] = 0.0;
         for(i=1;i<=ends;i++)
            edot[3][3] += 2.0*VV[1][i] * xk[1][1] * E->N.ppt[GNPINDEX(i,1)];  
         }
	
      for(p=1;p<=dims;p++)
        for(q=1;q<=dims;q++)
           edot[p][q] = dudx[p][q] + dudx[q][p];   

      
    
      edot[1][1] = ((1-alfa)*E->EVolder[e] + alfa*E->EVold[e])*edot[1][1] + E->S2xx[e]*E->Maxwelltime[e];
      edot[2][2] = ((1-alfa)*E->EVolder[e] + alfa*E->EVold[e])*edot[2][2] + E->S2zz[e]*E->Maxwelltime[e];
      edot[1][2] = ((1-alfa)*E->EVolder[e] + alfa*E->EVold[e])*edot[1][2] + E->S2xz[e]*E->Maxwelltime[e];
      if (E->control.AXI)   
         edot[3][3] = ((1-alfa)*E->EVolder[e] + alfa*E->EVold[e])*edot[3][3] + E->S2yy[e]*E->Maxwelltime[e];
  
      if (dims==2)   {
         EEDOT[e] = edot[1][1]*edot[1][1] + edot[2][2]*edot[2][2]
                  + edot[1][2]*edot[1][2]*2.0; 
         if (E->control.AXI)   
            EEDOT[e] += edot[3][3]*edot[3][3];
         }


      else if (dims==3)
         EEDOT[e] = edot[1][1]*edot[1][1] + edot[1][2]*edot[1][2]*2.0
                  + edot[2][2]*edot[2][2] + edot[2][3]*edot[2][3]*2.0
                  + edot[3][3]*edot[3][3] + edot[1][3]*edot[1][3]*2.0; 

      }

    if(SQRT)
	for(e=1;e<=nel;e++)
	    EEDOT[e] =  sqrt(0.5 *EEDOT[e]);
    else
	for(e=1;e<=nel;e++)
	    EEDOT[e] *=  0.5;

    return;
}


void strain_rate_2_inv(E,EEDOT,SQRT)
     struct All_variables *E;
     float *EEDOT;
     int SQRT;
{
    void get_global_shape_fn();
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
    
    double xk[3][5],edot[4][4],dudx[4][4];
    float VV[4][9];

    int e,i,p,q,n,nel,k;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->mesh.nno;
    const int vpts = vpoints[dims];
 
    nel = E->mesh.nel;

    for(e=1;e<=nel;e++) {
  
      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,xk,2,lev);

      for(i=1;i<=ends;i++)   {
        n=E->ien[e].node[i];
          VV[1][i] = E->V[1][n];
          VV[2][i] = E->V[2][n];
          if(dims==3) VV[3][i] = E->V[3][n];
        }
      
      for(p=1;p<=dims;p++)
        for(q=1;q<=dims;q++)
           dudx[p][q] = 0.0;  
      
      for(i=1;i<=ends;i++)
        for(p=1;p<=dims;p++)
           for(q=1;q<=dims;q++)
              dudx[p][q] += VV[p][i] * GNx.ppt[GNPXINDEX(q-1,i,1)];  

      if (E->control.AXI)   {
         edot[3][3] = 0.0;
         for(i=1;i<=ends;i++)
              edot[3][3] += 2.0*VV[1][i] * xk[1][1] * E->N.ppt[GNPINDEX(i,1)];  
         }
	
      for(p=1;p<=dims;p++)
        for(q=1;q<=dims;q++)
           edot[p][q] = dudx[p][q] + dudx[q][p];   

      if (dims==2)   {
         EEDOT[e] = edot[1][1]*edot[1][1] + edot[2][2]*edot[2][2]
                  + edot[1][2]*edot[1][2]*2.0; 
         if (E->control.AXI)   
            EEDOT[e] += edot[3][3]*edot[3][3];
         }

      else if (dims==3)
         EEDOT[e] = edot[1][1]*edot[1][1] + edot[1][2]*edot[1][2]*2.0
                  + edot[2][2]*edot[2][2] + edot[2][3]*edot[2][3]*2.0
                  + edot[3][3]*edot[3][3] + edot[1][3]*edot[1][3]*2.0; 

      }

    if(SQRT)
	for(e=1;e<=nel;e++)
	    EEDOT[e] =  sqrt(0.5 *EEDOT[e]);
    else
	for(e=1;e<=nel;e++)
	    EEDOT[e] *=  0.5;

    EEDOT[e] /= E->advection.timestep;
  
    return;
}
