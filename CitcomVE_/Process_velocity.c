/*  Here are the routines which process the results of each velocity solution, and call
    the relevant output routines. At this point, the velocity and pressure fields have
    been calculated and stored at the nodes. The only properties of the velocity field
    which are already known are those required to check convergence of the iterative
    scheme and so on. */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include <stdlib.h> /* for "system" command */

#include "element_definitions.h"
#include "global_defs.h"

void process_new_velocity(E,ii)
    struct All_variables *E;
    int ii;
{ 
    void output_velo_related();
    void output_stress();
    void get_STD_topo();
    void get_CBF_topo();
    void get_temperature_half_space();
    void print_surf_topo();
    void pickup_dt();

    static int been_here=0;

    if(been_here==0) {
 	   E->monitor.length_scale = E->data.layer_km/E->mesh.layer[2]; /* km */
	   E->monitor.time_scale = pow(E->data.layer_km*1000.0,2.0)/   /* Million years */
	    (E->data.therm_diff*3600.0*24.0*365.25*1.0e6);
	   been_here++;
       } 

//    pickup_dt(E);

    if (E->control.pure_visc==0)
       update_stress_strain(E,ii);  /* stress */

    deform_grid(E);

/*    get_temperature_half_space(E);  */

    print_surf_topo(E,ii);

    if ( ii % (10*E->control.record_every)  == 0)     {

/*      get_STD_topo(E,E->slice.tpg,E->slice.tpgb,ii);  */

      output_velo_related(E,ii);     

      }

    if ( ii % (2*E->control.record_every)  == 0)     {

      output_stress(E,ii);     

      }


    return;
}

/* ===============================================   */

void get_surface_velo(E, SV)
  struct All_variables *E;
  float *SV;
  {

  int el,els,i,m,node,lev;
  char output_file[255];
  FILE *fp;

  const int dims=E->mesh.nsd;
  const int ends=enodes[dims];
  const int nno=E->mesh.nno;

  lev = E->mesh.levmax;

  m = 0;

  for (node=1;node<=nno;node++)
    if ((node-1)%E->mesh.noz==0)   {
      i = (node-1)/E->mesh.noz + 1;
        SV[(i-1)*2+1] = E->V[1][node];
        SV[(i-1)*2+2] = E->V[3][node];
      }

  return;
  }

/* ===============================================   */

void get_ele_visc(E, EV)
  struct All_variables *E;
  float *EV;
  {

  int el,j,lev;

  const int nel=E->mesh.nel;
  const int vpts=vpoints[E->mesh.nsd];

  lev = E->mesh.levmax;

  for (el=1;el<=nel;el++)   {
    EV[el] = 0.0;
    for (j=1;j<=vpts;j++)
      EV[el] +=  E->EVI[lev][(el-1)*vpts+j];

    EV[el] /= vpts;
    }

  return;
  }


void get_surf_stress(E,SXX,SYY,SZZ,SXY,SXZ,SZY)
  struct All_variables *E;
  float *SXX,*SYY,*SZZ,*SXY,*SXZ,*SZY;
  {
  int i,node,stride;

  stride = E->mesh.nsf*6;

  for (node=1;node<=E->mesh.nno;node++)
     if ( ((node-1)%E->mesh.noz)==0 )  {
        i = (node-1)/E->mesh.noz+1;
        E->stress[(i-1)*6+1] = SXX[node];
        E->stress[(i-1)*6+2] = SZZ[node];
        E->stress[(i-1)*6+3] = SYY[node];
        E->stress[(i-1)*6+4] = SXY[node];
        E->stress[(i-1)*6+5] = SXZ[node];
        E->stress[(i-1)*6+6] = SZY[node];
        }
     else if ( ((node-2)%E->mesh.noz)==0 )  {
        i = (node-2)/E->mesh.noz+1;
        E->stress[stride+(i-1)*6+1] = SXX[node];
        E->stress[stride+(i-1)*6+2] = SZZ[node];
        E->stress[stride+(i-1)*6+3] = SYY[node];
        E->stress[stride+(i-1)*6+4] = SXY[node];
        E->stress[stride+(i-1)*6+5] = SXZ[node];
        E->stress[stride+(i-1)*6+6] = SZY[node];
        }

  return;
  }
