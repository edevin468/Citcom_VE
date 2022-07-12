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
/*   Functions which solve the heat transport equations using Petrov-Galerkin
     streamline-upwind methods. The process is basically as described in Alex
     Brooks PhD thesis (Caltech) which refers back to Hughes, Liu and Brooks.  */

#include <malloc.h>
#include <sys/types.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

struct el { double gpt[9]; };

/* ============================================
   Generic adv-diffusion for temperature field.
   ============================================ */


void advection_diffusion_parameters(E)
     struct All_variables *E;

{

    void std_timestep();
    /* Set intial values, defaults & read parameters*/
    
    E->advection.temp_iterations = 2; /* petrov-galerkin iterations: minimum value. */
    E->advection.total_timesteps = 1; 
    E->advection.sub_iterations = 1;
    E->advection.last_sub_iterations = 1;
    E->advection.gamma = 0.5;
    E->advection.dt_reduced = 1.0;         

    E->monitor.T_maxvaried = 1e05;
 
    input_boolean("ADV",&(E->advection.ADVECTION),"on");
   
    input_int("minstep",&(E->advection.min_timesteps),"1");
    input_int("maxstep",&(E->advection.max_timesteps),"1000");
    input_int("maxtotstep",&(E->advection.max_total_timesteps),"1000000");
    input_float("finetunedt",&(E->advection.fine_tune_dt),"0.9");
    input_float("fixed_timestep",&(E->advection.fixed_timestep),"0.0");
    input_int("adv_sub_iterations",&(E->advection.temp_iterations),"2,2,nomax");
    input_float("maxadvtime",&(E->advection.max_dimensionless_time),"10.0");
   
    input_float("sub_tolerance",&(E->advection.vel_substep_aggression),"0.005");  
    input_int("maxsub",&(E->advection.max_substeps),"25");
  
    input_float("liddefvel",&(E->advection.lid_defining_velocity),"0.01");
    input_float("sublayerfrac",&(E->advection.sub_layer_sample_level),"0.5");            
	      
  /* allocate memory */

  return;
}

void advection_diffusion_allocate_memory(E)
     struct All_variables *E;

{ int i;

  E->Tdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
  for(i=1;i<=E->mesh.nno;i++) 
    E->Tdot[i]=0.0;

  if (E->control.composition)   {
    E->Cdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
    for(i=1;i<=E->mesh.nno;i++)
      E->Cdot[i]=0.0;
    }
    

return;
}


void PG_timestep(E)
     struct All_variables *E;
{    
    void timestep();
    void predictor();
    void corrector();
    void pg_solver();
    void remove_horiz_ave();
    void std_timestep();
    void temperatures_conform_bcs();
    void thermal_buoyancy();
    float Tmax(),T_interior1;

    int i,j,psc_pass,count,steps,iredo;
    int keep_going;

    float *DTdot, *T1, *Tdot1;

    static int loops_since_new_eta = 0;
    static int been_here = 0;
   
    DTdot= (float *)malloc((E->mesh.nno+1)*sizeof(float));
    T1= (float *)malloc((E->mesh.nno+1)*sizeof(float));
    Tdot1= (float *)malloc((E->mesh.nno+1)*sizeof(float));

    if (been_here++ ==0)    {
	  E->advection.timesteps=0;
      }

    E->advection.timesteps++;
 
    std_timestep(E);


   if (E->control.composition)         {
      for (i=1;i<=E->mesh.nno;i++)   {
         T1[i] = E->C[i];
         Tdot1[i] = E->Cdot[i];
         }
                   /* get the max temperature for old T */
      T_interior1 = Tmax(E,E->C);

      E->advection.dt_reduced = 1.0;         
      E->advection.last_sub_iterations = 1;

      do  {

         E->advection.timestep *= E->advection.dt_reduced; 

         iredo = 0;
         if (E->advection.ADVECTION) {
	     predictor(E,E->C,E->Cdot,0);
	 
	     for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++) {
	 	   pg_solver(E,E->C,E->Cdot,DTdot,E->V,E->control.comp_diff,E->CB,E->node);
		   corrector(E,E->C,E->Cdot,DTdot,0);
	       }	     
	     }

                   /* get the max temperature for new T */
         E->monitor.T_interior = Tmax(E,E->C);

         if (E->monitor.T_interior/T_interior1 > E->monitor.T_maxvaried) {
           for (i=1;i<=E->mesh.nno;i++)   {
              E->C[i] = T1[i];
              E->Cdot[i] = Tdot1[i];
              }
           iredo = 1;
           E->advection.dt_reduced *= 0.5;         
           E->advection.last_sub_iterations ++;
           }
         
         }  while ( iredo==1 && E->advection.last_sub_iterations <= 5);


    /* update temperature    */
/*
       predictor(E,E->T,E->Tdot,1);
       for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++)   {
         pg_solver(E,E->T,E->Tdot,DTdot,E->V,1.0,E->TB,E->node);
         corrector(E,E->T,E->Tdot,DTdot,1);
         }	     
*/
       }

   else  if (!E->control.composition)   {

     for (i=1;i<=E->mesh.nno;i++)   {
       T1[i] = E->T[i];
       Tdot1[i] = E->Tdot[i];
       }
                   /* get the max temperature for old T */
     T_interior1 = Tmax(E,E->T);

     E->advection.dt_reduced = 1.0;         
     E->advection.last_sub_iterations = 1;

     do  {

       E->advection.timestep *= E->advection.dt_reduced; 

       iredo = 0;
       if (E->advection.ADVECTION) {
          predictor(E,E->T,E->Tdot,1);

          for(psc_pass=0;psc_pass<E->advection.temp_iterations;psc_pass++)   {
	     pg_solver(E,E->T,E->Tdot,DTdot,E->V,1.0,E->TB,E->node);
	     corrector(E,E->T,E->Tdot,DTdot,1);
	     }	     
	  }

                   /* get the max temperature for new T */
        E->monitor.T_interior = Tmax(E,E->T);

        if (E->monitor.T_interior/T_interior1 > E->monitor.T_maxvaried) {
           for (i=1;i<=E->mesh.nno;i++)   {
             E->T[i] = T1[i];
             E->Tdot[i] = Tdot1[i];
             }
           iredo = 1;
           E->advection.dt_reduced *= 0.5;         
           E->advection.last_sub_iterations ++;
           }
         
        }  while ( iredo==1 && E->advection.last_sub_iterations <= 5);

      }        /* end of if !composition */

    E->advection.total_timesteps++;
    E->monitor.elapsed_time += E->advection.timestep; 

    temperatures_conform_bcs(E,E->T); 


    thermal_buoyancy(E);
 
    if( (
	  (E->advection.timesteps < E->advection.max_timesteps) && 
	  (E->monitor.elapsed_time < E->advection.max_dimensionless_time) ) ||
      (E->advection.total_timesteps < E->advection.min_timesteps) )
        E->control.keep_going = 1;
    else
 	    E->control.keep_going = 0;

    if (E->advection.last_sub_iterations==5)
       E->control.keep_going = 0;
    
    free((void *) DTdot );  
    free((void *) T1 );   
    free((void *) Tdot1 );
    
return;  
}


/* ==============================
   predictor and corrector steps.
   ============================== */

void predictor(E,field,fielddot,ic)
     struct All_variables *E;
     float *field,*fielddot;
     int ic;

{ 
    int node;
    float multiplier;

   multiplier = (1.0-E->advection.gamma) * E->advection.timestep;

  if (ic==1)
    for(node=1;node<=E->mesh.nno;node++)  {
	if(!(E->node[node] & (OFFSIDE | TBX | TBZ | TBY))) 
	    field[node] += multiplier * fielddot[node] ;
	fielddot[node] = 0.0;
       }
  else
    for(node=1;node<=E->mesh.nno;node++)  {
	if(!(E->node[node] & OFFSIDE )) 
	    field[node] += multiplier * fielddot[node] ;
	fielddot[node] = 0.0;
       }

   return; }

void corrector(E,field,fielddot,Dfielddot,ic)
     struct All_variables *E;
     float *field,*fielddot,*Dfielddot;
     int ic;
    
{  int node;
   float multiplier;

   multiplier = E->advection.gamma * E->advection.timestep;

 if (ic==1)
   for(node=1;node<=E->mesh.nno;node++) {
       if(!(E->node[node] & (OFFSIDE | TBX | TBZ | TBY)))
	   field[node] += multiplier * Dfielddot[node];
       fielddot[node] +=  Dfielddot[node]; 
       }
 else
   for(node=1;node<=E->mesh.nno;node++) {
       if(!(E->node[node] & OFFSIDE ))
	   field[node] += multiplier * Dfielddot[node];
       fielddot[node] +=  Dfielddot[node]; 
       }
   
   return;  
 }

/* ===================================================
   The solution step -- determine residual vector from
   advective-diffusive terms and solve for delta Tdot
   Two versions are available -- one for Cray-style 
   vector optimizations etc and one optimized for 
   workstations.
   =================================================== */


void pg_solver(E,T,Tdot,DTdot,V,diff,TBC,FLAGS)
     struct All_variables *E;
     float *T,*Tdot,*DTdot;
     float **V;
     float diff;
     float **TBC;
     unsigned int *FLAGS;
{
    void get_global_shape_fn();
    void pg_shape_fn();
    void element_residual();
    int el,e,a,i,a1;
    double xk[3][5],Eres[9];  /* correction to the (scalar) Tdot field */

    struct Shape_function PG;
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;
 
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int ends=enodes[dims];
  
    for(i=1;i<=E->mesh.nno;i++)
 	  DTdot[i] = 0.0;

    for(el=1;el<=E->mesh.nel;el++)    {

	  get_global_shape_fn(E,el,&GN,&GNx,&dOmega,xk,0,E->mesh.levmax);
	  pg_shape_fn(E,el,&PG,&GNx,V,diff);
	  element_residual(E,el,PG,GNx,dOmega,V,T,Tdot,Eres,diff);

      for(a=1;a<=ends;a++) {
	    a1 = E->ien[el].node[a];
	    DTdot[a1] += Eres[a]; 
        }

      } /* next element */

    for(i=1;i<=E->mesh.nno;i++) {
 	  if(E->node[i] & OFFSIDE) continue;
	  DTdot[i] *= E->Mass[i];         /* lumped mass matrix */
      }
     
    return;    
}



/* ===================================================
   Petrov-Galerkin shape functions for a given element
   =================================================== */

void pg_shape_fn(E,el,PG,GNx,V,diffusion)
     struct All_variables *E;
     int el;
     struct Shape_function *PG;
     struct Shape_function_dx *GNx;
     float **V;
     float diffusion;

{ 
    int i,j,node;
    int *ienmatrix;

    double uc1,uc2,uc3;
    double u1,u2,u3,VV[4][9];
    double uxse,ueta,ufai,xse,eta,fai,dx1,dx2,dx3,adiff;

    double prod1,unorm,twodiff;
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int lev=E->mesh.levmax;
    const int nno=E->mesh.nno;
    const int ends=enodes[E->mesh.nsd];
    const int vpts=vpoints[E->mesh.nsd];
  
    ienmatrix=E->ien[el].node;

    twodiff = 2.0*diffusion;
 
    uc1 =  uc2 = uc3 = 0.0;

    for(i=1;i<=ends;i++)   {
        node = ienmatrix[i];
        VV[1][i] = V[1][node];
        VV[2][i] = V[2][node];
      }
         
    for(i=1;i<=ENODES2D;i++) {
      uc1 +=  E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
      uc2 +=  E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
      }
    dx1 = 0.5*(E->X[1][ienmatrix[3]]+E->X[1][ienmatrix[4]]
              -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[2]]);
    dx2 = 0.5*(E->X[2][ienmatrix[3]]+E->X[2][ienmatrix[4]]
              -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[2]]);
    uxse = fabs(uc1*dx1*E->eco[el].centre[2]+uc2*dx2); 


    dx1 = 0.5*(E->X[1][ienmatrix[2]]+E->X[1][ienmatrix[3]]
              -E->X[1][ienmatrix[1]]-E->X[1][ienmatrix[4]]);
    dx2 = 0.5*(E->X[2][ienmatrix[2]]+E->X[2][ienmatrix[3]]
              -E->X[2][ienmatrix[1]]-E->X[2][ienmatrix[4]]);
    ueta = fabs(uc1*dx1*E->eco[el].centre[2]+uc2*dx2); 

    xse = (uxse>twodiff)? (1.0-twodiff/uxse):0.0;
    eta = (ueta>twodiff)? (1.0-twodiff/ueta):0.0;

    unorm = uc1*uc1 + uc2*uc2;

    adiff = (unorm>0.000001)?( (uxse*xse+ueta*eta)/(2.0*unorm) ):0.0;

    for(i=1;i<=VPOINTS2D;i++) {
          u1 = u2 = 0.0;
          for(j=1;j<=ENODES2D;j++)  /* this line heavily used */ {
            u1 += VV[1][j] * E->N.vpt[GNVINDEX(j,i)]; 
            u2 += VV[2][j] * E->N.vpt[GNVINDEX(j,i)];
	    }
	    
	  for(j=1;j<=ENODES2D;j++) {
             prod1 = (u1 * GNx->vpt[GNVXINDEX(0,j,i)] +
                      u2 * GNx->vpt[GNVXINDEX(1,j,i)]);
             PG->vpt[GNVINDEX(j,i)] = E->N.vpt[GNVINDEX(j,i)] + adiff * prod1;
	     }
	  }

   return;
 }



/* ==========================================
   Residual force vector from heat-transport.
   Used to correct the Tdot term.
   =========================================  */

void element_residual(E,el,PG,GNx,dOmega,V,field,fielddot,Eres,diff)
     struct All_variables *E;
     int el;
     struct Shape_function PG;
     struct Shape_function_dA dOmega;
     struct Shape_function_dx GNx;
     float **V;
     float *field,*fielddot;
     double Eres[9];
     float diff;

{
    int i,j,a,k,node,nodes[4],d,aid,back_front,onedfns;
    double Q;
    double dT[9],VV[4][9];
    double tx1[9],tx2[9],tx3[9];
    double v1[9],v2[9],v3[9];
    double adv_dT,t2[4];
    double T,DT;
 
    register double prod,sfn;
    struct Shape_function1 GM;
    struct Shape_function1_dA dGamma;
    double temp;

    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    const int vpts=vpoints[dims];
    const int diffusion = (diff != 0.0); 
 
    for(i=1;i<=vpts;i++)	{ 
      dT[i]=0.0;
      v1[i] = tx1[i]=  0.0;
      v2[i] = tx2[i]=  0.0;
      }

	for(i=1;i<=ends;i++)   {
      node = E->ien[el].node[i];
        VV[1][i] = V[1][node];
        VV[2][i] = V[2][node];
      }
         
  
    for(j=1;j<=ends;j++)       {
      node = E->ien[el].node[j];
	  T = field[node];  
	  if(E->node[node] & (TBX | TBY | TBZ))
	    DT=0.0;
	  else
	    DT = fielddot[node];
	
	    for(i=1;i<=vpts;i++)  {
	 	  dT[i] += DT * E->N.vpt[GNVINDEX(j,i)];
		  tx1[i] +=  GNx.vpt[GNVXINDEX(0,j,i)] * T; 
		  tx2[i] +=  GNx.vpt[GNVXINDEX(1,j,i)] * T;   
		  sfn = E->N.vpt[GNVINDEX(j,i)];
		  v1[i] += VV[1][j] * sfn;
		  v2[i] += VV[2][j] * sfn;
	    }
      }

   Q = 0.0;

    /* construct residual from this information */


    if(diffusion){

	    for(j=1;j<=ends;j++) {
	 	  Eres[j]=0.0;
		  for(i=1;i<=vpts;i++) 
		    Eres[j] -= 
			PG.vpt[GNVINDEX(j,i)] * dOmega.vpt[i] * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i]) +
			diff*dOmega.vpt[i] * (GNx.vpt[GNVXINDEX(0,j,i)] * tx1[i] +
					      GNx.vpt[GNVXINDEX(1,j,i)] * tx2[i] ); 

	    }
      }

    else { /* no diffusion term */
	    for(j=1;j<=ends;j++) {
    	  Eres[j]=0.0;
		  for(i=1;i<=vpts;i++) 
		    Eres[j] -= PG.vpt[GNVINDEX(j,i)] * dOmega.vpt[i] * (dT[i] - Q + v1[i] * tx1[i] + v2[i] * tx2[i]);
		  }
      }	

	/* See brooks etc: the diffusive term is excused upwinding for 
	   rectangular elements  */
  
    return; 
}




/* =====================================================
   Obtain largest possible timestep (no melt considered)
   =====================================================  */


void std_timestep(E)
     struct All_variables *E;

{ 
    static int been_here = 0;
    static float diff_timestep,root3,root2;
    int i,d,n,nel,el,node;

    float adv_timestep;
    float ts,uc1,uc2,uc3,uc,size,step,VV[4][9];
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    
	nel=E->mesh.nel;

    if(E->advection.fixed_timestep != 0.0) {
      E->advection.timestep = E->advection.fixed_timestep;
      return;
    }

    if (been_here == 0)  {
	  diff_timestep = 1.0e8; 
	  for(el=1;el<=nel;el++)  { 
	    for(d=1;d<=dims;d++)    {
	 	  ts = E->eco[el].size[d] * E->eco[el].size[d];
              if (d==1) ts = ts*E->eco[el].centre[2]*E->eco[el].centre[2];
	      diff_timestep = min(diff_timestep,ts);
	      }
	    }
      diff_timestep = 0.5 * diff_timestep;
      }


    
  adv_timestep = 1.0e8;
  for(el=1;el<=nel;el++) {

	for(i=1;i<=ends;i++)   {
      node = E->ien[el].node[i];
        VV[1][i] = E->V[1][node];
        VV[2][i] = E->V[2][node];
        if(dims==3) VV[3][i] = E->V[3][node];
      }
         
    uc=uc1=uc2=uc3=0.0;	  
    if(3==dims) {
      for(i=1;i<=ENODES3D;i++) {
        uc1 += E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
        uc2 += E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
        uc3 += E->N.ppt[GNPINDEX(i,1)]*VV[3][i];
        }
      uc = fabs(uc1)/E->eco[el].size[1] + fabs(uc2)/E->eco[el].size[2] + fabs(uc3)/E->eco[el].size[3];

	  step = (0.5/uc);
	  adv_timestep = min(adv_timestep,step);
      }
    else {
	  for(i=1;i<=ENODES2D;i++) {
        uc1 += E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
        uc2 += E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
        }
      uc = fabs(uc1)/(E->eco[el].size[1]*E->eco[el].centre[1]) + fabs(uc2)/E->eco[el].size[2];
	  
	  step = (0.5/uc);
	  adv_timestep = min(adv_timestep,step);
      }
    }

    adv_timestep = E->advection.dt_reduced * adv_timestep;         

    adv_timestep =  1.0e-32+min(E->advection.fine_tune_dt*adv_timestep,diff_timestep);
    
    E->advection.timestep = adv_timestep;


    return; 
  }

void std1_timestep(E)
     struct All_variables *E;

{ 
    static int been_here = 0;
    static float diff_timestep,root3,root2;
    int i,d,n,nel,el,node;

    float adv_timestep;
    float ts,uc1,uc2,uc3,uc,size,step,VV[4][9];
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    
	nel=E->mesh.nel;

    if(E->advection.fixed_timestep != 0.0) {
      E->advection.timestep = E->advection.fixed_timestep;
      return;
    }

    if (been_here == 0)  {
	  diff_timestep = 1.0e8; 
	  for(el=1;el<=nel;el++)  { 
	    for(d=1;d<=dims;d++)    {
	 	  ts = E->eco[el].size[d] * E->eco[el].size[d];
              if (d==1) ts = ts*E->eco[el].centre[2]*E->eco[el].centre[2];
	      diff_timestep = min(diff_timestep,ts);
	      }
	    }
      diff_timestep = 0.5 * diff_timestep;
      }


    
  adv_timestep = 1.0e8;
  for(el=1;el<=nel;el++) {

	for(i=1;i<=ends;i++)   {
      node = E->ien[el].node[i];
        VV[1][i] = E->V[1][node];
        VV[2][i] = E->V[2][node];
        if(dims==3) VV[3][i] = E->V[3][node];
      }
         
    uc=uc1=uc2=uc3=0.0;	  
    for(i=1;i<=ENODES2D;i++) {
        uc1 += E->N.ppt[GNPINDEX(i,1)]*VV[1][i];
        uc2 += E->N.ppt[GNPINDEX(i,1)]*VV[2][i];
        }
      uc = fabs(uc1)/(E->eco[el].size[1]*E->eco[el].centre[1]) + fabs(uc2)/E->eco[el].size[2];
	  
	  step = (0.5/uc);
	  adv_timestep = min(adv_timestep,step);
      }

    adv_timestep = E->advection.fine_tune_dt * adv_timestep;         

    
    E->advection.timestep = adv_timestep;


    return; 
  }

void std2_timestep(E)
     struct All_variables *E;

{ 
    static int been_here = 0;
    static float init_u,diff_timestep,root3,root2;
    int i,d,n,nel,el,node;

    float adv_timestep;
    float ts,uc1,uc2,uc3,uc,size,step,VV[4][9];
    
    const int dims=E->mesh.nsd;
    const int dofs=E->mesh.dof;
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;
    const int ends=enodes[dims];
    
    if(E->advection.fixed_timestep != 0.0) {
      E->advection.timestep = E->advection.fixed_timestep;
      return;
    }

    uc2 = 0.0;
/*
    for(node=1;node<=nno;node++) {
       uc1 = sqrt(E->V[1][node]*E->V[1][node]+E->V[2][node]*E->V[2][node]); 
       if (uc1 > uc2) uc2 = uc1;
       }
         
      init_u = uc2;
*/

   been_here ++;


   if (been_here <6000)  {
      E->advection.timestep = E->data.delta_S;
      }
   else if (been_here <5000)  {
      E->advection.timestep = 2.*E->data.delta_S;
      }
   else  {
      E->advection.timestep = 1.0;
      E->advection.timestep = 5.*E->data.delta_S;
      }

    return; 
  }
