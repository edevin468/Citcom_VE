#include "global_defs.h"
#include "element_definitions.h"
#include <math.h>
#include <stdio.h>

/* =================================================== */
/* =================================================== */
void set_free_surfaces(E)
struct All_variables *E;
{

int lev,i,j,node1,node2;
 float total, length, stress_scale;
 double modified_plgndr_a(),con,t1;
 char output_file[255];
 FILE *fp1;
 char input_s[100];
  
             /* mark the vertical d.o.f. for nodes on density boundaries */

 for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)
   for (i=1;i<=E->mesh.NOX[lev];i++)    {
     node1=i*E->mesh.NOZ[lev];
     node2=(i-1)*E->mesh.NOZ[lev]+1;
     E->EQN[lev][E->ID[lev][node1].doff[2]] = E->EQN[lev][E->ID[lev][node1].doff[2]] | RESTORE;
     E->EQN[lev][E->ID[lev][node2].doff[2]] = E->EQN[lev][E->ID[lev][node2].doff[2]] | RESTORE;
     }

 if (E->control.pure_visc==0)   {
   E->data.botm_scaling = (E->data.density_below-E->data.density)*E->data.grav_acc
              *E->sphere.radius/E->data.shear_mod;
   E->data.surf_scaling = (E->data.density-E->data.density_above)*E->data.grav_acc*E->sphere.radius
              /E->data.shear_mod;
   }
 else   {         /* pure viscous */
   stress_scale = E->data.ref_viscosity*E->data.therm_diff/(E->sphere.radius*E->sphere.radius);
   E->data.botm_scaling = (E->data.density_below-E->data.density)*E->data.grav_acc
              *E->sphere.radius/stress_scale;
   E->data.surf_scaling = E->data.density*E->data.grav_acc*E->sphere.radius
              /stress_scale;
   }

       fp1=fopen(E->control.data_file1,"r");
       for (i=0;i<=E->load.stage;i++)  {
         fgets(input_s,100,fp1);
         sscanf(input_s,"%d %g %lf",&j,&E->load.time[i],&E->load.mag[i]);
         }
       fclose(fp1);
        

return;
}


/* ======================================================= 
  before assembling a force vector, determine topography
  at each density interface and convert it into stress 
  drho*g*h
 =======================================================  */

 void get_boundary_forces(E) 
 struct All_variables *E;
 {
 int ll,mm,i,j;
 double total, length;
 double modified_plgndr_a(),con,t1;
 void remove_average();


 for (j=1;j<=E->mesh.nox;j++)    {
     i = (j-1)*E->mesh.noz + 1;
     E->Xsurf[1][j] = E->XX[E->mesh.levmax][1][i];
     E->Xsurf[2][j] = E->XX[E->mesh.levmax][2][i];
     }

  remove_average(E,E->Xsurf[1],E->Xsurf[2]);

  for (j=1;j<=E->mesh.nox;j++)    {
     E->slice.load[3][j] = E->slice.load[2][j] + E->Xsurf[2][j]*E->data.botm_scaling;
     }

/*
  for (j=1;j<=E->mesh.nox;j++)    {
     i = j*E->mesh.noz;
     E->Xsurf[1][j] = E->XX[E->mesh.levmax][1][i];
     E->Xsurf[2][j] = E->XX[E->mesh.levmax][2][i];
     }
  remove_average(E,E->Xsurf[1],E->Xsurf[2]);
*/

  for (j=1;j<=E->mesh.nox;j++)    {
     E->slice.load[1][j] = E->slice.load[0][j] + E->slice.surf[1][j]*E->data.surf_scaling;
     }

  return;
  }

/* ======================================================= 
   for each element el on boundaries, determine contributions
   from topographic deflection
 =======================================================  */
void  construct_B_R(E)
  struct All_variables *E;
 {

  void get_global_1d_shape_fn();

  double con,slope;
  double x[3][3];
  int lev, elz, el, i, k, n1, n2, p, p2;
  static int been_here=0;

  struct Shape_function1 GM;
  struct Shape_function1_dA dGammax;

if ( been_here==0 && E->control.pure_visc==0)      {
 for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  { 
  elz = E->mesh.ELZ[lev];
  for (el=1;el<=E->mesh.NEL[lev];el++)   
    if (el%elz==0 || (el-1)%elz==0)   {

      get_global_1d_shape_fn(E,lev,el,&GM,&dGammax,0);

      if (el%elz==0)   {
        p = 2;     /*      top  */
        n1 = E->IEN[lev][el].node[4];
        n2 = E->IEN[lev][el].node[3];
        con  =E->data.surf_scaling;
        }
      else   {
        p = 1;     /*      bottom */
        n1 = E->IEN[lev][el].node[1];
        n2 = E->IEN[lev][el].node[2];
        con  =E->data.botm_scaling;
        }

      if (E->control.AXI)
        for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
          x[1][i] = E->M.vpt[GMVINDEX(1,i)]*E->XX[lev][1][n1] 
                  + E->M.vpt[GMVINDEX(2,i)]*E->XX[lev][1][n2];
          x[2][i] = E->M.vpt[GMVINDEX(1,i)]*E->XX[lev][2][n1] 
                  + E->M.vpt[GMVINDEX(2,i)]*E->XX[lev][2][n2];
        }
      else if (E->control.CART2D)
        for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
          x[1][i] = 1.0;
          x[2][i] = 1.0;
        }

      for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
        E->B_R[lev][E->ID[lev][n1].doff[2]] -= 
                     x[1][i]*E->M.vpt[GMVINDEX(1,i)]
                      * con * dGammax.vpt[GMVGAMMA(p,i)];
        E->B_R[lev][E->ID[lev][n2].doff[2]] -= 
                     x[1][i]*E->M.vpt[GMVINDEX(2,i)]
                      * con * dGammax.vpt[GMVGAMMA(p,i)];
        }

   }
 }     /* end for lev */
  been_here++; 
}     
else
 for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  
    for (i=0;i<=E->mesh.NEQ[lev];i++)
        E->B_R[lev][i] = 0.0;


 return;
 }

/* ======================================================= 
   for each element el on boundaries, determine contributions
   from topographic deflection to elt_f
 =======================================================  */

 void load_boundary_forces(E,el,elt_f)
  struct All_variables *E;
  int el;
  double elt_f[24];
  {

  void get_global_1d_shape_fn();

  double x[3][3],force[9],force_at_gs[9];
  int e, i, k, n1, n2, p, p2;

  struct Shape_function1 GM;
  struct Shape_function1_dA dGammax;

  get_global_1d_shape_fn(E,E->mesh.levmax,el,&GM,&dGammax,0);

  if (el%E->mesh.elz==0)   {
     e = el/E->mesh.elz;
     for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
       force[k] = E->slice.load[1][ E->sien[e].node[k] ];
     p = 2;     /*      top  */
     n1 = 4;
     n2 = 3;
     }
  else   {
     e = (el-1)/E->mesh.elz+1;
     for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
       force[k] = E->slice.load[3][ E->sien[e].node[k] ];
     p = 1;     /*      bottom */
     n1 = 1;
     n2 = 2;
     }

  if (E->control.AXI)
    for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
      x[1][i] = E->M.vpt[GMVINDEX(1,i)]*E->X[1][E->ien[el].node[n1]] 
              + E->M.vpt[GMVINDEX(2,i)]*E->X[1][E->ien[el].node[n2]];
      x[2][i] = E->X[2][E->ien[el].node[n2]];
      }
  else if (E->control.CART2D)
    for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
          x[1][i] = 1.0;
          x[2][i] = 1.0;
        }

  for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
      force_at_gs[i] = 0.0;
      for(k=1;k<=onedvpoints[E->mesh.nsd];k++)
         force_at_gs[i] += force[k] * E->M.vpt[GMVINDEX(k,i)];
      }

  for (i=1;i<=onedvpoints[E->mesh.nsd];i++)   {
      elt_f[2*n1-1] += x[1][i]*E->M.vpt[GMVINDEX(1,i)]
           * force_at_gs[i] * dGammax.vpt[GMVGAMMA(p,i)];
      elt_f[2*n2-1] += x[1][i]*E->M.vpt[GMVINDEX(2,i)]
           * force_at_gs[i] * dGammax.vpt[GMVGAMMA(p,i)];
      }

  return;
  }

/* =================================================== */
/* =================================================== */
 void get_surface_time_loads(E) 
 struct All_variables *E;
 {
 int ll,mm,i,j;
 float total, length, stress_scale;
 double con,x1;
 void remove_average();
 static int index=0;

  x1 = E->monitor.elapsed_time-0.5*E->advection.timestep;
  if(x1<0.0) x1=0.0;

    con = 1e-14+E->convection.perturb_mag[0]*sin(2.0*M_PI*x1/E->data.T_sol0);

    for (j=1;j<=E->mesh.nox;j++)    {
      i = j*E->mesh.noz;
      x1 = E->Xsurf[1][j] = E->XX[E->mesh.levmax][1][i];
      E->Xsurf[2][j] = con*cos(x1*M_PI/E->convection.perturb_ll[0]);
     }

  remove_average(E,E->Xsurf[1],E->Xsurf[2]);

  for (j=1;j<=E->mesh.nox;j++)    {
     E->slice.load[0][j] = E->Xsurf[2][j]*E->data.surf_scaling;
     E->slice.load[2][j] = 0.0;
     }

  return;

  }
