#include "global_defs.h"
#include "element_definitions.h"
#include <math.h>


void pickup_dt(E)
 struct All_variables *E;
 {

 int i,j;
 float dt,dth,dr;
 void std1_timestep();
 void std2_timestep();

 if (E->control.pure_visc==0)    {
/*   E->advection.timestep = E->data.delta_S;    */
   std2_timestep(E);
   } 
 else   {

   std1_timestep(E);

fprintf(E->fp,"dt=%.5e\n",E->advection.timestep); 
fflush(E->fp);

   }
 

 return;
 }


/*   ===================================================
  calculate the new positions of each nodes by adding 
  the displacements onto the old coordinates for each 
  nodes. And inject the new positions onto different 
  levels for multigrid
     =================================================== */

void deform_grid(E)
 struct All_variables *E;

  {
  
  int i,j,lev;
  void inject_scalar_d();
  void mass_matrix();


 if (E->control.NMULTIGRID||E->control.EMULTIGRID)   {
    for (j=1;j<=E->mesh.nsd;j++)   {

      if (E->control.pure_visc==0) 
        for (i=1;i<=E->mesh.nno;i++)  {
          E->X0[j][i] += E->U[E->id[i].doff[j]];
          E->XX[E->mesh.levmax][j][i] += E->U[E->id[i].doff[j]];
          }
      else
        for (i=1;i<=E->mesh.nno;i++)  {
          E->X0[j][i] += E->advection.timestep*E->U[E->id[i].doff[j]];
          E->XX[E->mesh.levmax][j][i] += E->advection.timestep*E->U[E->id[i].doff[j]];
          }

      for (lev=E->mesh.levmax;lev>E->mesh.levmin;lev--)
        inject_scalar_d(E,lev,E->XX[lev][j],E->XX[lev-1][j]);
    
      }
      for (i=1;i<=E->mesh.nsf;i++)  {
        E->slice.surf[1][i] +=E->U[E->id[i*E->mesh.noz].doff[2]];
        }
    }
  else  
    for (j=1;j<=E->mesh.nsd;j++)    {
      if (E->control.pure_visc==0) 
        for (i=1;i<=E->mesh.nno;i++)  {
          E->X0[j][i] += E->U[E->id[i].doff[j]];
          E->X[j][i] += E->U[E->id[i].doff[j]];
          }
      else
        for (i=1;i<=E->mesh.nno;i++)  {
          E->X0[j][i] += E->advection.timestep*E->U[E->id[i].doff[j]];
          E->X[j][i] += E->advection.timestep*E->U[E->id[i].doff[j]];
          }
    }

  mass_matrix(E);

  return;
  }

/*   ===================================================
    get the surface topography 
     =================================================== */

void print_surf_topo(E,ii)
 struct All_variables *E;
 int ii;
  {
  
  char outfile[255];
  FILE *fp;
  static FILE *fp1,*fp2;
  int node,i,j;
  static int been_here=0;
  static int iii=0;
  float error, tot, amp, tau, tau1,tau2,tau3,ana1,ana3,ana2, total,seconds_in_a_year,temp[601];
  int ll,mm;
  double modified_plgndr_a(),con,t1;

  if (been_here!=0 && E->mesh.nmx!=iii)   {
   fclose(fp1);
   fclose(fp2);
   }

  if (been_here==0 || E->mesh.nmx!=iii)   {
    sprintf(outfile,"%s/c%02d.topo_time",E->control.data_file,E->mesh.nmx);
    fp1 = fopen(outfile,"w");
    sprintf(outfile,"%s/c%02d.time",E->control.data_file,E->mesh.nmx);
    fp2 = fopen(outfile,"w");
    been_here++;
    }
 
  iii = E->mesh.nmx;


 
  seconds_in_a_year = 24*365*3600.0;
  tau = 1.0e21/(E->data.shear_mod*seconds_in_a_year);

    total = 0.0;
    tot = 0.0;
    for(i=1;i<=E->mesh.elx;i++)              {
        node = i*E->mesh.noz;
        total = total + (E->X[2][node]+E->X[2][node+E->mesh.noz])
                       *(E->X[1][node+E->mesh.noz]-E->X[1][node])*0.5;
        }
    total = total/(E->X[1][E->mesh.nno]-E->X[1][1]);


    for(i=1;i<=E->mesh.nox;i++)   {
       node = i*E->mesh.noz;
       temp[i] = E->X[2][node] - total;
       }


  if ( ((ii % E->control.record_every) == 0))    {

    sprintf(outfile,"%s/c%02d.topo_s_%05d.dat",E->control.data_file,E->mesh.nmx,ii);
    fp = fopen(outfile,"w");
    fprintf(fp2,"%05d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %5e\n",ii,E->monitor.elapsed_time,E->viscosity.dissipation_total,E->viscosity.elastic_total,E->viscosity.dissipation_total+E->viscosity.elastic_total,E->viscosity.work_done_total,E->viscosity.dissipation,E->viscosity.elastic,E->viscosity.work_done);

//    fprintf(fp2,"%05d %.5e %.5e %.5e %.5e %.5e\n",ii,E->monitor.elapsed_time,
//   E->viscosity.dissipation_total-E->viscosity.dissipation_total_old,
//   E->viscosity.elastic_total-E->viscosity.elastic,
//   E->viscosity.work_done,
//   fabs(E->viscosity.work_done)-fabs(E->viscosity.dissipation_total-E->viscosity.dissipation_total_old)-(E->viscosity.elastic_total-E->viscosity.elastic));

    total = 0;
    j = 0;
    for (i=1;i<=E->mesh.nno;i++)  
      if (i%E->mesh.noz==0)     {
        j++;
        E->Xsurf[1][j] = E->X[1][i];
        E->Xsurf[2][j] = E->X[2][i];
        fprintf(fp,"%.5e %.5e %.5e %.5e %.5e\n",E->X[1][i],E->slice.surf[1][j],E->U[E->id[i].doff[2]],E->U[E->id[i].doff[2]]/E->advection.timestep,E->slice.load[0][j]/E->data.surf_scaling);

        }
    fclose(fp);


    fflush(fp2);

    }

  return;
  }

/*========================================================= */
/*========================================================= */
void get_temperature_half_space(E)
 struct All_variables *E;
 {

 float timea;

 timea = E->monitor.elapsed_time;

 return;
 }
    

/*   ===================================================
  update stresses at each nodes by adding 
     =================================================== */

void update_stress_strain(E,ii)
 struct All_variables *E;
 int ii;
 {

    void get_global_shape_fn();
    int i,j,k,e,node;
   
    static double SXX[5],SXZ[5],SZZ[5],SYY[5];
    static double EXX[5],EXZ[5],EZZ[5],EYY[5];
    static int been_here=0;
    double VZ[9],VY[9],VX[9],Szz,Sxx,Sxz,Syy,Exx,Exz,Ezz,Eyy;
    double el_volume,visc[9],Visc,a,b,xk[3][5],Sxx1,Szz1,Sxz1,Syy1;
    double work,elastic,dissipation,dissipation1,stress[9],strain[9],stress1[9];
    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;

    void remove_average();

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int vpts=vpoints[dims];
    const int ppts=ppoints[dims];
    const int ends=enodes[dims];
    const int nno=E->mesh.nno;
    const int lev=E->mesh.levmax;


    a=1e10;
    E->viscosity.elastic = 0;
    dissipation = 0.0;

    for(e=1;e<=E->mesh.nel;e++)  {

      for(i=1;i<=vpts;i++)  {
      SZZ[i] = 0.0;
      SXX[i] = 0.0;
      SXZ[i] = 0.0;
      SYY[i] = 0.0;
      EZZ[i] = 0.0;
      EXX[i] = 0.0;
      EXZ[i] = 0.0;
      EYY[i] = 0.0;
      visc[i] = E->EVI[E->mesh.levmax][(e-1)*vpts+i];
      }

      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,xk,0,E->mesh.levmax);

      for(j=1;j<=ends;j++) {
         VX[j] = E->U[E->id[E->ien[e].node[j]].doff[1]];
         VZ[j] = E->U[E->id[E->ien[e].node[j]].doff[2]];
         }

      Szz = Sxx = Sxz = Syy = 0.0;
      E->S2zz[e] = E->S2xx[e] = E->S2xz[e] = E->S2yy[e] = 0.0;
      for(i=1;i<=vpts;i++)  {
        for(j=1;j<=ends;j++)  {
            EZZ[i] += 2.0 * VZ[j]* GNx.vpt[GNVXINDEX(1,j,i)];
            EXX[i] += 2.0 * VX[j]* GNx.vpt[GNVXINDEX(0,j,i)];
            EXZ[i] += VX[j]* GNx.vpt[GNVXINDEX(1,j,i)] + VZ[j]* GNx.vpt[GNVXINDEX(0,j,i)];
            }
        if (E->control.AXI)   {
          for(j=1;j<=ends;j++)  
            EYY[i] += 2.0 * (VX[j]* xk[1][i])* E->N.vpt[GNVINDEX(j,i)];
          E->Syy[(e-1)*vpts+i] = visc[i]*EYY[i]+E->Syy[(e-1)*vpts+i]*E->Maxwelltime[e];
          E->S2yy[e] += E->Syy[(e-1)*vpts+i];
          SYY[(e-1)*vpts+i] += E->Syy[(e-1)*vpts+i];
          }

// old stress
        SXX[i] = E->Sxx[(e-1)*vpts+i];
        SXZ[i] = E->Sxz[(e-1)*vpts+i];
        SZZ[i] = E->Szz[(e-1)*vpts+i];
        if (E->control.AXI)   {
          SYY[i] = E->Syy[(e-1)*vpts+i];
          }
       
// new stress
        E->Szz[(e-1)*vpts+i] = visc[i]*EZZ[i]+E->Szz[(e-1)*vpts+i]*E->Maxwelltime[e];
        E->Sxx[(e-1)*vpts+i] = visc[i]*EXX[i]+E->Sxx[(e-1)*vpts+i]*E->Maxwelltime[e];
        E->Sxz[(e-1)*vpts+i] = visc[i]*EXZ[i]+E->Sxz[(e-1)*vpts+i]*E->Maxwelltime[e];

        E->S2zz[e] += E->Szz[(e-1)*vpts+i];
        E->S2xx[e] += E->Sxx[(e-1)*vpts+i];
        E->S2xz[e] += E->Sxz[(e-1)*vpts+i];

        stress[i] = E->Sxx[(e-1)*vpts+i]*E->Sxx[(e-1)*vpts+i]
                  + E->Szz[(e-1)*vpts+i]*E->Szz[(e-1)*vpts+i]
                  + E->Sxz[(e-1)*vpts+i]*E->Sxz[(e-1)*vpts+i]*2.0;
        stress1[i] = E->Sxx[(e-1)*vpts+i]*SXX[i]
                   + E->Szz[(e-1)*vpts+i]*SZZ[i]
                   + E->Sxz[(e-1)*vpts+i]*SXZ[i]*2.0;
//        stress1[i] = SXX[i]*SXX[i] + SZZ[i]*SZZ[i] + SXZ[i]*SXZ[i]*2.0;

        strain[i] = EXX[i]*EXX[i] + EZZ[i]*EZZ[i] + EXZ[i]*EXZ[i]*2.0;

        if (E->control.AXI)  { 
          stress[i] += SYY[i]*SYY[i];
          strain[i] += EYY[i]*EYY[i];
          }

        stress[i] = (0.5*stress[i]);
        stress1[i] = (0.5*stress1[i]);
        strain[i] = sqrt(0.5*strain[i]);
        }
    
      E->S2zz[e] = E->S2zz[e]/vpts;
      E->S2xx[e] = E->S2xx[e]/vpts;
      E->S2xz[e] = E->S2xz[e]/vpts;
      E->S2yy[e] = E->S2yy[e]/vpts;

      E->PO[e] = E->P[e]; 

      E->Ezz[e] = E->S2zz[e]*E->S2zz[e] + E->S2xx[e]*E->S2xx[e]
                + E->S2xz[e]*E->S2xz[e]*2.0;
     if (E->control.AXI)   
        E->Ezz[e] += E->S2yy[e]*E->S2yy[e];

     E->Ezz[e] = sqrt(0.5*E->Ezz[e]);

     dissipation1 = 0.0;
     elastic = 0.0;
     el_volume = 0.0;
     for(j=1;j<=vpts;j++)  {
        el_volume += -dOmega.vpt[j]*g_point[j].weight[dims-1];
        dissipation1 += -stress[j]/E->viscosity.element[(e-1)*vpts+j]*dOmega.vpt[j]*g_point[j].weight[dims-1];	//compute power
        elastic += -(stress[j]-stress1[j])*dOmega.vpt[j]*g_point[j].weight[dims-1]/E->advection.timestep;	//compute power

//	dissipation1 += -stress[j]/E->viscosity.element[(e-1)*vpts+j]*dOmega.vpt[j]*g_point[j].weight[dims-1]*E->advection.timestep; 
//        elastic += -(stress[j]-stress1[j])*dOmega.vpt[j]*g_point[j].weight[dims-1]; 

//   the following gives almost the same results 
//        dissipation += -E->Ezz[e]*E->Ezz[e]/E->viscosity.element[(e-1)*vpts+j]*dOmega.vpt[j]*g_point[j].weight[dims-1]; 
        }

     E->viscosity.dissipation_e[e] = dissipation1/el_volume;
     dissipation += (dissipation1);

     E->viscosity.elastic += (elastic);
     }

// time integral of dissipation for this delta_t using dissipation from 
// the current step and the previous step
   if (ii>0) {
      E->viscosity.dissipation_total_old = E->viscosity.dissipation_total;
      E->viscosity.dissipation_total += (fabs(dissipation)+E->viscosity.dissipation)*0.5;
      }
   E->viscosity.dissipation = fabs(dissipation);

   E->viscosity.elastic_total += E->viscosity.elastic;

   for (j=1;j<=E->mesh.nox;j++)    {
     i = j*E->mesh.noz;
     E->Xsurf[1][j] = E->XX[E->mesh.levmax][1][i];
     E->Xsurf[2][j] = E->U[E->id[i].doff[2]];
     }

  remove_average(E,E->Xsurf[1],E->Xsurf[2]);

   work = 0;
   for (j=2;j<=E->mesh.nox;j++)   {
      i = j*E->mesh.noz;

      work+=(E->slice.load[1][j]  *E->Xsurf[2][j]+
             E->slice.load[1][j-1]*E->Xsurf[2][j-1])
             *0.5*(E->X[1][i]-E->X[1][i-E->mesh.noz]);

//      work+=(E->Xsurf[2][j]*E->Xsurf[2][j]+
//             E->Xsurf[2][j-1]*E->Xsurf[2][j-1])
 //            *0.5*E->data.surf_scaling*(E->X[1][i]-E->X[1][i-E->mesh.noz]);
   }
work = work/E->advection.timestep;	

// for CMB deformation
   for (j=1;j<=E->mesh.nox;j++)    {
     i = (j-1)*E->mesh.noz+1;
     E->Xsurf[1][j] = E->XX[E->mesh.levmax][1][i];
     E->Xsurf[2][j] = E->U[E->id[i].doff[2]];
     }

  remove_average(E,E->Xsurf[1],E->Xsurf[2]);

   for (j=2;j<=E->mesh.nox;j++)   {
      i = (j-1)*E->mesh.noz+1;

//    work = work/E->advection.timestep;

      work+=(E->slice.load[3][j]  *E->Xsurf[2][j]+
             E->slice.load[3][j-1]*E->Xsurf[2][j-1])
             *0.5*(E->X[1][i]-E->X[1][i-E->mesh.noz]);
   }




   E->viscosity.work_done_total += work;

   E->viscosity.work_done = work;


 return;
 }

/*========================================================= */
/*========================================================= */

  void add_viscoelasticity(E, evisc)
  struct All_variables *E;
  float *evisc;
  {

  float visc, alpha;
  const int vpts = vpoints[E->mesh.nsd];
  int i,j;


  for (i=1;i<=E->mesh.nel;i++)   {
    visc = 0.0;
    for(j=1;j<=vpts;j++) 
       visc += evisc[(i-1)*vpts+j];
    visc = visc/vpts;

    alpha = visc/E->Maxwelltime[i]; 
    E->Maxwelltime[i] = (alpha-E->advection.timestep/2.)/(alpha+E->advection.timestep/2.);
    for(j=1;j<=vpts;j++) 
      evisc[(i-1)*vpts+j]=evisc[(i-1)*vpts+j]/(E->advection.timestep/2.+alpha);
    }

                           /* Maxwelltime is now (alpha-dt)/(alpha+dt/2) */
                           /* evisc is now evisc/(dt/2+alpha) */

  return;
  }
