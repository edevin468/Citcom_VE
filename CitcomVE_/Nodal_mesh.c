/* Functions relating to the building and use of mesh locations ... */


#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

/* =================================================
   Standard node positions including mesh refinement 

   =================================================  */

void node_locations(E)
     struct All_variables *E;
{ 
  int lev,i,j,k,ijk[4],ii,d,node,ixbasin;
  double *XX[4],*XG[4],dx[4],*dxx[2],dx1,dx2,dc,rc,dr1,dr2;
  int nox,noz,noy,fn,step,ncr,ndc,ll,mm;
  double modified_plgndr_a(),con,t1;
  void inject_scalar_d();
  static int iii=0;

  const int dims = E->mesh.nsd;

  if (iii==0) {
         E->viscosity.z410 = E->viscosity.z410/E->sphere.radius;
         E->viscosity.zlm = E->viscosity.zlm/E->sphere.radius;
         iii++;
       }

     for(d=1;d<=E->mesh.nsd;d++) {
       XX[d] = (double *)malloc((2+E->mesh.nnx[d])*sizeof(double)); 
       XG[d] = (double *)malloc ((2+1)*sizeof(double)); 
       }

     for(d=1;d<=E->mesh.nsd;d++) {     /* for even space */
        dx[d]                 = E->mesh.layer[d]/(E->mesh.nnx[d]-1);
        XX[d][1]              = 0.0;
        XX[d][E->mesh.nnx[d]] = E->mesh.layer[d];
        for(i=2;i<E->mesh.nnx[d];i++)
              XX[d][i] = XX[d][i-1]+dx[d];
        }

/*
      dx1 = 0.0030;
      fn = E->mesh.nnx[1];

      dc = 2.0*(XX[1][fn]-(fn-1)*dx1)/((fn-1)*(fn-2));

      for(i=2;i<fn;i++)
            XX[1][i] = XX[1][i-1]+dx1+(i-2)*dc;
*/

    lev = E->mesh.levmax;
    nox=E->mesh.NOX[lev];
    noy=E->mesh.NOY[lev];
    noz=E->mesh.NOZ[lev];
    

/* comment it out for uniform z mesh */

  if (E->viscosity.nlm==0) {
     rc = E->viscosity.z410/(E->viscosity.n410-1);
     dc = (1.0-E->viscosity.z410)/(noz-E->viscosity.n410);

fprintf(stderr,"ok1, %g %g %g\n",E->viscosity.z410,rc,dc);

     XX[2][1] = 0.0;
     for(i=2;i<=E->mesh.nnx[2];i++)   {
        if (i<=E->viscosity.n410)
           XX[2][i] = XX[2][i-1]+rc;
        else if (i>E->viscosity.n410)
           XX[2][i] = XX[2][i-1]+dc;
        }
     }
  else {
     dr1 = E->viscosity.zlm/(E->viscosity.nlm-1);
     rc = (E->viscosity.z410-E->viscosity.zlm)/(E->viscosity.n410-E->viscosity.nlm);
     dc = (1.0-E->viscosity.z410)/(noz-E->viscosity.n410);

fprintf(stderr,"ok2, %g %g %g %g %g\n",E->viscosity.z410,E->viscosity.zlm,dr1,rc,dc);

     XX[2][1] = 0.0;
     for(i=2;i<=E->mesh.nnx[2];i++)   {
        if (i<=E->viscosity.nlm)
           XX[2][i] = XX[2][i-1]+dr1;
        else if (i<=E->viscosity.n410)
           XX[2][i] = XX[2][i-1]+rc;
        else if (i>E->viscosity.n410)
           XX[2][i] = XX[2][i-1]+dc;
        }
     }
/**/

     for(i=1;i<=nox;i++)   
     for(j=1;j<=noz;j++)   {
           node=j+(i-1)*noz;
           E->XX[lev][1][node] = XX[1][i]; 
           E->XX[lev][2][node] = XX[2][j]; 
           }

    if (E->control.NMULTIGRID||E->control.EMULTIGRID)
      for (lev=E->mesh.levmax;lev>E->mesh.levmin;lev--)  {
        inject_scalar_d(E,lev,E->XX[lev][1],E->XX[lev-1][1]);
        inject_scalar_d(E,lev,E->XX[lev][2],E->XX[lev-1][2]);
        }

     for(d=1;d<=E->mesh.nsd;d++) {
       free((void *)XX[d]);
       free((void *)XG[d]);
       }

  if (E->control.verbose)    
    for (lev=E->mesh.levmin;lev<=E->mesh.levmax;lev++) {
      fprintf(E->fp,"output_coordinates %d\n",lev);
      if (dims==2)  {
         for (i=1;i<=E->mesh.NNO[lev];i++)
             fprintf(E->fp,"%d %g %g\n",i,E->XX[lev][1][i],E->XX[lev][2][i]);
         }
      else if (dims==3)  {
         for (i=1;i<=E->mesh.NNO[lev];i++)
             fprintf(E->fp,"%d %g %g %g\n",i,E->XX[lev][1][i],E->XX[lev][2][i],E->XX[lev][3][i]);
         }
      }

return;  }




void flogical_mesh_to_real(E,data,level)
     struct All_variables *E;
     float *data;
     int level;

{ int i,j,n1,n2;

  return;
}

void dp_to_nodes(E,P,PN,lev)
     struct All_variables *E;
     double *P;
     float *PN;
     int lev;

{ int e,element,node,j;

  for(node=1;node<=E->mesh.NNO[lev];node++)
    PN[node] =  0.0;
	  
  for(element=1;element<=E->mesh.NEL[lev];element++) {

      for(j=1;j<=enodes[E->mesh.nsd];j++)  {
     	  node = E->IEN[lev][element].node[j];
    	  PN[node] += P[element] * E->TW[lev][node] ; 
    	  }

      } 

     return; }

void p_to_nodes(E,P,PN,lev)
     struct All_variables *E;
     float *P,*PN;
     int lev;

{ int e,element,node,j;

  for(node=1;node<=E->mesh.NNO[lev];node++)
    PN[node] =  0.0;


  for(element=1;element<=E->mesh.NEL[lev];element++) {
      for(j=1;j<=enodes[E->mesh.nsd];j++)  {
     	  node = E->IEN[lev][element].node[j];
    	  PN[node] += P[element] * E->TWW[lev][element].node[j] ; 
    	  }
      } 

  for(node=1;node<=E->mesh.NNO[lev];node++)
      PN[node] *= E->MASS[lev][node];

     return; }


void p_to_centres(E,PN,P,lev)
     struct All_variables *E;
     float *PN,*P;
     int lev;

{  int p,element,node,j;
   double weight;

   for(p=1;p<=E->mesh.NEL[lev];p++)
     P[p] = 0.0;

   weight=1.0/((double)enodes[E->mesh.nsd]) ;
   
   for(p=1;p<=E->mesh.NEL[lev];p++)
     for(j=1;j<=enodes[E->mesh.nsd];j++)
       P[p] +=  PN[E->IEN[lev][p].node[j]] * weight;

   return;  
   }


void v_to_intpts(E,VN,VE,lev)
  struct All_variables *E;
  float *VN,*VE;
  int lev;
  {

   int e,i,j,k;
   const int nsd=E->mesh.nsd;
   const int vpts=vpoints[nsd];
   const int ends=enodes[nsd];

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)                 {
        VE[(e-1)*vpts + i] = 0.0;
        for(j=1;j<=ends;j++)
          VE[(e-1)*vpts + i] += VN[E->IEN[lev][e].node[j]] *  E->N.vpt[GNVINDEX(j,i)];
        }

   return;
  }

void v_to_nodes(E,VE,VN,lev)
   struct All_variables *E;
   float *VE,*VN;
   int lev;
   {
    int e,i,j,k,n;
    const int nsd=E->mesh.nsd;
    const int vpts=vpoints[nsd];
    const int ends=enodes[nsd];

    for(i=1;i<=E->mesh.NNO[lev];i++)
	    VN[i] = 0.0;

    for(e=1;e<=E->mesh.NEL[lev];e++)
      for(j=1;j<=ends;j++) {
        n = E->IEN[lev][e].node[j];
        for(i=1;i<=vpts;i++)
          VN[n] += E->N.vpt[GNVINDEX(j,i)] * E->TW[lev][n] * VE[(e-1)*vpts + i];
        }

    return;
    }

void visc_to_intpts(E,VN,VE,lev)
   struct All_variables *E;
   float *VN,*VE;
   int lev;
   {

   int e,i,j,k;
   const int nsd=E->mesh.nsd;
   const int vpts=vpoints[nsd];
   const int ends=enodes[nsd];

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++) {
        VE[(e-1)*vpts + i] = 0.0;
	for(j=1;j<=ends;j++)
          VE[(e-1)*vpts + i] += log(VN[E->IEN[lev][e].node[j]]) *  E->N.vpt[GNVINDEX(j,i)];
        VE[(e-1)*vpts + i] = exp(VE[(e-1)*vpts + i]);
        }

  }


void visc_to_nodes(E,VE,VN,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {
  int e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

  for(i=1;i<=E->mesh.NNO[lev];i++)
    VN[i] = 0.0;

  for(e=1;e<=E->mesh.NEL[lev];e++)
    for(j=1;j<=ends;j++) {
      n = E->IEN[lev][e].node[j];
      temp_visc=0.0;
      for(i=1;i<=vpts;i++)
	temp_visc += E->TW[lev][n] * log(E->N.vpt[GNVINDEX(j,i)] * VE[(e-1)*vpts + i]);
      VN[n] += exp(temp_visc);
      }
   return;
}

void visc_from_gint_to_nodes(E,VE,VN,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {
  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

   for(i=1;i<=E->mesh.NNO[lev];i++)
     VN[i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(j=1;j<=ends;j++)                {
       n = E->IEN[lev][e].node[j];
       temp_visc=0.0;
       for(i=1;i<=vpts;i++)
         temp_visc += E->TW[lev][n] * E->N.vpt[GNVINDEX(j,i)] * VE[(e-1)*vpts + i];
       VN[n] += temp_visc;
       }

   return;
}


 void visc_from_nodes_to_gint(E,VN,VE,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {

  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;
   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)
       VE[(e-1)*vpts+i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)      {
       temp_visc=0.0;
       for(j=1;j<=ends;j++)
         temp_visc += E->N.vpt[GNVINDEX(j,i)]*VN[E->IEN[lev][e].node[j]];

       VE[(e-1)*vpts+i] = temp_visc;
       }

   return;
   }

void visc_from_gint_to_ele(E,VE,VN,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {
  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

   for(i=1;i<=E->mesh.NEL[lev];i++)
     VN[i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)   {
     temp_visc=0.0;
     for(i=1;i<=vpts;i++)
        temp_visc += VE[(e-1)*vpts + i];
     temp_visc = temp_visc/vpts;

     VN[e] = temp_visc;
    }

   return;
}


 void visc_from_ele_to_gint(E,VN,VE,lev)
  struct All_variables *E;
  float *VE,*VN;
  int lev;
  {

  int m,e,i,j,k,n;
  const int nsd=E->mesh.nsd;
  const int vpts=vpoints[nsd];
  const int ends=enodes[nsd];
  double temp_visc;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)
       VE[(e-1)*vpts+i] = 0.0;

   for(e=1;e<=E->mesh.NEL[lev];e++)
     for(i=1;i<=vpts;i++)      {

       VE[(e-1)*vpts+i] = VN[e];
       }

   return;
 }

