#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"

void phase_change(E,B,B_b)
  struct All_variables *E;
  float *B,*B_b;
{
  int i,j,k,n,ns;
  float e_pressure,pt5,one;

  pt5 = 0.5; one=1.0;

  for(i=1;i<=E->mesh.nno;i++)  {
    e_pressure = E->X[2][i]-E->viscosity.zlm-
            E->control.clapeyron670*(E->T[i]-E->control.transT670);
    B[i] = pt5*(one+tanh(E->control.width670*e_pressure));
    }

  ns = 0;
  for (k=1;k<=E->mesh.noy;k++)
    for (j=1;j<=E->mesh.nox;j++)  {
      ns = ns + 1;
      B_b[ns]=0.0;
      for (i=1;i<E->mesh.noz;i++)   {
        n = (k-1)*E->mesh.noz*E->mesh.nox + (j-1)*E->mesh.noz + i;
        if (B[n]<=pt5&&B[n+1]>=pt5)
          B_b[ns]=(E->X[2][n+1]-E->X[2][n])*(pt5-B[n])/(B[n+1]-B[n])+E->X[2][n];
        }
      }


  return;
  }
