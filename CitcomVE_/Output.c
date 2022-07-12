/* Routine to process the output of the finite element cycles 
   and to turn them into a coherent suite  files  */


#include <fcntl.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>             /* for "system" command */
#ifndef __sunos__               /* string manipulations */
#include <strings.h>
#else
#include <string.h>
#endif

#include "element_definitions.h"
#include "global_defs.h"

void output_velo_related(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  int el,els,i,j,k,ii,m,node,fd;
  int nox,noz,noy,nfx,nfz,nfy1,nfy2,size1,size2;
  char output_file[255];
  static float *SV,*EV;
  static int been_here=0;

  void get_surface_velo ();
  void get_ele_visc ();
  void coordinates_dx(); 
  void frames_for_dx(); 

  const int nno = E->mesh.nno;

  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s/velo",E->control.data_file);
    E->filed[10]=fopen(output_file,"w");
    been_here++;

    }

   fprintf(E->filed[10],"%6d %6d %.5e\n",E->mesh.nno,E->advection.timesteps,E->monitor.elapsed_time);
    for (i=1;i<=E->mesh.nno;i++)
      fprintf(E->filed[10],"%6d %.3e %.3e %.5e %.5e %.5e %.5e\n",i,E->X[1][i],E->X[2][i],E->V[1][i],E->V[2][i],E->T[i],E->C[i]);
  fflush(E->filed[10]);

return;
}

void output_stress(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  int el,els,i,j,k,ii,m,node,fd;
  int nox,noz,noy,nfx,nfz,nfy1,nfy2,size1,size2;
  float visc;
  char output_file[255];
  static float *SV,*EV;
  static int been_here=0;
  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int vpts=vpoints[dims];
  FILE *fp;

  void get_surface_velo ();
  void get_ele_visc ();
  void coordinates_dx(); 
  void frames_for_dx(); 
  void p_to_nodes();

  const int nno = E->mesh.nno;

  if (been_here==0)  {
    SV = (float *) malloc ((E->mesh.nno+2)*sizeof(float));
    been_here ++;
    }


  p_to_nodes(E,E->Ezz,SV,E->mesh.levmax);

  p_to_nodes(E,E->viscosity.dissipation_e,E->viscosity.dissipation_n,E->mesh.levmax);

   sprintf(output_file,"%s/c%02d.stress.%d",E->control.data_file,E->mesh.nmx,file_number);
   fp=fopen(output_file,"w");

   fprintf(fp,"%6d %6d %.5e %.5e\n",E->mesh.nno,E->advection.timesteps,E->monitor.elapsed_time,E->viscosity.dissipation_total);
  for (i=1;i<=E->mesh.noz;i++)   {
      fprintf(fp,"%.4e %.4e %.4e\n",E->sphere.radius*(1-E->X[2][i])/1000,E->T[i]*E->data.ref_temperature+273,E->Vi[i]*E->data.ref_viscosity);
      }
   for (i=1;i<=E->mesh.nno;i++)   {
      fprintf(fp,"%.4e %.4e %.4e\n",SV[i],E->viscosity.dissipation_n[i],E->Vi[i]);
      }
   fclose(fp);

return;
}

/*      ----------------------------------- */
 void coordinates_dx(E)
  struct All_variables *E;
 {

 int i;

  E->ibm_dx.x1 = (float *) malloc((E->mesh.nno+1)*sizeof(float));
  E->ibm_dx.x2 = (float *) malloc((E->mesh.nno+1)*sizeof(float));

  E->ibm_dx.nox = E->mesh.nox;
  E->ibm_dx.noz = E->mesh.noz;

   for (i=1;i<=E->mesh.nno;i++)   {
      E->ibm_dx.x1[i] = E->X[2][i] * sin(E->X[1][i]);
      E->ibm_dx.x2[i] = E->X[2][i] * cos(E->X[1][i]);
      }

 return;
 }

/*      ----------------------------------- */
  void  frames_for_dx(E,file_number)
  struct All_variables *E;
  int file_number; 
  {

  int i;
  static int nframe=0;
  FILE *fp;
  char output_file[255];
  const float offset1 = -1.2;
  const float offset2 = 0.2;

   nframe ++;

   sprintf(output_file,"%s/mv.%03d.dx",E->control.data_file,nframe);
   fp=fopen(output_file,"w");
   for (i=1;i<=E->mesh.nno;i++)
     fprintf(fp,"%g %g %g\n",E->ibm_dx.x1[i]+offset1,E->ibm_dx.x2[i],E->T[i]);
   fclose(fp);

   sprintf(output_file,"%s/nv.%03d.dx",E->control.data_file,nframe);
   fp=fopen(output_file,"w");
   for (i=1;i<=E->mesh.nno;i++)
     fprintf(fp,"%g %g %g\n",E->ibm_dx.x1[i]+offset2,E->ibm_dx.x2[i],E->C[i]);
   fclose(fp);

   sprintf(output_file,"%s/mv.%03d.general",E->control.data_file,nframe);
   fp=fopen(output_file,"w");
   fprintf(fp,"file = mv.%03d.dx\n",nframe);
   fprintf(fp,"grid = %2d x %2d\n",E->ibm_dx.nox,E->ibm_dx.noz);
   fprintf(fp,"format = ascii\n");
   fprintf(fp,"interleaving = field\n");
   fprintf(fp,"majority = row\n");
   fprintf(fp,"field = locations, field0\n");
   fprintf(fp,"structure = 2-vector, scalar\n");
   fprintf(fp,"type = float, float\n");
   fprintf(fp,"\n");
   fprintf(fp,"end\n");
   fclose(fp);

   sprintf(output_file,"%s/nv.%03d.general",E->control.data_file,nframe);
   fp=fopen(output_file,"w");
   fprintf(fp,"file = nv.%03d.dx\n",nframe);
   fprintf(fp,"grid = %2d x %2d\n",E->ibm_dx.nox,E->ibm_dx.noz);
   fprintf(fp,"format = ascii\n");
   fprintf(fp,"interleaving = field\n");
   fprintf(fp,"majority = row\n");
   fprintf(fp,"field = locations, field0\n");
   fprintf(fp,"structure = 2-vector, scalar\n");
   fprintf(fp,"type = float, float\n");
   fprintf(fp,"\n");
   fprintf(fp,"end\n");
   fclose(fp);



   return;
   }

void output_velo_related_binary(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  int el,els,i,j,k,ii,m,node,fd;
  int nox,noz,noy,nfx,nfz,nfy1,nfy2,size1,size2;
  char output_file[255];
  static float *SV,*EV;
  static int been_here=0;

  void get_surface_velo ();
  void get_ele_visc ();
  const int nno = E->mesh.nno;
/*
  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s.velo",E->control.data_file);
    E->filed[10]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_t",E->control.data_file);
    E->filed[11]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_b",E->control.data_file);
    E->filed[12]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.visc",E->control.data_file);
    E->filed[13]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.fas670",E->control.data_file);
    E->filed[14]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.stress",E->control.data_file);
    E->filed[9]=open(output_file,O_RDWR | O_CREAT, 0644);
    }

  if (been_here==0)  {
    ii = E->mesh.nsf;
    SV = (float *) malloc ((2*ii+2)*sizeof(float));

    size2 = (E->mesh.nel+1)*sizeof(float);
    EV = (float *) malloc (size2);
    been_here++;
    }

  ii = E->mesh.nsf;
  size2 = 2*(ii+2)*sizeof(float);
    get_surface_velo (E,SV);
  write(E->filed[10],SV,size2);

  size2 = (E->mesh.nsf+1)*sizeof(float);
  write(E->filed[11],E->slice.tpg,size2);
  write(E->filed[12],E->slice.tpgb,size2);

  size2 = (E->mesh.nel+1)*sizeof(float);
    get_ele_visc (E,EV);
  write(E->filed[13],EV,size2);

  size2 = (E->mesh.nsf+1)*sizeof(float);
  write(E->filed[14],E->Fas670_b,size2);

  size2 = (2*E->mesh.nsf+1)*sizeof(float);
  write(E->filed[9],E->stress,size2);
*/
  return;
  }

/* ====================================================================== */

void output_temp(E,file_number)
  struct All_variables *E;
  int file_number; 
{
  int nno,i,j,fd;
  static int *temp1;
  static int been_here=0;
  static int size2,size1;
  char output_file[255];
/*
  if (been_here==0 && E->control.restart==0) {
    sprintf(output_file,"%s.temp",E->control.data_file);
    E->filed[5]=open(output_file,O_RDWR | O_CREAT, 0644);
    }

  if (been_here==0) {
    temp1 = (int *) malloc ((E->mesh.noy*6)*sizeof(int));

    sprintf(output_file,"%s.mesh",E->control.data_file);
    E->filed[1]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.x",E->control.data_file);
    E->filed[2]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.z",E->control.data_file);
    E->filed[3]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.y",E->control.data_file);
    E->filed[4]=open(output_file,O_RDWR | O_CREAT, 0644);

    size1 = (E->mesh.noy*6)*sizeof(int);
    size2= (E->mesh.nno+1)*sizeof(float);

    temp1[1] = E->mesh.nno;
    temp1[3] = size2;
    temp1[5] = E->mesh.nsf;
    temp1[6] = E->mesh.nel;

        write(E->filed[1],temp1,size1);
        write(E->filed[2],E->X[1],size2);
        write(E->filed[3],E->X[2],size2);
        write(E->filed[4],E->X[3],size2);

    close(E->filed[1]);
    close(E->filed[2]);
    close(E->filed[3]);
    close(E->filed[4]);

    been_here++;
    }


    write(E->filed[5],E->T,size2);
*/
  return; 
}


/* ====================================================================== */

void process_restart(E)
    struct All_variables *E;
{
    int fileid[20];
    int i,j,k,ii,size2;
    char output_file[255],in_file[255];
/*
    sprintf(output_file,"%s.temp",E->control.data_file);
    E->filed[5]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.velo",E->control.data_file);
    E->filed[10]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_t",E->control.data_file);
    E->filed[11]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.topo_b",E->control.data_file);
    E->filed[12]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.visc",E->control.data_file);
    E->filed[13]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.fas670",E->control.data_file);
    E->filed[14]=open(output_file,O_RDWR | O_CREAT, 0644);
    sprintf(output_file,"%s.stress",E->control.data_file);
    E->filed[9]=open(output_file,O_RDWR | O_CREAT, 0644);

    sprintf(in_file,"%s.temp",E->control.data_file1);
    fileid[5]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.velo",E->control.data_file1);
    fileid[10]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.topo_t",E->control.data_file1);
    fileid[11]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.topo_b",E->control.data_file1);
    fileid[12]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.visc",E->control.data_file1);
    fileid[13]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.fas670",E->control.data_file1);
    fileid[14]=open(in_file,O_RDONLY,0);
    sprintf(in_file,"%s.stress",E->control.data_file1);
    fileid[9]=open(in_file,O_RDONLY,0);


  for (i=1;i<=E->control.restart_frame;i++)  {
    ii = E->mesh.nsf;
    size2 = (2*ii+2)*sizeof(float);
    read(fileid[10],E->T,size2);
    write(E->filed[10],E->T,size2);

    size2 = (E->mesh.nsf+1)*sizeof(float);
    read(fileid[11],E->T,size2);
    write(E->filed[11],E->T,size2);
    read(fileid[12],E->T,size2);
    write(E->filed[12],E->T,size2);

    size2 = (E->mesh.nel+1)*sizeof(float);
    read(fileid[13],E->T,size2);
    write(E->filed[13],E->T,size2);

    size2 = (E->mesh.nsf+1)*sizeof(float);
    read(fileid[14],E->T,size2);
    write(E->filed[14],E->T,size2);

    size2 = (12*E->mesh.nsf+1)*sizeof(float);
    read(fileid[9],E->T,size2);
    write(E->filed[9],E->T,size2);

    size2= (E->mesh.nno+1)*sizeof(float);
    read(fileid[5],E->T,size2);
    write(E->filed[5],E->T,size2);

    }

    close(fileid[5]);
    close(fileid[10]);
    close(fileid[11]);
    close(fileid[12]);
    close(fileid[13]);
    close(fileid[14]);
    close(fileid[9]);
*/
  return;
  }
