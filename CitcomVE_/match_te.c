/* find the proper Te for a data file */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[])
{
  char *int2str(char *nchar, int jin);
  int nor,nox,i,j,n,p,Te,ti,visci,E,t1;
  double *topo_s=malloc(1001*sizeof(double)),*r=malloc(1001*sizeof(double));
  double *topo_s1=malloc(1001*sizeof(double));
  double blah,r1,dr,top,bottom,error,topo_ste,topo_ste1,boxl;
  char *first,*second=0,*age=0,*time=0,*visc=0;
  char filename[20]={'o','u','t','_','h','a','n','k','e','l','.','\0'};
  char **filename1 = argv+1;
  FILE *infile;
  FILE *telog = fopen("Te_log.dat", "a");
  FILE *tefile;

  first = strchr(filename,'\0');
  nor = 501;
  boxl=1500;
  nox=192;
  error = 100000;
  Te = 0;

  if(topo_s == NULL || r == NULL)
    {
      printf("Couldn't allocate memory at %lf",topo_s);
      return 0;
    }

  infile = fopen(*filename1, "r");
  printf("Filename is %s\n",*filename1);
  /* read data from files */
  for(i=1;i<=nor;i++)
    {
      fscanf(infile,"%i%lE%lE%lE",&n,&r[i],&blah,&topo_s1[i]);
    }

  /* Need to set up FE data in coords of SA data */
  topo_s[1]=topo_s1[1];
  j=2;
  for(i=2;i<=nor;i++)
    {
      if(((i-1)/(nor-1)*boxl)<=r[i])
	{
	  topo_s[i] = (topo_s1[j]-topo_s1[j-1])/(r[j]-r[j-1])
	    *(((i-1)/(nor-1)*boxl)-r[j-1])+topo_s1[j-1];
	}
      else
	{
	  j++;
	  topo_s[i] = (topo_s1[j]-topo_s1[j-1])/(r[j]-r[j-1])
	    *(((i-1)/(nor-1)*boxl)-r[j-1])+topo_s1[j-1];
	}
    }
  

  printf("Read %s\n",*filename1);
  dr = r[2] - r[1];
  
  for(j=1;j<=110;j=j+1)
    {
      first = int2str(first,j);
      tefile = fopen(filename, "r");  
      top = 0;
      bottom = 0;
      for(i=1;i<=nor;i++)
	{
	  fscanf(tefile,"%i%lE%lE%lE",&n,&r1,&blah,&topo_ste);
	  if(r1 != r[i])
	    {
	      printf("r's do not match for i = %i in %s\nquitting\n"
		     ,i,filename);
	      return 0;
	    }
	  if(i > 1)
	    /* integrate non-weighted error by trap method */
	    /*
	    {
	      top = top + 
		0.5*dr*(fabs(topo_s[i-1] - topo_ste1) 
			+ fabs(topo_s[i] - topo_ste));
	      bottom = bottom + 
		0.5*dr*(fabs(topo_s[i-1]) + fabs(topo_s[i]));
	    }
	    */
	    /* integrate area weighted error by trap method */
	    {
	      top = top + 
		0.5*dr*(r[i]-dr/2)*(fabs(topo_s[i-1] - topo_ste1) 
			+ fabs(topo_s[i] - topo_ste));
	      bottom = bottom + 
		0.5*dr*(r[i]-dr/2)*(fabs(topo_s[i-1]) + fabs(topo_s[i]));
	    }
	  topo_ste1 = topo_ste;
	}
      fclose(tefile);
      if((top/bottom) < error)
	{
	  error = top/bottom;
	  Te = j;
	}
    }
  printf("%s Te is %i with error %lf\n",*filename1,Te,error);
  fprintf(telog,"%s Te = %i error = %lf\n",*filename1,Te,error);
fclose(infile);
fclose(telog);
return 0;
}

/* int2str converts a 3 digit int to a string and returns pointer to string 
 * this was a real pain in the ass, but it seems to work 
 */
char *int2str(char *nchar, int jin)
{
  if(jin>=100)
    {
      *nchar = 48+floor(jin/100);
      *(nchar+1) = 48+floor((jin-floor(jin/100)*100)/10);
      *(nchar+2) = 48+fmod(jin,10);
      *(nchar+3) = '\0';
    }
  else if(jin>=10)
    {
      *nchar = 48+floor(jin/10);
      *(nchar+1) = 48+fmod(jin,10);
      *(nchar+2) = '\0';
    }
  else
    {
      *nchar = 48+jin;
      *(nchar+1) = '\0';
    }
  return nchar;
}
