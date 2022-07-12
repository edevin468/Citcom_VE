 /*CITCOM: A finite element convection program written at Caltech 1992 */
 /*Aims to include an iterative matrix solver based on Multigrid techniques */
 /*To do this requires the use of a mixed method and a conjugate-gradient */
 /*approach to determining the */

#include <math.h>
#include <malloc.h>
#include <sys/types.h>

#include "element_definitions.h"
#include "global_defs.h"

extern int Emergency_stop;

main(argc,argv)
     int argc;
     char **argv;
     
{	/* Functions called by main*/
  void general_stokes_solver();
  void read_instructions();
  void solve_constrained_flow();
  void solve_derived_velocities();
  void process_temp_field(); 
  void initial_again();

  float dot();

  int i,k, *temp;
  double CPU_time0(),time,time_a,time_b,time_c,initial_time,start_time;
 
  struct All_variables E;

 E.monitor.solution_cycles=0;

    start_time = time = CPU_time0();
 
  E.advection.timestep = 0.0;

  read_instructions(&E,argc,argv);


  fprintf(stderr,"Input parameters taken from file '%s'\n",argv[1]);
  fprintf(stderr,"Initialization complete after %g seconds\n\n",CPU_time0()-time); fflush(E.fp);
  initial_time = CPU_time0()-time;
  E.monitor.cpu_time_on_vp_it = CPU_time0();

 i=0;
 time_a = 320;

 do {

  i ++;

  E.mesh.nmx = i;

  initial_again(&E);

  E.monitor.solution_cycles=0;
  E.monitor.elapsed_time = 0.0;
  E.advection.timestep = 0.0;
  E.control.keep_going=1;

  time_b=time_a/pow((double)2.0,(double)(i-1)); 
  
  E.data.T_sol0=time_b;

  fprintf(stderr,"ttttt %lf %lf %d\n",time_a,time_b,i);
  fprintf(E.fp,"ttttt %lf %lf %d\n",time_a,time_b,i);

  if (time_b<time_a*0.4 && time_b>time_a*0.1) {
    E.advection.max_timesteps=6001;
    E.control.record_every=60;  
    }
  else if (time_b<=time_a*0.1 && time_b>time_a*0.04) {
    E.advection.max_timesteps=3001;
    E.control.record_every=30; 
    }
  else if (time_b<=time_a*0.04 && time_b>time_a*0.004) {
    E.advection.max_timesteps=1501;
    E.control.record_every=15;  
    }
  else if (time_b<=time_a*0.004) {
    E.advection.max_timesteps=601;
    E.control.record_every=6;  
    }

  time_c = time_b/(E.advection.max_timesteps-1);
  E.data.delta_S = time_c;
  fprintf(stderr,"tttttaaa %g %d\n",time_c,E.advection.max_timesteps);
  fprintf(E.fp,"tttttaaa %g %d\n",time_c,E.advection.max_timesteps);


  general_stokes_solver(&E);

  if (E.control.stokes)  {
     process_new_velocity(&E,E.monitor.solution_cycles);
     E.control.keep_going=0;
     E.monitor.solution_cycles++; 
     }

  else    {

     E.advection.timestep = E.data.delta_S;
     process_new_velocity(&E,E.monitor.solution_cycles);

    while ( E.control.keep_going   &&  (Emergency_stop == 0) )   {

      E.monitor.solution_cycles++; 
      if(E.monitor.solution_cycles>E.control.print_convergence)
         E.control.print_convergence=1;

      if(E.monitor.solution_cycles>E.advection.max_timesteps)
         E.control.keep_going = 0;

/*      (E.next_buoyancy_field)(&E);    */

      E.monitor.elapsed_time += E.advection.timestep;

/*      process_temp_field(&E,E.monitor.solution_cycles); 
*/
           /* determine MatProp and solve incremental displacement 
                          for this time step */
      general_stokes_solver(&E);

           /* determine stress at t and update cumulated stress, and update grid */
      process_new_velocity(&E,E.monitor.solution_cycles);

        fprintf(E.fp,"CPU total = %g & CPU = %g for step %d time = %.4e dt = %.4e \n",CPU_time0()-start_time,CPU_time0()-time,E.monitor.solution_cycles,E.monitor.elapsed_time,E.advection.timestep);
        time = CPU_time0();

        }

     }
   }while (time_b>0.00008);
  // }while (i<1);  /* use this to run for only one period (specified by time_a). timesteps will be those specified in the input file.  */
  
     E.monitor.cpu_time_on_vp_it=CPU_time0()-E.monitor.cpu_time_on_vp_it;
     fprintf(E.fp,"Initialization overhead = %f\n",initial_time);
     fprintf(E.fp,"Average cpu time taken for velocity step = %f\n",
	 E.monitor.cpu_time_on_vp_it/((float)(E.monitor.solution_cycles)));
     fprintf(stderr,"Average cpu time taken for velocity step = %f\n",
	 E.monitor.cpu_time_on_vp_it/((float)(E.monitor.solution_cycles)));

  fclose(E.fp);

  return;  

  } 
