 	/* This file contains the definitions of variables which are passed as arguments */
	/* to functions across the whole filespace of CITCOM. #include this file everywhere !*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(__osf__) 
void *Malloc1();
#endif

#define Malloc0(a) Malloc1((a),__FILE__,__LINE__)


/* #define Malloc0 malloc */

#define LIDN 0x1
#define VBX 0x2
#define VBZ 0x4
#define VBY 0x8
#define TBX 0x10
#define TBZ 0x20
#define TBY 0x40
#define TZEDGE 0x80
#define TXEDGE 0x100
#define TYEDGE 0x200
#define VXEDGE 0x400
#define VZEDGE 0x800
#define VYEDGE 0x1000
#define INTX 0x2000
#define INTZ 0x4000
#define INTY 0x8000
#define SBX 0x10000
#define SBZ 0x20000
#define SBY 0x40000
#define FBX 0x80000
#define FBZ 0x100000
#define FBY 0x200000

#define CBX 0x400000
#define CBZ 0x800000
#define CBY 0x1000000
#define HBX 0x2000000
#define HBZ 0x4000000
#define HBY 0x8000000

#define OFFSIDE 0x10000000

#define SKIP 0x1000000
#define SKIPID 0x1
#define ZEROID 0x2

#define RESTORE 0x1
#define REFINE1 0x1
#define REFINE2 0x2
#define GREFINE1 0x4
#define GREFINE2 0x8

#define LIDE 1

#ifndef COMPRESS_BINARY
#define COMPRESS_BINARY "/usr/bin/compress"
#endif

#define MAX_LEVELS 12
#define MAX_F    10
#define MAX_S    30

/* Macros */

#define max(A,B) (((A) > (B)) ? (A) : (B))
#define min(A,B) (((A) < (B)) ? (A) : (B))
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

typedef double higher_precision;  /* matrix coeffs etc */
typedef double higher_precision1; /* intermediate calculations for finding above coeffs */


/* Common structures */

struct Rect {
    int numb;
    char overlay[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float y1[40];
    float y2[40];
    float halfw[40];
    float mag[40];
} ;
 

struct Circ {
    int numb;
    char overlay[40];
    float x[40];
    float z[40];
    float y[40];
    float rad[40];
    float mag[40];
    float halfw[40];
};


struct Harm {
    int numb;
    int harms;
    char overlay[40];
    float off[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float y1[40];
    float y2[40];
    float kx[20][40];
    float kz[20][40];
    float ky[20][40];
    float ka[20][40];
    float phx[20][40];
    float phz[20][40];
    float phy[20][40];
};

struct Erfc {
 int numb;


};

struct RectBc {
    int numb;
    char norm[40];
    float intercept[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float halfw[40];
    float mag[40];
} ;
 

struct CircBc {
    int numb;
    char norm[40];
    float intercept[40];
    float x[40];
    float z[40];
    float rad[40];
    float mag[40];
    float halfw[40];
};


struct PolyBc {
    int numb;
    int order;
    char norm[40];
    float intercept[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float ax[20][40];
    float az[20][40];
};


struct HarmBc {
    int numb;
    int harms;
    char norm[40];
    float off[40];
    float intercept[40];
    float x1[40];
    float x2[40];
    float z1[40];
    float z2[40];
    float kx[20][40];
    float kz[20][40];
    float ka[20][40];
    float phx[20][40];
    float phz[20][40];
 };


struct Shape_function_dA  {
  double vpt[8];
  double spt[4];
  double ppt[1]; };

struct Shape_function1_dA  {
  double vpt[6*4];
  double ppt[6*1]; };

struct Shape_function1 	{ 
    double vpt[4*4];  /* node & gauss pt */
    double ppt[4*1];  };

struct Shape_function 	{ 
    double vpt[8*8];  /* node & gauss pt */
    double spt[8*4];  /* node & gauss pt */
    double ppt[8*1];  };

struct Shape_function_dx 	{ 
    double vpt[3*8*8]; /* dirn & node & gauss pt */
    double spt[3*8*4]; /* dirn & node & gauss pt */
    double ppt[3*8*1];  };

struct Shape_function1_dx 	{ 
    double vpt[2*4*4]; /* dirn & node & gauss pt */
    double ppt[2*4*1];  };

struct EG { 
    higher_precision g[24][1]; };

struct EK2 { 
    double k[8*8]; };

struct EK { 
    double k[24*24]; };

struct MEK { 
    double nint[9]; };
 
struct NK {
    higher_precision *k;
    int *map;
};

struct COORD {
    double centre[4];
    double size[4];
    double recip_size[4];
    double area;   } ;

struct SUBEL { 
    int sub[9];   };
			
struct ID  { 
    int doff[6];	}; /* can  be 1 or 2 or 3 */
struct IEN {
    int node[9];	};
struct FNODE {
    double node[9];	};
struct SIEN {
    int node[5];	};
struct LM  { 
    struct { int doff[4]; } node[9]; } ;

struct NEI {
    int *nels;
    int *lnode;
    int *element; };

struct BOUND  { 
    int bound[8];	}; 

struct IBM_DX   {
    float *x1;
    float *x2;
    int nox;
    int noz;
    };

struct MESH_DATA {/* general information concerning the fe mesh */ 
    int nsd;        /* Spatial extent 1,2,3d*/
    int dof;        /* degrees of freedom per node */
    int levmax;     
    int levmin;
    int levels;
    int mgunitx;
    int mgunitz;
    int mgunity;
    int NEQ[MAX_LEVELS];	/* All other values refer to the biggest mesh (& lid)  */
    int NNO[MAX_LEVELS];
    int NNOV[MAX_LEVELS];
    int NLNO[MAX_LEVELS];
    int NPNO[MAX_LEVELS];
    int NEL[MAX_LEVELS];
    int NOX[MAX_LEVELS];
    int NOZ[MAX_LEVELS];
    int NOY[MAX_LEVELS];
    int NNX[MAX_LEVELS][4];
    int ELX[MAX_LEVELS];
    int ELZ[MAX_LEVELS];
    int ELY[MAX_LEVELS];
    int LNDS[MAX_LEVELS];
    int LELS[MAX_LEVELS];
    int neqd;
    int neq;
    int nno;
    int nnov;
    int nlno;
    int npno;
    int nel;
    int snel;
    int elz;
    int ely;
    int nnx[4]; /* general form of ... */
    int nox;
    int elx;
    int noz;
    int noy;
    int *exs;
    int ezs;
    int eys;
    int nxs;
    int nzs;
    int nys;
    int nmx;
    int nsf; /* nodes for surface observables */
    int toptbc,topcbc;
    int bottbc,botcbc;
    int topvbc;
    int botvbc;
    int sidevbc;

    char topvbc_file[100];
    char botvbc_file[100];
    char sidevbc_file[100];
    char gridfile[4][100];


    int periodic_x;
    int periodic_y;
    double layer[4];			/* dimensionless dimensions */
    float lidz;
    float bl1width[4],bl2width[4],bl1mag[4],bl2mag[4];
    float hwidth[4],magnitude[4],offset[4],width[4]; /* grid compression information */ 
    int fnodal_malloc_size;
    int dnodal_malloc_size;
    int feqn_malloc_size;
    int deqn_malloc_size;
    int bandwidth;
    int null_source;
    int null_sink;
    int matrix_size[MAX_LEVELS];

} ;

struct HAVE {    /* horizontal averages */
    float *T;
    float *Vi;
    float *Rho;
    float *f;
    float *F;
    float *vrms;
    float *V[4];
};

struct SLICE {    /* horizontally sliced data, including topography */
    double *load[10];
    double *surf[10];
    float *tpg;
    float *tpgb;
    float *grv;
    float *geo;
    float *geok;
    float *grvk;
    float *grvb;
    float *geob;
    float *geobk;
    float *grvbk;
    float *tpgk;
    float *tpgbk;
    float *shflux;
    float *bhflux;
    float *cen_mflux;
    float *vxsurf[3];    /* surface velocity vectors */
    float *vline;        /* for kernels, velocity at force term */
    float *vlinek;
    float *tpglong;
    float *tpgblong;
    float *grvlong;
    float *geolong;
    float *grvblong;
    float *geoblong;

    int minlong;
    int maxlong;
  };


struct SPHERE  {
    double radius;
    double ro;
    double ri;
    double rcomp;

  int nox;
  int nsf;
  int elx;
  int snel;
  int llmax;
  int output_llmax;
  int hindex[1000];
  int hindice;
  int spectl_field;
  int axis_symmetry;

  double *sphc;
  double *sphs;
  double *field;

  double *sx[4];
  double *con;
  double *tableplm[964];
  double tablesint[964];

  double *tableplm_n[482];
  double tablesint_n[482];
  struct SIEN *sien;

    };

struct BAVE {
    float T;
    float Vi;
    double V[4]; };


struct TOTAL {
    float melt_prod;  };

struct MONITOR {
    char node_output[100][6];  /* recording the format of the output data */
    char sobs_output[100][6];  /* recording the format of the output data */
    int node_output_cols;
    int sobs_output_cols;

    int solution_cycles;

    float  time_scale;
    float  length_scale;
    float  viscosity_scale;
    float  geoscale;
    float  tpgscale;
    float  grvscale;
  
    float  delta_v_last_soln;
    float  elapsed_time;
    float  elapsed_time_vsoln;
    float  elapsed_time_vsoln1;
    float  reference_stress;
    float  incompressibility;
    float  vdotv;
    float  nond_av_heat_fl;  
    float  nond_av_adv_hfl;  
    float  cpu_time_elapsed;
    float  cpu_time_on_vp_it;
    float  cpu_time_on_forces;
    float  cpu_time_on_mg_maps;
    float  tpgkmag;
    float  grvkmag;
   
    float  Nusselt;
    float  Vmax;
    float  Vsrms;
    float  Vrms;
    float  Vrms_surface;
    float  Vrms_base;
    float  F_surface;
    float  F_base;
    float  Frat_surface;
    float  Frat_base;
    float  T_interior;
    float  T_maxvaried;
    float  Sigma_max;
    float  Sigma_interior;
    float  Vi_average;
   
};

struct CONTROL {
    int PID;

    char output_written_external_command[500];   /* a unix command to run when output files have been created */

    int ORTHO,ORTHOZ;   /* indicates levels of mesh symmetry */
    char B_is_good[MAX_LEVELS];  /* general information controlling program flow */
    char Ahat_is_good[MAX_LEVELS];  /* general information controlling program flow */
    char old_P_file[100];
    char data_file[100];
    char data_file1[100];

    char which_data_files[1000];
    char which_horiz_averages[1000];
    char which_running_data[1000];
    char which_observable_data[1000];
  
    char PROBLEM_TYPE[20]; /* one of ... */
    int KERNEL;
    int stokes;
    int pure_visc;
    int CONVECTION;
    int SLAB;	
    char GEOMETRY[20]; /* one of ... */
    int CART2D;
    int CART2pt5D;
    int CART3D;
    int AXI;	 
    char SOLVER_TYPE[20]; /* one of ... */
    int DIRECT;
    int CONJ_GRAD;
    int NMULTIGRID;
    int EMULTIGRID;
    int DIRECTII;
    char NODE_SPACING[20]; /* turns into ... */
    int GRID_TYPE;
    int COMPRESS;
    int DX;
    int CONMAN;
 
    int composition;
    float comp_diff;
 
  

    int dfact;
    double penalty;
    int augmented_Lagr;
    double augmented;
    int macroele;
    int faults;
    int NASSEMBLE;
    int comparison;
    int crust;
    int restart;
    int restart_frame;
    float plate_vel;
    float plate_age;

    float tole_comp;
  
    float sob_tolerance;
 
    float Ra_temp,Ra_comp; 
    float Ra_670,clapeyron670,transT670,width670; 
    float Ra_410,clapeyron410,transT410,width410; 
    float VBXtopval;
    float VBXbotval;
    float VBYtopval;
    float VBYbotval;
    float VBZtopval;
    float VBZbotval;

    float TBCtopval,CBCtopval;
    float TBCbotval,CBCbotval;

    float Q0;
    
    int precondition;
    int vprecondition;
    int keep_going;
    int v_steps_low;
    int v_steps_high;
    int v_steps_upper;
    int max_vel_iterations;
    int p_iterations;
    int max_same_visc;
    float max_res_red_each_p_mg;
    float sub_stepping_factor;
    int mg_cycle;
    int true_vcycle;
    int down_heavy;
    int up_heavy;
    int depth_dominated;
    int eqn_viscosity;
    int eqn_zigzag;
    int verbose;
    double accuracy;
    double vaccuracy; 
   
    int total_iteration_cycles;
    int total_v_solver_calls;
    
    int record_every;
    int record_all_until;

    int print_convergence;
    int sdepv_print_convergence;

     /* modules */
    int MELTING_MODULE;
    int CHEMISTRY_MODULE;
};

struct LOADA {
    double *mag;
    float *time;
    int stage;
    };

struct DATA {  
    float  rotation_rate;
    float  layer_km;
    float   grav_acc;
    float   therm_exp;
    float   Cp;
    float  therm_diff;
    float  therm_cond;
    float   density;
    float  res_density;
    float  res_density_X;
    float   melt_density;
    float   density_above;
    float   density_below;
    float   gas_const;
    float   surf_heat_flux;
    float  ref_viscosity;
    float   melt_viscosity;
    float   permeability;
    float   grav_const;
    float  surf_temp;
    float   youngs_mod; 
    float   shear_mod; 
    float   Te;
    float   ref_temperature;
    float   Tsurf;
    float   T_sol0;
    float   botm_scaling;
    float   surf_scaling;
    float   delta_S;
    float   dTsol_dz;
    float   dTsol_dF;
    float   dT_dz; };
	      
struct All_variables {     
#include "Convection_variables.h"
#include "viscosity_descriptions.h"
#include "temperature_descriptions.h"
#include "advection.h"

    FILE *fp;
    FILE *filed[20];
    struct HAVE Have;
    struct BAVE Bulkave;
    struct TOTAL Total;
    struct MESH_DATA mesh;
    struct CONTROL control;
    struct MONITOR monitor;
    struct DATA data;
    struct SLICE slice;
    struct COORD *eco;
    struct LOADA load;
    struct IBM_DX ibm_dx;
    struct IEN *ien;  /* global */
    struct SIEN *fault_ien;
    struct SIEN *sien;
    struct ID *id;
    struct COORD *ECO[MAX_LEVELS];
    struct IEN *IEN[MAX_LEVELS]; /* global at each level */
    struct FNODE *TWW[MAX_LEVELS]; /* */
    struct ID *ID[MAX_LEVELS];
    struct NEI NEI[MAX_LEVELS];
    struct SUBEL *EL[MAX_LEVELS];
    struct EG *elt_del[MAX_LEVELS];
    struct EK *elt_k[MAX_LEVELS];
    struct SPHERE sphere;

    higher_precision *Eqn_k[MAX_LEVELS];  
    int *Node_map[MAX_LEVELS];
    int *Node_eqn[MAX_LEVELS];
    int *Node_k_id[MAX_LEVELS];


    double *BI[MAX_LEVELS];      /* inv of  diagonal elements of K matrix */
    double *BPI[MAX_LEVELS];
    float *V[4];                 /* velocity X[dirn][node] can save memory */
    float *V1[4];                 /* velocity X[dirn][node] can save memory */
    float *Vest[4];                 /* velocity X[dirn][node] can save memory */

    double *P,*PO,*F,*H,*S,*U;
    float *Psi;
    float *NP;
    float *edot;               /* strain rate invariant */
    double *MASS[MAX_LEVELS],*Mass;               /* lumped mass matrix (diagonal elements) for p-g solver etc. */  
    float *tw;
    float *stress;
    double *X[4],*XX[MAX_LEVELS][4],*X0[4],*Interp[MAX_LEVELS][4];
    float *ZZ;
    float *T,*C;		
    double *buoyancy;		
    float *Tdot,*Cdot;		
    float *Fas670,*Fas410;		
    float *Fas670_b,*Fas410_b;		

    double *Maxwelltime,*MaxwelltimeO;

    float *Exx,*Ezz,*Exz,*Eyy;
    double *Sxx,*Szz,*Sxz,*Syy;
    double *Xsurf[3];
    double *ShearMO;

    float *EViO,*EVold,*EVolder,*S2xx,*S2zz,*S2xz,*S2yy;

    float *Vi,*EVi;
    float *VI[MAX_LEVELS];	/* viscosity has to soak down to all levels */
    float *EVI[MAX_LEVELS];	/* element viscosity has to soak down to all levels */
    float *VB[4],*TB[4],*CB[4];/* boundary conditions for V,T defined everywhere */
    float *TW[MAX_LEVELS];	/* nodal weightings */

    int num_zero_resid[MAX_LEVELS];	       
    int *zero_resid[MAX_LEVELS];	       
    int *surf_element;	       
    int *surf_node;	       
    int *mat;	        /* properties of mat */
    unsigned int *node;	        /* properties of node */
    unsigned int *NODE[MAX_LEVELS];
    unsigned int *ELEMENT[MAX_LEVELS];
    unsigned int *eqn;

    unsigned int *EQN[MAX_LEVELS];
    higher_precision *B_R[MAX_LEVELS];  

    double **global_K;  /* direct solver stuff */
    double **factor_K;
    double *global_F;
    struct LM *lmd;
    struct LM *lm;
    struct LM *LMD[MAX_LEVELS];
  
    struct Shape_function1 M; /* master-element shape funtions */
    struct Shape_function1_dx Mx; 
    struct Shape_function N;
    struct Shape_function_dx Nx;
    struct Shape_function1 L; /* master-element shape funtions */
    struct Shape_function1_dx Lx; 
 
    void (* build_forcing_term)();
    void (* iterative_solver)();
    void (* next_buoyancy_field)();
    void (* obtain_gravity)();
    void (* problem_settings)();
    void (* problem_derived_values)();
    void (* problem_allocate_vars)();
    void (* problem_boundary_conds)();
    void (* problem_node_positions)();
    void (* problem_update_node_positions)();
    void (* problem_initial_fields)();
    void (* problem_update_bcs)();
    void (* special_process_new_velocity)();
    void (* special_process_new_buoyancy)();
    void (* solve_stokes_problem)(); 
    void (* solver_allocate_vars)(); 
    void (* transform)();

    float (* node_space_function[3])();
 
};
