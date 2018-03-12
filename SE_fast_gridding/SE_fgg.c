#include "SE_fgg.h"
#include "x86intrin.h"
#include "string.h"
#include "malloc.h"

#include "fgg_thrd.h"
#include "fgg_thrd.c"

#include "SE_fkg.c"

// Dag Lindbo, dag@kth.se

// =============================================================================
// Internal routines ===========================================================

// -----------------------------------------------------------------------------
static double randnum(double min, double L)
{
    double q = ( (double) rand() )/RAND_MAX;
    return L*q+min;
}

// -----------------------------------------------------------------------------
int is_aligned(double* h, int alignment)
{
    if ((((unsigned long)h) & (alignment-1)) == 0)
        return 1;
    else
        return 0;
}

// =============================================================================
// SE FGG Utility routines =====================================================

// -----------------------------------------------------------------------------
// Set array elements to double-precision zero
void SE_fp_set_zero(double* restrict x, const int N)
{
    for(int i=0; i<N; i++)
	x[i] = 0.0;
}

// -----------------------------------------------------------------------------
int SE_prod3(const int v[3])
{
    return v[0]*v[1]*v[2];
}

// -----------------------------------------------------------------------------
double SE_gettime(void) 
{
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec + 1e-6*tv.tv_usec;
}

// -----------------------------------------------------------------------------
void SE_FGG_pack_params(SE_FGG_params* params, int N, int M0, int M1, int M2, 
			int P, double c, double h)
{
    params->N = N;
    params->P = P;
    params->P_half=half(P);
    params->c = c;
    params->d = pow(c/PI,1.5);
    params->h = h;
    params->a = -1;
    params->b = -1;

    params->dims[0] = M0;
    params->dims[1] = M1;
    params->dims[2] = M2;

    params->npdims[0] = M0+P;
    params->npdims[1] = M1+P;
    params->npdims[2] = M2+P;
}

// -----------------------------------------------------------------------------
void SE2P_FGG_pack_params(SE_FGG_params* params, int N, int M0, int M1, int M2, 
			  int P, double c, double h, double a)
{
    params->N = N;
    params->P = P;
    params->P_half=half(P);
    params->c = c;
    params->d = pow(c/PI,1.5);
    params->h = h;
    params->a = a;
    params->b = -1;

    params->dims[0] = M0;
    params->dims[1] = M1;
    params->dims[2] = M2;

    params->npdims[0] = M0+P;
    params->npdims[1] = M1+P;
    params->npdims[2] = M2+P;
}

// -----------------------------------------------------------------------------
void SE1P_FGG_pack_params(SE_FGG_params* params, int N, int M0, int M1, int M2, 
			  int P, double c, double h, double a, double b)
{
    params->N = N;
    params->P = P;
    params->P_half=half(P);
    params->c = c;
    params->d = pow(c/PI,1.5);
    params->h = h;
    params->a = a;
    params->b = b;

    params->dims[0] = M0;
    params->dims[1] = M1;
    params->dims[2] = M2;

    params->npdims[0] = M0+P;
    params->npdims[1] = M1+P;
    params->npdims[2] = M2+P;
}

// -----------------------------------------------------------------------------
void SE_FGG_allocate_workspace(SE_FGG_work* work, const SE_FGG_params* params, 
			       int allocate_zs, int allocate_fgg_expa)
{
    const int P=params->P;
    int numel = SE_prod3(params->npdims);
    work->H = SE_FGG_MALLOC(numel*sizeof(double));

    SE_fp_set_zero(work->H, numel);

    if(allocate_zs)
	work->zs = SE_FGG_MALLOC(P*P*P*sizeof(double));
    else
	work->zs = NULL;

    work->free_zs=allocate_zs;
    
    if(allocate_fgg_expa)
    {
        numel = (params->N)*(params->P);
	work->zx = (double*) SE_FGG_MALLOC(numel*sizeof(double));
	work->zy = (double*) SE_FGG_MALLOC(numel*sizeof(double));
	work->zz = (double*) SE_FGG_MALLOC(numel*sizeof(double));
	work->idx= (int*) SE_FGG_MALLOC(params->N*sizeof(int));
#ifdef FORCE
        work->zfx = (double*) SE_FGG_MALLOC(numel*sizeof(double));
        work->zfy = (double*) SE_FGG_MALLOC(numel*sizeof(double));
        work->zfz = (double*) SE_FGG_MALLOC(numel*sizeof(double));
#endif
    }
    else
    {
	work->zx=NULL;
	work->zy=NULL;
	work->zz=NULL;
#ifdef FORCE
        work->zfx=NULL;
        work->zfy=NULL;
        work->zfz=NULL;
#endif
	work->idx=NULL;
    }
    work->free_fgg_expa=allocate_fgg_expa;
}


// -----------------------------------------------------------------------------
double* SE_FGG_allocate_grid(const SE_FGG_params* params)
{
    int numel = SE_prod3(params->dims);
    double* H_per = SE_FGG_MALLOC(numel*sizeof(double));
    SE_fp_set_zero(H_per, numel);
    return H_per;
}

// -----------------------------------------------------------------------------
double* SE_FGG_allocate_vec(const int Nm)
{
    double* phi = SE_FGG_MALLOC(Nm*sizeof(double));
    SE_fp_set_zero(phi, Nm);
    return phi;
}


// -----------------------------------------------------------------------------
void SE_FGG_free_workspace(SE_FGG_work* work)
{
    SE_FGG_FREE(work->H);

    if(work->free_zs)
	SE_FGG_FREE(work->zs);

    if(work->free_fgg_expa)
    {
	SE_FGG_FREE(work->zx);
	SE_FGG_FREE(work->zy);
	SE_FGG_FREE(work->zz);
#ifdef FORCE
        SE_FGG_FREE(work->zfx);
        SE_FGG_FREE(work->zfy);
        SE_FGG_FREE(work->zfz);
#endif
	SE_FGG_FREE(work->idx);
    }
}

// -----------------------------------------------------------------------------
void SE_init_unit_system(SE_state* s, const SE_FGG_params* params)
{
    const int N=params->N;

    if(N%2!=0)
	return;

    s->x = SE_FGG_MALLOC(3*N*sizeof(double));
    s->q = SE_FGG_MALLOC(  N*sizeof(double));
    s->phi = SE_FGG_MALLOC(N*sizeof(double));
   
    FILE *fp;
    fp = fopen("atoms.txt","r");
    int err=0,i;

    for(i=0; i<N; i++)
    {
        err = fscanf(fp,"%lf",&(s->x[i    ]));
        err|= fscanf(fp,"%lf",&(s->x[i+  N]));
        err|= fscanf(fp,"%lf",&(s->x[i+2*N]));
        err|= fscanf(fp,"%lf",&(s->q[i]));
    }

    fclose(fp);
    if(err==EOF)
	printf("Not successful!!!\n");
}

// -----------------------------------------------------------------------------
void SE_init_system(SE_state* s, const SE_FGG_params* params)
{
    const int N=params->N;

    s->x = SE_FGG_MALLOC(3*N*sizeof(double));
    s->q = SE_FGG_MALLOC(  N*sizeof(double));
    s->phi = SE_FGG_MALLOC(N*sizeof(double));


    for(int i=0; i<N; i++)
    {
	s->q[i] = randnum(-1,2);

	s->x[i    ] = randnum(0,1);
	s->x[i+N  ] = randnum(0,1);
	s->x[i+2*N] = randnum(0,1);
    }
}

// -----------------------------------------------------------------------------
void SE_free_system(SE_state* s)
{
    SE_FGG_FREE(s->x);
    SE_FGG_FREE(s->q);
#ifdef CALC_ENERGY
    SE_FGG_FREE(s->phi);
#endif
}

// -----------------------------------------------------------------------------
// Wrap H to produce periodicity
// OUTPUT IN FORTRAN/MATLAB-STYLE COLUMN MAJOR LAYOUT!
void SE_FGG_wrap_fcn(double* H_per, 
		     const SE_FGG_work* work, const SE_FGG_params* params)
{
    int idx;
    int widx[3];
    const int p_half = half(params->P);

    // can not openMP here, race to += on H_per beacuse indices wrap around
    for(int i=0; i<params->npdims[0]; i++)
    {
	for(int j=0; j<params->npdims[1]; j++)
	{
	    for(int k=0; k<params->npdims[2]; k++)
	    {
		widx[0] = vmod(i-p_half,params->dims[0]);
		widx[1] = vmod(j-p_half,params->dims[1]);
		widx[2] = vmod(k-p_half,params->dims[2]);
		idx = __IDX3_CMAJ(widx[0], widx[1], widx[2], 
				  params->dims[0], params->dims[1]);
		H_per[idx] += work->H[ __IDX3_RMAJ(i,j,k,
						   params->npdims[1],
						   params->npdims[2])];
	    }
	}
    }
}

// -----------------------------------------------------------------------------
// Wrap H to produce 2-periodicity
// OUTPUT IN FORTRAN/MATLAB-STYLE COLUMN MAJOR LAYOUT!
void SE2P_FGG_wrap_fcn(double* restrict H_per, 
		       const SE_FGG_work* work, 
		       const SE_FGG_params* params)
{
    int idx;
    int widx[2];
    const int p_half = half(params->P);

    // can not openMP here, race to += on H_per beacuse indices wrap around
    for(int i=0; i<params->npdims[0]; i++)
    {
	for(int j=0; j<params->npdims[1]; j++)
	{
	    for(int k=0; k<params->npdims[2]; k++)
	    {
		widx[0] = vmod(i-p_half,params->dims[0]);
		widx[1] = vmod(j-p_half,params->dims[1]);
		idx = __IDX3_CMAJ(widx[0], widx[1], k, 
				  params->dims[0], params->dims[1]);
		H_per[idx] += work->H[ __IDX3_RMAJ(i,j,k,
						   params->npdims[1],
						   params->npdims[2]) ];
	    }
	}
    }
}

// -----------------------------------------------------------------------------
// Wrap H to produce 1-periodicity
// OUTPUT IN FORTRAN/MATLAB-STYLE COLUMN MAJOR LAYOUT!
void SE1P_FGG_wrap_fcn(double* restrict H_per, 
		       const SE_FGG_work* work, 
		       const SE_FGG_params* params)
{
    int idx;
    int widx;
    const int p_half = half(params->P);

    // can not openMP here, race to += on H_per beacuse indices wrap around
    for(int i=0; i<params->npdims[0]; i++)
    {
      widx= vmod(i-p_half,params->dims[0]);
	for(int j=0; j<params->npdims[1]; j++)
	{
	    for(int k=0; k<params->npdims[2]; k++)
	    {
		idx = __IDX3_CMAJ(widx, j, k, 
				  params->dims[0], params->dims[1]);
		H_per[idx] += work->H[ __IDX3_RMAJ(i,j,k,
						   params->npdims[1],
						   params->npdims[2]) ];
	    }
	}
    }
}

// -----------------------------------------------------------------------------
// Extend periodic function larger box
// INPUT IN FORTRAN/MATLAB-STYLE COLUMN MAJOR LAYOUT!
void SE_FGG_extend_fcn(SE_FGG_work* work, const double* H_per, 
		       const SE_FGG_params* params)
{
    int idx;
    int widx[3];
    const int p_half = params->P_half;

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int i=0; i<params->npdims[0]; i++)
    {
	for(int j=0; j<params->npdims[1]; j++)
	{
	    for(int k=0; k<params->npdims[2]; k++)
	    {
		widx[0] = vmod(i-p_half,params->dims[0]);
		widx[1] = vmod(j-p_half,params->dims[1]);
		widx[2] = vmod(k-p_half,params->dims[2]);
		idx = __IDX3_CMAJ(widx[0], widx[1], widx[2], 
				  params->dims[0], params->dims[1]);
		work->H[__IDX3_RMAJ(i,j,k,params->npdims[1],params->npdims[2])]
		    = H_per[idx];
	    }
	}
    }
}

// -----------------------------------------------------------------------------
// Extend 2-periodic function larger box
// INPUT IN FORTRAN/MATLAB-STYLE COLUMN MAJOR LAYOUT!
void SE2P_FGG_extend_fcn(SE_FGG_work* work, const double* H_per, 
			 const SE_FGG_params* params)
{
    int idx;
    int widx[2];
    const int p_half = half(params->P);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int i=0; i<params->npdims[0]; i++)
    {
	for(int j=0; j<params->npdims[1]; j++)
	{
	    for(int k=0; k<params->npdims[2]; k++)
	    {
		widx[0] = vmod(i-p_half,params->dims[0]);
		widx[1] = vmod(j-p_half,params->dims[1]);
		idx = __IDX3_CMAJ(widx[0], widx[1], k, 
				  params->dims[0], params->dims[1]);
		work->H[__IDX3_RMAJ(i,j,k,params->npdims[1],params->npdims[2])]
		    = H_per[idx];
	    }
	}
    }
}

// -----------------------------------------------------------------------------
// Extend 1-periodic function larger box
// INPUT IN FORTRAN/MATLAB-STYLE COLUMN MAJOR LAYOUT!
void SE1P_FGG_extend_fcn(SE_FGG_work* work, const double* H_per, 
			 const SE_FGG_params* params)
{
    int idx;
    int widx;
    const int p_half = half(params->P);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int i=0; i<params->npdims[0]; i++)
    {
      widx = vmod(i-p_half,params->dims[0]);
	for(int j=0; j<params->npdims[1]; j++)
	{
	    for(int k=0; k<params->npdims[2]; k++)
	    {
		idx = __IDX3_CMAJ(widx, j, k, 
				  params->dims[0], params->dims[1]);
		work->H[__IDX3_RMAJ(i,j,k,params->npdims[1],params->npdims[2])]
		    = H_per[idx];
	    }
	}
    }
}

// =============================================================================
// Core SE FGG routines ========================================================

// -----------------------------------------------------------------------------
void SE_FGG_base_gaussian(SE_FGG_work* work, const SE_FGG_params* params)
{
    int idx;
    
    // unpack parameters
    const int p=params->P;
    const int p_half = half(p);
    const int p_from = (is_odd(p) ? p_half:p_half-1);
    const double h=params->h;
    const double c=params->c;
    const double d=params->d;

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int i = -p_from; i<=p_half; i++)
    {
	// hoisting this index calculation (more) breaks omp-parallel code
	idx = __IDX3_RMAJ(i+p_from, 0, 0, p, p);
	for(int j = -p_from; j<=p_half; j++)
	    for(int k = -p_from; k<=p_half; k++)
	    {
		work->zs[idx++] = d*exp(-c*((i*h)*(i*h) + 
					    (j*h)*(j*h) + 
					    (k*h)*(k*h)));
	    }
    }
}

// -----------------------------------------------------------------------------
#ifdef THREE_PERIODIC
void fgg_offset_3p(const double x[3],
		   const SE_FGG_params* params,
		   double t0[3],
		   int idx_from[3])
{
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    int idx;
    if(is_odd(p))
    {
	for(int j=0; j<3; j++)
	{
	    idx = (int) round(x[j]/h);
	    idx_from[j] = idx - p_half;
	    t0[j] = x[j]-h*idx;
	}
    }
    else
    {
	for(int j=0; j<3; j++)
	{
	    idx = (int) floor(x[j]/h);
	    idx_from[j] = idx - (p_half-1);
	    t0[j] = x[j]-h*idx;
	}
    }
}

static inline
int fgg_expansion_3p(const double x[3], const double q,
		     const SE_FGG_params* params,
		     double z2_0[P_MAX], 
		     double z2_1[P_MAX], 
		     double z2_2[P_MAX])
{
    // unpack params
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double c=params->c;
    
    double t0[3];
    int idx_from[3];

    // compute index range and centering
    fgg_offset_3p(x, params, t0, idx_from);

    // compute third factor 
    double z3 = exp(-c*(t0[0]*t0[0] + t0[1]*t0[1] + t0[2]*t0[2]) )*q;

    // compute second factor by induction
    double z_base0 = exp(2*c*h*t0[0]);
    double z_base1 = exp(2*c*h*t0[1]);
    double z_base2 = exp(2*c*h*t0[2]);

    double z0, z1, z2;
    if(is_odd(p))
    {
	z0 = pow(z_base0,-p_half);
	z1 = pow(z_base1,-p_half);
	z2 = pow(z_base2,-p_half);
    }	
    else
    {
    	z0 = pow(z_base0,-p_half+1);
    	z1 = pow(z_base1,-p_half+1);
    	z2 = pow(z_base2,-p_half+1);
    }

    z2_0[0] = z0;
    z2_1[0] = z1;
    z2_2[0] = z2;
    for(int i=1; i<p; i++)
    {
	z0 *=z_base0;
	z1 *=z_base1;
	z2 *=z_base2;

	z2_0[i] = z0;
	z2_1[i] = z1;
	z2_2[i] = z2;
    }

    // save some flops by multiplying one vector with z3 factor
    for(int i=0; i<p; i++)
    {
	z2_0[i] *= z3;
    }

    return __IDX3_RMAJ(idx_from[0]+p_half, 
		       idx_from[1]+p_half, 
		       idx_from[2]+p_half, 
		       params->npdims[1], params->npdims[2]);
}

static inline
int fgg_index_3p(const double x[3],
		 const SE_FGG_params* params)
{
    const int p_half = params->P_half;
    double t0[3];
    int idx_from[3];
    fgg_offset_3p(x, params, t0, idx_from);
    return __IDX3_RMAJ(idx_from[0]+p_half, 
		       idx_from[1]+p_half, 
		       idx_from[2]+p_half, 
		       params->npdims[1], params->npdims[2]);
}

#endif
// -----------------------------------------------------------------------------
#ifdef TWO_PERIODIC
// -----------------------------------------------------------------------------
static inline
int fgg_expansion_2p(const double x[3], const double q,
		     const SE_FGG_params* params,
		     double z2_0[P_MAX], 
		     double z2_1[P_MAX], 
		     double z2_2[P_MAX])
{
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double c=params->c;
    const double a=params->a;

    double t0[3];
    int idx;
    int idx_from[3];

    // compute index range and centering
    if(is_odd(p))
    {
	idx = (int) round(x[0]/h);
	idx_from[0] = idx - p_half;
	t0[0] = x[0]-h*idx;
	
	idx = (int) round(x[1]/h);
	idx_from[1] = idx - p_half;
	t0[1] = x[1]-h*idx;

	idx = (int) round((x[2]-(a+h/2))/h);
	idx_from[2] = idx - p_half;
	t0[2] = x[2] - (idx*h + (a+h/2));
    }
    else
    {
	idx = (int) floor(x[0]/h);
	idx_from[0] = idx - (p_half-1);
	t0[0] = x[0]-h*idx;

	idx = (int) floor(x[1]/h);
	idx_from[1] = idx - (p_half-1);
	t0[1] = x[1]-h*idx;

	idx = (int) floor((x[2]-(a+h/2))/h);
	idx_from[2] = idx - (p_half-1);
	t0[2] = x[2] - (idx*h + (a+h/2));
    }

    // compute third factor 
    double z3 = exp(-c*(t0[0]*t0[0] + t0[1]*t0[1] + t0[2]*t0[2]) )*q;

    // compute second factor by induction
    double z_base0 = exp(2*c*h*t0[0]);
    double z_base1 = exp(2*c*h*t0[1]);
    double z_base2 = exp(2*c*h*t0[2]);

    double z0, z1, z2;
    if(is_odd(p))
    {
	z0 = pow(z_base0,-p_half);
	z1 = pow(z_base1,-p_half);
	z2 = pow(z_base2,-p_half);
    }	
    else
    {
    	z0 = pow(z_base0,-p_half+1);
    	z1 = pow(z_base1,-p_half+1);
    	z2 = pow(z_base2,-p_half+1);
    }

    z2_0[0] = z0;
    z2_1[0] = z1;
    z2_2[0] = z2;
    for(int i=1; i<p; i++)
    {
	z0 *=z_base0;
	z1 *=z_base1;
	z2 *=z_base2;

	z2_0[i] = z0;
	z2_1[i] = z1;
	z2_2[i] = z2;
    }

    // save some flops by multiplying one vector with z3 factor
    for(int i=0; i<p; i++)
    {
	z2_0[i] *= z3;
    }

    return __IDX3_RMAJ(idx_from[0]+p_half, 
		       idx_from[1]+p_half, 
		       idx_from[2], 
		       params->npdims[1], params->npdims[2]);
}

static inline
int fgg_index_2p(const double x[3],
		 const SE_FGG_params* params)
{
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double a=params->a;

    int idx;
    int idx_from[3];

    // compute index range and centering
    if(is_odd(p))
    {
	idx = (int) round(x[0]/h);
	idx_from[0] = idx - p_half;
	
	idx = (int) round(x[1]/h);
	idx_from[1] = idx - p_half;

	idx = (int) round((x[2]-(a+h/2))/h);
	idx_from[2] = idx - p_half;
    }
    else
    {
	idx = (int) floor(x[0]/h);
	idx_from[0] = idx - (p_half-1);

	idx = (int) floor(x[1]/h);
	idx_from[1] = idx - (p_half-1);

	idx = (int) floor((x[2]-(a+h/2))/h);
	idx_from[2] = idx - (p_half-1);
    }

    return __IDX3_RMAJ(idx_from[0]+p_half, 
		       idx_from[1]+p_half, 
		       idx_from[2], 
		       params->npdims[1], params->npdims[2]);
}
#endif


// -----------------------------------------------------------------------------
#ifdef ONE_PERIODIC
static 
int fgg_expansion_1p(const double x[3], const double q,
		     const SE_FGG_params* params,
		     double z2_0[P_MAX], 
		     double z2_1[P_MAX], 
		     double z2_2[P_MAX])
{
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double c=params->c;
    const double a=params->a;
    const double b=params->b;

    double t0[3];
    int idx;
    int idx_from[3];

    // compute index range and centering
    if(is_odd(p))
    {
	idx = (int) round((x[0]-(a+h/2))/h);
	idx_from[0] = idx - p_half;
	t0[0] = x[0] - (idx*h + (a+h/2));

	idx = (int) round((x[1]-(b+h/2))/h);
	idx_from[1] = idx - p_half;
	t0[1] = x[1] - (idx*h + (b+h/2));

	idx = (int) round(x[2]/h);
	idx_from[2] = idx - p_half;
	t0[2] = x[2]-h*idx;
    }
    else
    {
	idx = (int) floor((x[0]-(a+h/2))/h);
	idx_from[0] = idx - (p_half-1);
	t0[0] = x[0] - (idx*h + (a+h/2));

	idx = (int) floor((x[1]-(b+h/2))/h);
	idx_from[1] = idx - (p_half-1);
	t0[1] = x[1] - (idx*h + (b+h/2));

	idx = (int) floor(x[2]/h);
	idx_from[2] = idx - (p_half-1);
	t0[2] = x[2]-h*idx;
    }

    // compute third factor 
    double z3 = exp(-c*(t0[0]*t0[0] + t0[1]*t0[1] + t0[2]*t0[2]) )*q;

    // compute second factor by induction
    double z_base0 = exp(2*c*h*t0[0]);
    double z_base1 = exp(2*c*h*t0[1]);
    double z_base2 = exp(2*c*h*t0[2]);

    double z0, z1, z2;
    if(is_odd(p))
    {
	z0 = pow(z_base0,-p_half);
	z1 = pow(z_base1,-p_half);
	z2 = pow(z_base2,-p_half);
    }	
    else
    {
    	z0 = pow(z_base0,-p_half+1);
    	z1 = pow(z_base1,-p_half+1);
    	z2 = pow(z_base2,-p_half+1);
    }

    z2_0[0] = z0;
    z2_1[0] = z1;
    z2_2[0] = z2;
    for(int i=1; i<p; i++)
    {
	z0 *=z_base0;
	z1 *=z_base1;
	z2 *=z_base2;

	z2_0[i] = z0;
	z2_1[i] = z1;
	z2_2[i] = z2;
    }

    // save some flops by multiplying one vector with z3 factor
    for(int i=0; i<p; i++)
    {
	z2_0[i] *= z3;
    }

    return __IDX3_RMAJ(idx_from[0], 
		       idx_from[1], 
		       idx_from[2]+p_half, 
		       params->npdims[1], params->npdims[2]);
}

// -----------------------------------------------------------------------------
static inline
int fgg_index_1p(const double x[3],
		 const SE_FGG_params* params)
{
    const int p = params->P;
    const int p_half = params->P_half;
    const double h = params->h;
    const double a=params->a;
    const double b=params->b;

    int idx;
    int idx_from[3];

    // compute index range and centering
    if(is_odd(p))
    {
	idx = (int) round((x[0]-(a+h/2))/h);
	idx_from[0] = idx - p_half;

	idx = (int) round((x[1]-(b+h/2))/h);
	idx_from[1] = idx - p_half;

	idx = (int) round(x[2]/h);
	idx_from[2] = idx - p_half;
    }
    else
    {
	idx = (int) floor((x[0]-(a+h/2))/h);
	idx_from[0] = idx - (p_half-1);

	idx = (int) floor((x[1]-(b+h/2))/h);
	idx_from[1] = idx - (p_half-1);

	idx = (int) floor(x[2]/h);
	idx_from[2] = idx - (p_half-1);
    }

    return __IDX3_RMAJ(idx_from[0], 
		       idx_from[1], 
		       idx_from[2]+p_half, 
		       params->npdims[1], params->npdims[2]);
}
#endif

// -----------------------------------------------------------------------------
void SE_FGG_expand_all(SE_FGG_work* work, 
		       const SE_state* st, 
		       const SE_FGG_params* params)
{
    double xn[3] MEM_ALIGNED;
    const int N = params->N;
    const int P = params->P;

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int n=0; n<N; n++)
    {
	// compute index and expansion vectors
	xn[0] = st->x[n]; xn[1] = st->x[n+N]; xn[2] = st->x[n+2*N];
	
	*(work->idx+n) = __FGG_EXPA(xn,1,params, 
				    work->zx+n*P, 
				    work->zy+n*P, 
				    work->zz+n*P);
    }
}

// ----------------------------------------------------------------------
// vanilla grid gather
void SE_FGG_int(double* restrict phi,  
		const SE_FGG_work* work, 
		const SE_state* st, 
		const SE_FGG_params* params)
{
    double z2_0[P_MAX] MEM_ALIGNED;
    double z2_1[P_MAX] MEM_ALIGNED;
    double z2_2[P_MAX] MEM_ALIGNED;

    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const int p = params->P;
    const int N = params->N;
    const double h=params->h;

    double xm[3];
    int i,j,k,idx, zidx;
    double phi_m, cij;

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	xm[0] = st->x[m]; xm[1] = st->x[m+N]; xm[2] = st->x[m+2*N];
       
	idx = __FGG_EXPA(xm, 1, params, z2_0, z2_1, z2_2);
	
	phi_m = 0;
	zidx = 0;

	for(i = 0; i<p; i++)
	{
	    for(j = 0; j<p; j++)
	    {
		cij = z2_0[i]*z2_1[j];
		for(k = 0; k<p; k++)
		{
		    phi_m += H[idx]*zs[zidx]*z2_2[k]*cij;
		    idx++; zidx++;
		}
		idx += incrj;
	    }
	    idx += incri;
	}
	phi[m] = (h*h*h)*phi_m;
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_dispatch(double* restrict phi,  
				   const SE_FGG_work* work, 
				   const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#if 0
    // THIS BYPASSES THE FAST SSE KERNELS.
    // 
    // THEY ARE PLATFORM DEPENDENT, AND THUS MAY NOT WORK OUT OF THE BOX.
    // REMOVE THIS BLOCK ONLY IF YOU ARE FAMILIAR WITH BASIC DEBUGGING, 
    // THE BASICS OF SSE INTRINSICS, AND ARE WILLING TO UNDERSTAND WHERE
    // THE (DATA ALIGNMENT) PRECONDITIONS OF SSE INSTRUCTIONS MAY BREAK
    // IN THE SSE CODE BELOW.
    __DISPATCHER_MSG("[FGG INT SSE] SSE Disabled\n");
    SE_FGG_int_split(phi, work, params);
    return;
#endif

    // if P is odd, or if either increment is odd, fall back on vanilla
    if( is_odd(p) || is_odd(incri) || is_odd(incrj) )
    {
	__DISPATCHER_MSG("[FGG INT SSE] SSE Abort (PARAMS)\n");
	SE_FGG_int_split(phi, work, params);
	return;
    }

#if 0
    // If the work arrays zs or zX are misaligned, fall back on vanilla.
    // These arrays are dynamically allocated, so getting this alignment
    // is really the compilers job! Once you trust it, remove this 
    // check, because the long integer modulus operation is not fast.
    if( ( (unsigned long) work->zs)%16 != 0 || 
	( (unsigned long) work->zx)%16 != 0 || 
	( (unsigned long) work->zy)%16 != 0 ||
	( (unsigned long) work->zz)%16 != 0 )
    {
	__DISPATCHER_MSG("[FGG INT SSE] SSE Abort (DATA)\n");
	SE_FGG_int_split(phi, work, params);
	return;
    }
#endif

    // otherwise the preconditions for SSE codes are satisfied. 
    if(p==8)
    {
	// specific for p=8
	__DISPATCHER_MSG("[FGG INT SSE] P=8\n");
	SE_FGG_int_split_SSE_P8(phi, work, params);
    } 
    else if(p==16)
    {
	// specific for p=16
	__DISPATCHER_MSG("[FGG INT SSE] P=16\n");
	SE_FGG_int_split_SSE_P16(phi, work, params); 
    }
    else if(p%8==0)
    {
	// for p divisible by 8
	__DISPATCHER_MSG("[FGG INT SSE] P unroll 8\n");
	SE_FGG_int_split_SSE_u8(phi, work, params); 
    }
    else
    {
	// vanilla SSE code (any even p)
	__DISPATCHER_MSG("[FGG INT SSE] Vanilla\n");
	SE_FGG_int_split_SSE(phi, work, params);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split(double* restrict phi,  
		      const SE_FGG_work* work, 
		      const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    const int p = params->P;
    const int N = params->N;
    const double h=params->h;

    int i,j,k,idx,idx_zs,idx_zz;
    double phi_m, cij;

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];	
	phi_m = 0;
	idx_zs = 0;

	for(i = 0; i<p; i++)
	{
	    for(j = 0; j<p; j++)
	    {
		cij = zx[m*p+i]*zy[m*p+j];
		idx_zz=m*p;
		for(k = 0; k<p; k++)
		{
		    phi_m += H[idx]*zs[idx_zs]*zz[idx_zz]*cij;
		    idx++; idx_zs++; idx_zz++;
		}
		idx += incrj;
	    }
	    idx += incri;
	}
	phi[m] = (h*h*h)*phi_m;
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE(double* restrict phi,  
			  const SE_FGG_work* work, 
			  const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    const int p = params->P;
    const int N = params->N;
    const double h=params->h;

    int i,j,k,idx,idx_zs,idx_zz;
    double s[2] MEM_ALIGNED;

    __m128d rH0, rZZ0, rZS0, rC, rP;

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];	
	idx_zs = 0;
	rP=_mm_setzero_pd();

	if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm_set1_pd( zx[m*p+i]*zy[m*p+j]);
		    idx_zz=m*p;
		    for(k = 0; k<p; k+=2)
		    {
			rH0  = _mm_load_pd( H+idx );
			rZZ0 = _mm_load_pd( zz + idx_zz);
			rZS0 = _mm_load_pd( zs + idx_zs);
			rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
			
			idx+=2; 
			idx_zs+=2; 
			idx_zz+=2;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 16-aligned, so use non-aligned loads
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm_set1_pd( zx[m*p+i]*zy[m*p+j]);
		    idx_zz=m*p;
		    for(k = 0; k<p; k+=2)
		    {
			rH0  = _mm_loadu_pd( H+idx );
			rZZ0 = _mm_load_pd( zz + idx_zz);
			rZS0 = _mm_load_pd( zs + idx_zs);
			rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
			
			idx+=2; 
			idx_zs+=2; 
			idx_zz+=2;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }

	}
	_mm_store_pd(s,rP);
	phi[m] = (h*h*h)*(s[0]+s[1]);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_P8(double* restrict phi,  
			     const SE_FGG_work* work, 
			     const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    /* ASSUME P=8 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;

    int i,j,idx,idx_zs;
    double s[2] MEM_ALIGNED;

    // hold entire zz vector
    __m128d rZZ0, rZZ1, rZZ2, rZZ3; 
    __m128d rC, rP;
    __m128d rH0, rH1, rH2, rH3; 
    __m128d rZS0, rZS1, rZS2, rZS3;

    const int incrj = params->npdims[2]-8;
    const int incri = params->npdims[2]*(params->npdims[1]-8);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];
	idx_zs = 0;
	rP=_mm_setzero_pd();

	/* hoist load of ZZ vector */
	rZZ0 = _mm_load_pd(zz + m*8     );
	rZZ1 = _mm_load_pd(zz + m*8 + 2 );
	rZZ2 = _mm_load_pd(zz + m*8 + 4 );
	rZZ3 = _mm_load_pd(zz + m*8 + 6 );

	if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	    for(i = 0; i<8; i++)
	    {
		for(j = 0; j<8; j++)
		{
		    rC = _mm_set1_pd( zx[m*8+i]*zy[m*8+j]);

		    rH0  = _mm_load_pd( H+idx    );
		    rH1  = _mm_load_pd( H+idx + 2);
		    rH2  = _mm_load_pd( H+idx + 4);
		    rH3  = _mm_load_pd( H+idx + 6);

		    rZS0 = _mm_load_pd( zs + idx_zs    );
		    rZS1 = _mm_load_pd( zs + idx_zs + 2);
		    rZS2 = _mm_load_pd( zs + idx_zs + 4);
		    rZS3 = _mm_load_pd( zs + idx_zs + 6);
		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));

		    idx_zs +=8;
		    idx += incrj + 8;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 16-aligned, so use non-aligned loads
	{
	    for(i = 0; i<8; i++)
	    {
		for(j = 0; j<8; j++)
		{
		    rC = _mm_set1_pd( zx[m*8+i]*zy[m*8+j]);

		    rH0  = _mm_loadu_pd( H+idx    );
		    rH1  = _mm_loadu_pd( H+idx + 2);
		    rH2  = _mm_loadu_pd( H+idx + 4);
		    rH3  = _mm_loadu_pd( H+idx + 6);

		    rZS0 = _mm_load_pd( zs + idx_zs    );
		    rZS1 = _mm_load_pd( zs + idx_zs + 2);
		    rZS2 = _mm_load_pd( zs + idx_zs + 4);
		    rZS3 = _mm_load_pd( zs + idx_zs + 6);
		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));

		    idx_zs +=8;
		    idx += incrj + 8;
		}
		idx += incri;
	    }
	}
	_mm_store_pd(s,rP);
	phi[m] = (h*h*h)*(s[0]+s[1]);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_P16(double* restrict phi,  
			      const SE_FGG_work* work, 
			      const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    /* ASSUME P=16 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;

    int i,j,idx,idx_zs;
    double s[2] MEM_ALIGNED;

    // hold entire zz vector
    __m128d rZZ0, rZZ1, rZZ2, rZZ3, rZZ4, rZZ5, rZZ6, rZZ7; 
    __m128d rC, rP;
    __m128d rH0, rZS0;

    const int incrj = params->npdims[2]-16;
    const int incri = params->npdims[2]*(params->npdims[1]-16);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];
	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);

	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0);

	rP=_mm_setzero_pd();

	/* hoist load of ZZ vector */
	rZZ0 = _mm_load_pd(zz + m*16     );
	rZZ1 = _mm_load_pd(zz + m*16 + 2 );
	rZZ2 = _mm_load_pd(zz + m*16 + 4 );
	rZZ3 = _mm_load_pd(zz + m*16 + 6 );
	rZZ4 = _mm_load_pd(zz + m*16 + 8 );
	rZZ5 = _mm_load_pd(zz + m*16 + 10);
	rZZ6 = _mm_load_pd(zz + m*16 + 12);
	rZZ7 = _mm_load_pd(zz + m*16 + 14);

	if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm_set1_pd( zx[m*16+i]*zy[m*16+j]);

		    /* 0 */ 
		    rH0  = _mm_load_pd( H+idx );
		    rZS0 = _mm_load_pd( zs + idx_zs);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));

		    /* 1 */ 
		    rH0  = _mm_load_pd( H+idx + 2);
		    rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0)));

		    /* 2 */ 
		    rH0  = _mm_load_pd( H+idx + 4);
		    rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0)));

		    /* 3 */ 
		    rH0  = _mm_load_pd( H+idx + 6);
		    rZS0 = _mm_load_pd( zs + idx_zs + 6);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0)));

		    /* 4 */ 
		    rH0  = _mm_load_pd( H+idx + 8);
		    rZS0 = _mm_load_pd( zs + idx_zs + 8);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0)));

		    /* 5 */ 
		    rH0  = _mm_load_pd( H+idx + 10);
		    rZS0 = _mm_load_pd( zs + idx_zs + 10);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0)));

		    /* 6 */ 
		    rH0  = _mm_load_pd( H+idx + 12);
		    rZS0 = _mm_load_pd( zs + idx_zs + 12);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0)));

		    /* 7 */ 
		    rH0  = _mm_load_pd( H+idx + 14);
		    rZS0 = _mm_load_pd( zs + idx_zs + 14);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0)));

		    idx_zs +=16;
		    idx += incrj + 16;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 16-aligned, so use non-aligned loads
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm_set1_pd( zx[m*16+i]*zy[m*16+j]);

		    /* 0 */ 
		    rH0  = _mm_loadu_pd( H+idx );
		    rZS0 = _mm_load_pd( zs + idx_zs);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));

		    /* 1 */ 
		    rH0  = _mm_loadu_pd( H+idx + 2);
		    rZS0 = _mm_load_pd( zs + idx_zs + 2);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0)));

		    /* 2 */ 
		    rH0  = _mm_loadu_pd( H+idx + 4);
		    rZS0 = _mm_load_pd( zs + idx_zs + 4);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0)));

		    /* 3 */ 
		    rH0  = _mm_loadu_pd( H+idx + 6);
		    rZS0 = _mm_load_pd( zs + idx_zs + 6);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0)));

		    /* 4 */ 
		    rH0  = _mm_loadu_pd( H+idx + 8);
		    rZS0 = _mm_load_pd( zs + idx_zs + 8);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0)));

		    /* 5 */ 
		    rH0  = _mm_loadu_pd( H+idx + 10);
		    rZS0 = _mm_load_pd( zs + idx_zs + 10);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0)));

		    /* 6 */ 
		    rH0  = _mm_loadu_pd( H+idx + 12);
		    rZS0 = _mm_load_pd( zs + idx_zs + 12);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0)));

		    /* 7 */ 
		    rH0  = _mm_loadu_pd( H+idx + 14);
		    rZS0 = _mm_load_pd( zs + idx_zs + 14);		    
		    rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0)));

		    idx_zs +=16;
		    idx += incrj + 16;
		}
		idx += incri;
	    }
	}
	_mm_store_pd(s,rP);
	phi[m] = (h*h*h)*(s[0]+s[1]);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_SSE_u8(double* restrict phi,  
			     const SE_FGG_work* work, 
			     const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    const int p = params->P;
    const int N = params->N;
    const double h=params->h;

    int i,j,k,idx,idx_zs,idx_zz;
    double s[2] MEM_ALIGNED;

    __m128d rH0, rZZ0, rZS0, rC, rP;
    __m128d rH1, rZZ1, rZS1;
    __m128d rH2, rZZ2, rZS2;
    __m128d rH3, rZZ3, rZS3;

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];
	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
	
	idx_zs = 0;
	_mm_prefetch( (void*) zs, _MM_HINT_T0);

	rP=_mm_setzero_pd();

	if(idx%2==0) // H[idx] is 16-aligned so vectorization simple
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm_set1_pd( zx[m*p+i]*zy[m*p+j] );
		    idx_zz=m*p;

		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm_load_pd( H+idx    );
			rH1  = _mm_load_pd( H+idx + 2);
			rH2  = _mm_load_pd( H+idx + 4);
			rH3  = _mm_load_pd( H+idx + 6);

			rZZ0 = _mm_load_pd( zz + idx_zz    );
			rZZ1 = _mm_load_pd( zz + idx_zz + 2);
			rZZ2 = _mm_load_pd( zz + idx_zz + 4);
			rZZ3 = _mm_load_pd( zz + idx_zz + 6);

			rZS0 = _mm_load_pd( zs + idx_zs    );
			rZS1 = _mm_load_pd( zs + idx_zs + 2);
			rZS2 = _mm_load_pd( zs + idx_zs + 4);
			rZS3 = _mm_load_pd( zs + idx_zs + 6);

			rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));
			
			idx+=8; 
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 16-aligned, so use non-aligned load from H
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm_set1_pd( zx[m*p+i]*zy[m*p+j] );
		    idx_zz=m*p;

		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm_loadu_pd( H+idx    );
			rH1  = _mm_loadu_pd( H+idx + 2);
			rH2  = _mm_loadu_pd( H+idx + 4);
			rH3  = _mm_loadu_pd( H+idx + 6);

			rZZ0 = _mm_load_pd( zz + idx_zz    );
			rZZ1 = _mm_load_pd( zz + idx_zz + 2);
			rZZ2 = _mm_load_pd( zz + idx_zz + 4);
			rZZ3 = _mm_load_pd( zz + idx_zz + 6);

			rZS0 = _mm_load_pd( zs + idx_zs    );
			rZS1 = _mm_load_pd( zs + idx_zs + 2);
			rZS2 = _mm_load_pd( zs + idx_zs + 4);
			rZS3 = _mm_load_pd( zs + idx_zs + 6);

			rP = _mm_add_pd(rP,_mm_mul_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2)));
			rP = _mm_add_pd(rP,_mm_mul_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3)));
			
			idx+=8; 
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}

	// done accumulating
	_mm_store_pd(s,rP);
	phi[m] = (h*h*h)*(s[0]+s[1]);
    }
}

// -----------------------------------------------------------------------------
#ifdef __AVX__
void SE_FGG_int_split_AVX_dispatch(double* restrict phi,
                                   const SE_FGG_work* work,
                                   const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#ifdef AVX_FMA
    __DISPATCHER_MSG("[FGG INT AVX-FMA] ");
#else
    __DISPATCHER_MSG("[FGG INT AVX] ");
#endif
#if 0
    // THIS BYPASSES THE FAST AVX KERNELS.
    __DISPATCHER_MSG("AVX Disabled\n");
    SE_FGG_int_split(phi, work, params);
    return;
#endif
    
    // if either P or increments are not divisible by 4, fall back to SSE
    if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) )
    {
        __DISPATCHER_MSG("AVX Abort (PARAMS)\n");
        SE_FGG_int_split_SSE_dispatch(phi, work, params);
        return;
    }

    // otherwise the preconditions for AVX codes are satisfied. 
    if(p==16)
    {
        // specific for p=16
        __DISPATCHER_MSG("P=16\n");
        SE_FGG_int_split_AVX_P16(phi, work, params);
    }
    else if(p==8)
    {
        // specific for p=8
        __DISPATCHER_MSG("P=8\n");
        SE_FGG_int_split_AVX_P8(phi, work, params);
    }
    else if(p%8==0)
    {
      // specific for p divisible by 8
      __DISPATCHER_MSG("P unroll 8\n");
      SE_FGG_int_split_AVX_u8(phi, work, params);
    }
    else if(p%4==0)
    {
        // specific for p divisible by 4
        __DISPATCHER_MSG("P unroll 4\n");
        SE_FGG_int_split_AVX(phi, work, params);
    }
}


// -----------------------------------------------------------------------------
void SE_FGG_int_split_AVX(double* restrict phi,
                          const SE_FGG_work* work,
                          const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    const int p = params->P;
    const int N = params->N;
    const double h=params->h;

    int i,j,k,idx,idx_zs,idx_zz;
    double s[4] MEM_ALIGNED;

    __m256d rH0, rZZ0, rZS0, rC, rP;

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
        idx = work->idx[m];
        idx_zs = 0;
        rP=_mm256_setzero_pd();

        if(idx%4==0) // H[idx] is 32-aligned so vectorization simple
        {
            for(i = 0; i<p; i++)
            {
                for(j = 0; j<p; j++)
                {
                    rC = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j]);
                    idx_zz=m*p;
                    for(k = 0; k<p; k+=4)
                    {
                        rH0  = _mm256_load_pd( H+idx );
                        rZZ0 = _mm256_load_pd( zz + idx_zz);
                        rZS0 = _mm256_load_pd( zs + idx_zs);
#ifdef AVX_FMA
                        rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0),rP);
#else
                        rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
#endif
                        idx+=4;
                        idx_zs+=4;
                        idx_zz+=4;
                    }
                    idx += incrj;
                }
                idx += incri;
            }
        }
        else // H[idx] not 32-aligned, so use non-aligned loads
        {
            for(i = 0; i<p; i++)
            {
                for(j = 0; j<p; j++)
                {
                    rC = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j]);
                    idx_zz=m*p;
                    for(k = 0; k<p; k+=4)
                    {
                        rH0  = _mm256_loadu_pd( H+idx );
                        rZZ0 = _mm256_load_pd( zz + idx_zz);
                        rZS0 = _mm256_load_pd( zs + idx_zs);
#ifdef AVX_FMA
                        rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0),rP);
#else
                        rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
#endif
                        idx+=4;
                        idx_zs+=4;
                        idx_zz+=4;
                    }
                    idx += incrj;
                }
                idx += incri;
            }

        }
        _mm256_store_pd(s,rP);
        phi[m] += (h*h*h)*(s[0]+s[1]+s[2]+s[3]);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_AVX_P8(double* restrict phi,
                             const SE_FGG_work* work,
                             const SE_FGG_params* params)
{  
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    // ASSUME P=8 const int p = params->P; //
    const int N = params->N;
    const double h=params->h;

    int i,j,idx,idx_zs;
    double s[4] MEM_ALIGNED;
    
    // hold entire zz vector
    __m256d rZZ0, rZZ1;
    __m256d rC, rP;
    __m256d rH0, rH1;
    __m256d rZS0, rZS1;
    
    const int incrj = params->npdims[2]-8;
    const int incri = params->npdims[2]*(params->npdims[1]-8);
    
#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
        idx = work->idx[m];
        idx_zs = 0;
        rP = _mm256_setzero_pd();

        // hoist load of ZZ vector //
        rZZ0 = _mm256_load_pd(zz + m*8     );
        rZZ1 = _mm256_load_pd(zz + m*8 + 4 );

        if(idx%4==0) // H[idx] is 32-aligned so vectorization simple
        {
            for(i = 0; i<8; i++)
            {
                for(j = 0; j<8; j++)
                {
                    rC = _mm256_set1_pd( zx[m*8+i]*zy[m*8+j]);

                    rH0  = _mm256_load_pd( H+idx    );
                    rH1  = _mm256_load_pd( H+idx + 4);

                    rZS0 = _mm256_load_pd( zs + idx_zs    );
                    rZS1 = _mm256_load_pd( zs + idx_zs + 4);
#ifdef AVX_FMA
                    rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0),rP);
                    rP = _mm256_fmadd_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1),rP);
#else
                    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
                    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
#endif
                    idx_zs +=8;
                    idx += incrj + 8;
                }
                idx += incri;
            }
        }
        else // H[idx] not 32-aligned, so use non-aligned loads
        {
            for(i = 0; i<8; i++)
            {
                for(j = 0; j<8; j++)
                {
                    rC = _mm256_set1_pd( zx[m*8+i]*zy[m*8+j]);

                    rH0  = _mm256_loadu_pd( H+idx    );
                    rH1  = _mm256_loadu_pd( H+idx + 4);

                    rZS0 = _mm256_load_pd( zs + idx_zs    );
                    rZS1 = _mm256_load_pd( zs + idx_zs + 4);
#ifdef AVX_FMA
                    rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0),rP);
                    rP = _mm256_fmadd_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1),rP);
#else
                    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
                    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
#endif
                    idx_zs +=8;
                    idx += incrj + 8;
                }
                idx += incri;
            }
        }
        _mm256_store_pd(s,rP);
        _mm256_store_pd(s,rP);
	phi[m] = (h*h*h)*(s[0]+s[1]+s[2]+s[3]);
    }
 
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_AVX_P16(double* restrict phi,  
			     const SE_FGG_work* work, 
			     const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    /* ASSUME P=16 const int p = params->P; */
    const int N = params->N;
    const double h=params->h;

    int i,j,idx,idx_zs;
    double s[4] MEM_ALIGNED;

    // hold entire zz vector
    __m256d rZZ0, rZZ1, rZZ2, rZZ3; 
    __m256d rC, rP;
    __m256d rH0, rH1, rH2, rH3; 
    __m256d rZS0, rZS1, rZS2, rZS3;

    const int incrj = params->npdims[2]-16;
    const int incri = params->npdims[2]*(params->npdims[1]-16);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];
	idx_zs = 0;
	rP=_mm256_setzero_pd();

	/* hoist load of ZZ vector */
	rZZ0 = _mm256_load_pd(zz + m*16     );
	rZZ1 = _mm256_load_pd(zz + m*16 + 4 );
	rZZ2 = _mm256_load_pd(zz + m*16 + 8 );
	rZZ3 = _mm256_load_pd(zz + m*16 + 12);

	if(idx%4==0) // H[idx] is 32-aligned so vectorization simple
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm256_set1_pd( zx[m*16+i]*zy[m*16+j]);

		    rH0  = _mm256_load_pd( H+idx     );
		    rH1  = _mm256_load_pd( H+idx + 4 );
		    rH2  = _mm256_load_pd( H+idx + 8 );
		    rH3  = _mm256_load_pd( H+idx + 12);

		    rZS0 = _mm256_load_pd( zs + idx_zs     );
		    rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
		    rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
		    rZS3 = _mm256_load_pd( zs + idx_zs + 12);
#ifdef AVX_FMA
		    rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0),rP);
		    rP = _mm256_fmadd_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1),rP);
		    rP = _mm256_fmadd_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2),rP);
		    rP = _mm256_fmadd_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3),rP);
#else		    
		    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
		    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
		    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2)));
		    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3)));
#endif
		    idx_zs +=16;
		    idx += incrj + 16;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 32-aligned, so use non-aligned loads
	{
	    for(i = 0; i<16; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm256_set1_pd( zx[m*16+i]*zy[m*16+j]);

		    rH0  = _mm256_loadu_pd( H+idx     );
		    rH1  = _mm256_loadu_pd( H+idx + 4 );
		    rH2  = _mm256_loadu_pd( H+idx + 8 );
		    rH3  = _mm256_loadu_pd( H+idx + 12);

		    rZS0 = _mm256_load_pd( zs + idx_zs     );
		    rZS1 = _mm256_load_pd( zs + idx_zs + 4 );
		    rZS2 = _mm256_load_pd( zs + idx_zs + 8 );
		    rZS3 = _mm256_load_pd( zs + idx_zs + 12);
#ifdef AVX_FMA
		    rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0),rP);
		    rP = _mm256_fmadd_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1),rP);
		    rP = _mm256_fmadd_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2),rP);
		    rP = _mm256_fmadd_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3),rP);
#else		    
		    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
		    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
		    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2)));
		    rP = _mm256_add_pd(rP,_mm256_mul_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3)));
#endif
		    idx_zs +=16;
		    idx += incrj + 16;
		}
		idx += incri;
	    }
	}
	_mm256_store_pd(s,rP);
	phi[m] = (h*h*h)*(s[0]+s[1]+s[2]+s[3]);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_int_split_AVX_u8(double* restrict phi,  
			     const SE_FGG_work* work, 
			     const SE_FGG_params* params)
{
    // unpack params
    const double* restrict H = work->H;
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    const int p = params->P;
    const int N = params->N;
    const double h=params->h;

    int i,j,k,idx,idx_zs,idx_zz;
    double s[4] MEM_ALIGNED;

    __m256d rH0, rZZ0, rZS0, rC, rP;
    __m256d rH1, rZZ1, rZS1;

    const int incrj = params->npdims[2]-p;
    const int incri = params->npdims[2]*(params->npdims[1]-p);

#ifdef _OPENMP
#pragma omp for // work-share over OpenMP threads here
#endif
    for(int m=0; m<N; m++)
    {
	idx = work->idx[m];
	idx_zs = 0;
	rP=_mm256_setzero_pd();

	if(idx%4==0) // H[idx] is 32-aligned so vectorization simple
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j] );
		    idx_zz=m*p;

		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm256_load_pd( H+idx    );
			rH1  = _mm256_load_pd( H+idx + 4);

			rZZ0 = _mm256_load_pd( zz + idx_zz    );
			rZZ1 = _mm256_load_pd( zz + idx_zz + 4);

			rZS0 = _mm256_load_pd( zs + idx_zs    );
			rZS1 = _mm256_load_pd( zs + idx_zs + 4);
#ifdef AVX_FMA
			rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0),rP);
			rP = _mm256_fmadd_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1),rP);
#else
			rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
			rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
#endif			
			idx+=8; 
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}
	else // H[idx] not 32-aligned, so use non-aligned load from H
	{
	    for(i = 0; i<p; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm256_set1_pd( zx[m*p+i]*zy[m*p+j] );
		    idx_zz=m*p;

		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm256_loadu_pd( H+idx    );
			rH1  = _mm256_loadu_pd( H+idx + 4);

			rZZ0 = _mm256_load_pd( zz + idx_zz    );
			rZZ1 = _mm256_load_pd( zz + idx_zz + 4);

			rZS0 = _mm256_load_pd( zs + idx_zs    );
			rZS1 = _mm256_load_pd( zs + idx_zs + 4);
#ifdef AVX_FMA
			rP = _mm256_fmadd_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0),rP);
			rP = _mm256_fmadd_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1),rP);
#else
			rP = _mm256_add_pd(rP,_mm256_mul_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0)));
			rP = _mm256_add_pd(rP,_mm256_mul_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1)));
#endif			
			idx+=8; 
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx += incrj;
		}
		idx += incri;
	    }
	}

	// done accumulating
	_mm256_store_pd(s,rP);
	phi[m] = (h*h*h)*(s[0]+s[1]+s[2]+s[3]);
    }
}
#endif // AVX

// -----------------------------------------------------------------------------
void SE_FGG_grid(SE_FGG_work* work, const SE_state* st, 
		 const SE_FGG_params* params)
{
    // vectors for FGG expansions
    double zx0[P_MAX] MEM_ALIGNED;
    double zy0[P_MAX] MEM_ALIGNED;
    double zz0[P_MAX] MEM_ALIGNED;

    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const int p = params->P;
    
    double cij0;
    double xn[3];
    double qn;
    int idx0, zidx, i,j,k, i_end;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, p);

    for(int n=0; n<N; n++)
    {
	// compute index and expansion vectors
	xn[0] = st->x[n]; xn[1] = st->x[n+N]; xn[2] = st->x[n+2*N];
	idx0 = __FGG_INDEX(xn, params);
	grid_thrd_slice(&grid_thrd_ws, &idx0, &i, &i_end, &zidx);
	if (grid_thrd_ws.skip) 
	    continue;
	qn = st->q[n];		
	__FGG_EXPA(xn, qn, params, zx0, zy0, zz0);
	for(; i < i_end; i++)
	{
	    for(j = 0; j<p; j++)
	    {
		cij0 = zx0[i]*zy0[j];
		for(k = 0; k<p; k++)
		{
		    H[idx0] += zs[zidx]*zz0[k]*cij0;
		    idx0++; 
		    zidx++;
		}
		idx0 += incrj; 
	    }
	    idx0 += incri; 
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE_dispatch(SE_FGG_work* work, const SE_state* st, 
				    const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#if 0
    // THIS BYPASSES THE FAST SSE KERNELS.
    // 
    // THEY ARE PLATFORM DEPENDENT, AND THUS MAY NOT WORK OUT OF THE BOX.
    // REMOVE THIS BLOCK ONLY IF YOU ARE FAMILIAR WITH BASIC DEBUGGING, 
    // THE BASICS OF SSE INTRINSICS, AND ARE WILLING TO UNDERSTAND WHERE
    // THE (DATA ALIGNMENT) PRECONDITIONS OF SSE INSTRUCTIONS MAY BREAK
    // IN THE SSE CODE BELOW.
    __DISPATCHER_MSG("[FGG GRID SSE] SSE Disabled\n");
    SE_FGG_grid_split(work, st, params);
    return;
#endif

    // if P is odd, or if either increment is odd, fall back on vanilla
    if( is_odd(p) || is_odd(incri) || is_odd(incrj) )
    {
	__DISPATCHER_MSG("[FGG GRID SSE] SSE Abort (PARAMS)\n");
	SE_FGG_grid_split(work, st, params);
	return;
    }

#if 0
    // If the work arrays zs or zX are misaligned, fall back on vanilla.
    // These arrays are dynamically allocated, so getting this alignment
    // is really the compilers job! Once you trust it, remove this 
    // check, because the long integer modulus operation is not fast.
    if( ( (unsigned long) work->zs)%16 != 0 || 
	( (unsigned long) work->zx)%16 != 0 || 
	( (unsigned long) work->zy)%16 != 0 ||
	( (unsigned long) work->zz)%16 != 0 )
    {
	__DISPATCHER_MSG("[FGG GRID SSE] SSE Abort (DATA)\n");
	SE_FGG_grid_split(work, st, params);
	return;
    }
#endif

    // otherwise the preconditions for SSE codes are satisfied. 
    if(p==16)
    {
	// specific for p=16
	__DISPATCHER_MSG("[FGG GRID SSE] P=16\n");
	SE_FGG_grid_split_SSE_P16(work, st, params); 
    }
    else if(p%8==0)
    {
	// specific for p divisible by 8
	__DISPATCHER_MSG("[FGG GRID SSE] P unroll 8\n");
	SE_FGG_grid_split_SSE_u8(work, st, params); 
    }
    else
    {
	// vanilla SSE code (any even p)
	__DISPATCHER_MSG("[FGG GRID SSE] Vanilla\n");
	SE_FGG_grid_split_SSE(work, st, params);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split(SE_FGG_work* work, const SE_state* st, 
		       const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    double cij0;
    double qn;
    int idx0, zidx, idxzz, i, j, k, i_end;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment
    
    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, p);

    for(int n=0; n<N; n++)
    {
	idx0 = work->idx[n];
	grid_thrd_slice(&grid_thrd_ws, &idx0, &i, &i_end, &zidx);
	if (grid_thrd_ws.skip) 
	    continue;
	qn = st->q[n];
	// inline vanilla loop
	for(; i < i_end; i++)
	{
	    for(j = 0; j<p; j++)
	    {
		cij0 = qn*zx[p*n+i]*zy[p*n+j];
		idxzz=p*n;
		for(k = 0; k<p; k++)
		{
		    H[idx0] += zs[zidx]*zz[idxzz]*cij0;
		    idx0++; zidx++; idxzz++;
		}
		idx0 += incrj; 
	    }
	    idx0 += incri; 
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE_P16(SE_FGG_work* work, const SE_state* st, 
			       const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    double qn;
    int idx, idx_zs, i, j;
    const int incrj = params->npdims[2]-16; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-16);// outer increment

    __m128d rZZ0, rZZ1, rZZ2, rZZ3, rZZ4, rZZ5, rZZ6, rZZ7; 
    __m128d rH0, rH1, rH2, rH3;
    __m128d rC, rZS0;

    int i_end;
    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, 16);

    for(int n=0; n<N; n++)
    {
	idx = work->idx[n];
	idx_zs = 0;
	grid_thrd_slice(&grid_thrd_ws, &idx, &i, &i_end, &idx_zs);
	if (grid_thrd_ws.skip) 
	    continue;

	_mm_prefetch( (void*) (H+idx), _MM_HINT_T0);
	_mm_prefetch( (void*) zs, _MM_HINT_T0);
	qn = st->q[n];

        rZZ0 = _mm_load_pd(zz + n*16     );
        rZZ1 = _mm_load_pd(zz + n*16 + 2 );
        rZZ2 = _mm_load_pd(zz + n*16 + 4 );
        rZZ3 = _mm_load_pd(zz + n*16 + 6 );
        rZZ4 = _mm_load_pd(zz + n*16 + 8 );
        rZZ5 = _mm_load_pd(zz + n*16 + 10);
        rZZ6 = _mm_load_pd(zz + n*16 + 12);
        rZZ7 = _mm_load_pd(zz + n*16 + 14);

	if(idx%2 == 0) // H[idx0] is 16-aligned
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

                    /* 0 - 3 */ 
                    rH0  = _mm_load_pd( H+idx    );
                    rH1  = _mm_load_pd( H+idx + 2);
                    rH2  = _mm_load_pd( H+idx + 4);
                    rH3  = _mm_load_pd( H+idx + 6);

                    rZS0 = _mm_load_pd( zs + idx_zs);
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 2);                   
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 4);                   
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 6);                   
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0));

		    _mm_store_pd(H + idx, rH0);
		    _mm_store_pd(H + idx + 2, rH1);
		    _mm_store_pd(H + idx + 4, rH2);
		    _mm_store_pd(H + idx + 6, rH3);

                    /* 4 - 7*/ 
		    rH0  = _mm_load_pd( H+idx + 8 );
                    rH1  = _mm_load_pd( H+idx + 10);
                    rH2  = _mm_load_pd( H+idx + 12);
                    rH3  = _mm_load_pd( H+idx + 14);

                    rZS0 = _mm_load_pd( zs + idx_zs + 8);
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 10);                   
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 12);                   
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 14);                   
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0));

		    _mm_store_pd(H + idx + 8 , rH0);
		    _mm_store_pd(H + idx + 10, rH1);
		    _mm_store_pd(H + idx + 12, rH2);
		    _mm_store_pd(H + idx + 14, rH3);

		    idx += incrj + 16;
		    idx_zs += 16;
		}
		idx += incri;
	    }
	}
	else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

                    /* 0 - 3 */ 
                    rH0  = _mm_loadu_pd( H+idx    );
                    rH1  = _mm_loadu_pd( H+idx + 2);
                    rH2  = _mm_loadu_pd( H+idx + 4);
                    rH3  = _mm_loadu_pd( H+idx + 6);

		    // if zs does not have 16-byte alignment, this will core.
		    // PLATFORM AND COMPILER DEPENDENT (FIXME)
                    rZS0 = _mm_load_pd( zs + idx_zs);
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 2);                   
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 4);                   
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 6);                   
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS0));

		    _mm_storeu_pd(H + idx, rH0);
		    _mm_storeu_pd(H + idx + 2, rH1);
		    _mm_storeu_pd(H + idx + 4, rH2);
		    _mm_storeu_pd(H + idx + 6, rH3);

                    /* 4 - 7*/ 
		    rH0  = _mm_loadu_pd( H+idx + 8 );
                    rH1  = _mm_loadu_pd( H+idx + 10);
                    rH2  = _mm_loadu_pd( H+idx + 12);
                    rH3  = _mm_loadu_pd( H+idx + 14);

                    rZS0 = _mm_load_pd( zs + idx_zs + 8);
		    rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ4,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 10);                   
		    rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ5,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 12);                   
		    rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ6,rC),rZS0));

                    rZS0 = _mm_load_pd( zs + idx_zs + 14);                   
		    rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ7,rC),rZS0));

		    _mm_storeu_pd(H + idx + 8 , rH0);
		    _mm_storeu_pd(H + idx + 10, rH1);
		    _mm_storeu_pd(H + idx + 12, rH2);
		    _mm_storeu_pd(H + idx + 14, rH3);

		    idx += incrj + 16;
		    idx_zs += 16;
		}
		idx += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE_u8(SE_FGG_work* work, const SE_state* st, 
			      const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    double qn;
    int idx0, idx_zs, idx_zz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    __m128d rH0, rZZ0, rZS0, rC;
    __m128d rH1, rZZ1, rZS1;
    __m128d rH2, rZZ2, rZS2;
    __m128d rH3, rZZ3, rZS3;

    int i_end;
    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, p);

    for(int n=0; n<N; n++)
    {
	idx0 = work->idx[n];
	idx_zs = 0;
	grid_thrd_slice(&grid_thrd_ws, &idx0, &i, &i_end, &idx_zs);
	if (grid_thrd_ws.skip) 
	    continue;

	qn = st->q[n];
	_mm_prefetch( (void*) (H+idx0), _MM_HINT_T0);
	_mm_prefetch( (void*) zs, _MM_HINT_T0);

	if(idx0%2 == 0) // H[idx0] is 16-aligned
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		    idx_zz=p*n;

		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm_load_pd( H+idx0     );
			rH1  = _mm_load_pd( H+idx0 + 2 );
			rH2  = _mm_load_pd( H+idx0 + 4 );
			rH3  = _mm_load_pd( H+idx0 + 6 );

			rZZ0 = _mm_load_pd( zz + idx_zz     );
			rZZ1 = _mm_load_pd( zz + idx_zz + 2 );
			rZZ2 = _mm_load_pd( zz + idx_zz + 4 );
			rZZ3 = _mm_load_pd( zz + idx_zz + 6 );

			rZS0 = _mm_load_pd( zs + idx_zs    );
			rZS1 = _mm_load_pd( zs + idx_zs + 2);
			rZS2 = _mm_load_pd( zs + idx_zs + 4);
			rZS3 = _mm_load_pd( zs + idx_zs + 6);

			rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));
			rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1));
			rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2));
			rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3));
			
			_mm_store_pd( H+idx0    , rH0 );
			_mm_store_pd( H+idx0 + 2, rH1 );
			_mm_store_pd( H+idx0 + 4, rH2 );
			_mm_store_pd( H+idx0 + 6, rH3 );

			idx0  +=8;
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx0 += incrj;
		}
		idx0 += incri;
	    }
	}
	else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	    for(; i<i_end; i++)
	    {
	    	for(j = 0; j<p; j++)
	    	{
	    	    rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	    	    idx_zz=p*n;
		    
	    	    for(k = 0; k<p; k+=8)
	    	    {
	    		rH0  = _mm_loadu_pd( H+idx0     );
	    		rH1  = _mm_loadu_pd( H+idx0 + 2 );
	    		rH2  = _mm_loadu_pd( H+idx0 + 4 );
	    		rH3  = _mm_loadu_pd( H+idx0 + 6 );

	    		rZZ0 = _mm_load_pd( zz + idx_zz     );
	    		rZZ1 = _mm_load_pd( zz + idx_zz + 2 );
	    		rZZ2 = _mm_load_pd( zz + idx_zz + 4 );
	    		rZZ3 = _mm_load_pd( zz + idx_zz + 6 );

	    		rZS0 = _mm_load_pd( zs + idx_zs    );
	    		rZS1 = _mm_load_pd( zs + idx_zs + 2);
	    		rZS2 = _mm_load_pd( zs + idx_zs + 4);
	    		rZS3 = _mm_load_pd( zs + idx_zs + 6);

			rH0 = _mm_add_pd(rH0,_mm_mul_pd(_mm_mul_pd(rZZ0,rC),rZS0));
			rH1 = _mm_add_pd(rH1,_mm_mul_pd(_mm_mul_pd(rZZ1,rC),rZS1));
			rH2 = _mm_add_pd(rH2,_mm_mul_pd(_mm_mul_pd(rZZ2,rC),rZS2));
			rH3 = _mm_add_pd(rH3,_mm_mul_pd(_mm_mul_pd(rZZ3,rC),rZS3));

	    		_mm_storeu_pd( H+idx0    , rH0 );
	    		_mm_storeu_pd( H+idx0 + 2, rH1 );
	    		_mm_storeu_pd( H+idx0 + 4, rH2 );
	    		_mm_storeu_pd( H+idx0 + 6, rH3 );

	    		idx0  +=8;
	    		idx_zs+=8;
	    		idx_zz+=8;
	    	    }
	    	    idx0 += incrj;
	    	}
	    	idx0 += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_SSE(SE_FGG_work* work, const SE_state* st, 
			   const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    double qn;
    int idx0, idx_zs, idx_zz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    __m128d rH0, rZZ0, rZS0, rC;

    int i_end;
    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, p);

    for(int n=0; n<N; n++)
    {
	idx0 = work->idx[n];
	idx_zs = 0;

	grid_thrd_slice(&grid_thrd_ws, &idx0, &i, &i_end, &idx_zs);
	if (grid_thrd_ws.skip) 
	    continue;

	qn = st->q[n];
	if(idx0%2 == 0) // H[idx0] is 16-aligned
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		    idx_zz=p*n;
		    for(k = 0; k<p; k+=2)
		    {
			rH0  = _mm_load_pd( H+idx0     );
			rZZ0 = _mm_load_pd( zz + idx_zz     );
			rZS0 = _mm_load_pd( zs + idx_zs    );

			rZZ0 = _mm_mul_pd(rZZ0,rC);
			rZZ0 = _mm_mul_pd(rZZ0,rZS0);
			rH0  = _mm_add_pd(rH0,rZZ0);

			_mm_store_pd( H+idx0    , rH0 );

			idx0  +=2;
			idx_zs+=2; 
			idx_zz+=2;
		    }
		    idx0 += incrj; 
		}
		idx0 += incri; 
	    }
	}
	else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	    for(; i<i_end; i++)
	    {
	    	for(j = 0; j<p; j++)
	    	{
	    	    rC = _mm_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	    	    idx_zz=p*n;
	    	    for(k = 0; k<p; k+=2)
	    	    {
	    		rH0  = _mm_loadu_pd( H+idx0 );
	    		rZZ0 = _mm_load_pd( zz + idx_zz );
	    		rZS0 = _mm_load_pd( zs + idx_zs );
	    		rZZ0 = _mm_mul_pd(rZZ0,rC);
	    		rZZ0 = _mm_mul_pd(rZZ0,rZS0);
	    		rH0  = _mm_add_pd(rH0,rZZ0);
	    		_mm_storeu_pd( H+idx0, rH0 );

	    		idx0  +=2;
	    		idx_zs+=2;
	    		idx_zz+=2;
	    	    }
	    	    idx0 += incrj;
	    	}
	    	idx0 += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
#ifdef __AVX__
void SE_FGG_grid_split_AVX_dispatch(SE_FGG_work* work, const SE_state* st, 
				    const SE_FGG_params* params)
{
    const int p = params->P;
    const int incrj = params->dims[2]; // middle increment
    const int incri = params->npdims[2]*(params->dims[1]);// outer increment

#ifdef AVX_FMA
    __DISPATCHER_MSG("[FGG GRID AVX-FMA");
#else
    __DISPATCHER_MSG("[FGG GRID AVX");
#endif
#ifdef FGG_THRD
    __DISPATCHER_MSG(" THRD] ");
#else
    __DISPATCHER_MSG("] ");
#endif

#if 0
    // THIS BYPASSES THE FAST AVX KERNELS.
    // 
    // THEY ARE PLATFORM DEPENDENT, AND THUS MAY NOT WORK OUT OF THE BOX.
    // REMOVE THIS BLOCK ONLY IF YOU ARE FAMILIAR WITH BASIC DEBUGGING, 
    // THE BASICS OF AVX INTRINSICS, AND ARE WILLING TO UNDERSTAND WHERE
    // THE (DATA ALIGNMENT) PRECONDITIONS OF AVX INSTRUCTIONS MAY BREAK
    // IN THE AVX CODE BELOW.
    __DISPATCHER_MSG("AVX Disabled\n");
    SE_FGG_grid_split(work, st, params);
    return;
#endif

    // if either P or increments are not divisible by 4, fall back to SSE
    if( isnot_div_by_4(p) || isnot_div_by_4(incri) || isnot_div_by_4(incrj) )
    {
	__DISPATCHER_MSG("AVX Abort (PARAMS)\n");
	SE_FGG_grid_split_SSE_dispatch(work, st, params);
	return;
    }

#if 0
    // If the work arrays zs or zX are misaligned, fall back on vanilla.
    // These arrays are dynamically allocated, so getting this alignment
    // is really the compilers job! Once you trust it, remove this 
    // check, because the long integer modulus operation is not fast.
    if( ( (unsigned long) work->zs)%32 != 0 || 
	( (unsigned long) work->zx)%32 != 0 || 
	( (unsigned long) work->zy)%32 != 0 ||
	( (unsigned long) work->zz)%32 != 0 )
    {
	__DISPATCHER_MSG("AVX Abort (DATA)\n");
	SE_FGG_grid_split(work, st, params);
	return;
    }
#endif

    // otherwise the preconditions for AVX codes are satisfied. 
    if(p==16)
    {
	// specific for p=16
	__DISPATCHER_MSG("P=16\n");
	SE_FGG_grid_split_AVX_P16(work, st, params); 
    }
    else if(p==8)
    {
      // specific for p=8
      __DISPATCHER_MSG("P=8\n");
      SE_FGG_grid_split_AVX_P8(work, st, params);
    }
    else if(p%8==0)
    {
	// specific for p divisible by 8
	__DISPATCHER_MSG("P unroll 8\n");
	SE_FGG_grid_split_AVX_u8(work, st, params); 
    }
    else if(p%4==0)
    {
      // specific for p divisible by 4
      __DISPATCHER_MSG("P unroll 4\n");
      SE_FGG_grid_split_AVX(work, st, params);
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_AVX(SE_FGG_work* work, const SE_state* st, 
			   const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    double qn;
    int idx0, idx_zs, idx_zz, i, j, k;

    __m256d rH0, rZZ0, rZS0, rC;

    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    int i_end;
    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, p);

    for(int n=0; n<N; n++)
    {
	idx0 = work->idx[n];
	idx_zs = 0;

	grid_thrd_slice(&grid_thrd_ws, &idx0, &i, &i_end, &idx_zs);
	if (grid_thrd_ws.skip) 
	    continue;

	qn = st->q[n];

	if(idx0%4 == 0) // H[idx0] is 16-aligned
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		    idx_zz=p*n;
		    for(k = 0; k<p; k+=4)
		    {
			rH0  = _mm256_load_pd( H+idx0     );
			rZZ0 = _mm256_load_pd( zz + idx_zz     );
			rZS0 = _mm256_load_pd( zs + idx_zs    );
			rZZ0 = _mm256_mul_pd(rZZ0,rC);
#ifdef AVX_FMA
			rH0 = _mm256_fmadd_pd(rZZ0, rZS0, rH0);
#else
			rZZ0 = _mm256_mul_pd(rZZ0,rZS0);
			rH0  = _mm256_add_pd(rH0,rZZ0);
#endif
			_mm256_store_pd(H + idx0, rH0);
			idx0  +=4;
			idx_zs+=4; 
			idx_zz+=4;
		    }
		    idx0 += incrj; 
		}
		idx0 += incri; 
	    }
	}
	else // H[idx0] is 8-aligned, preventing nice vectorization
	{
	    for(; i<i_end; i++)
	    {
	    	for(j = 0; j<p; j++)
	    	{
	    	    rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	    	    idx_zz=p*n;
	    	    for(k = 0; k<p; k+=4)
	    	    {
	    		rH0  = _mm256_loadu_pd( H+idx0 );
	    		rZZ0 = _mm256_load_pd( zz + idx_zz );
	    		rZS0 = _mm256_load_pd( zs + idx_zs );
	    		rZZ0 = _mm256_mul_pd(rZZ0,rC);
#ifdef AVX_FMA
			rH0 = _mm256_fmadd_pd(rZZ0, rZS0, rH0);
#else
			rZZ0 = _mm256_mul_pd(rZZ0,rZS0);
			rH0  = _mm256_add_pd(rH0,rZZ0);
#endif
	    		_mm256_storeu_pd( H+idx0, rH0 );

	    		idx0  +=4;
	    		idx_zs+=4;
	    		idx_zz+=4;
	    	    }
	    	    idx0 += incrj;
	    	}
	    	idx0 += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_AVX_P16(SE_FGG_work* work, const SE_state* st, 
			       const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;


    double qn;
    int idx, idx_zs, i, j;

    __m256d rZZ0, rZZ1, rZZ2, rZZ3; 
    __m256d rH0, rH1, rH2, rH3;
    __m256d rC, rZS0, rZS1, rZS2, rZS3;

    const int incrj = params->npdims[2]-16; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-16);// outer increment

    int i_end;
    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, 16);

    for(int n=0; n<N; n++)
    {
	idx = work->idx[n];
	idx_zs = 0;

	grid_thrd_slice(&grid_thrd_ws, &idx, &i, &i_end, &idx_zs);
	if (grid_thrd_ws.skip) 
	    continue;

	qn = st->q[n];

        rZZ0 = _mm256_load_pd(zz + n*16     );
        rZZ1 = _mm256_load_pd(zz + n*16 + 4 );
        rZZ2 = _mm256_load_pd(zz + n*16 + 8 );
        rZZ3 = _mm256_load_pd(zz + n*16 + 12);

	if(idx%4 == 0) // H[idx0] is 32-aligned
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

                    rH0  = _mm256_load_pd( H+idx    );
                    rH1  = _mm256_load_pd( H+idx + 4);
		    rH2  = _mm256_load_pd( H+idx + 8 );
                    rH3  = _mm256_load_pd( H+idx + 12);

                    rZS0 = _mm256_load_pd( zs + idx_zs);
                    rZS1 = _mm256_load_pd( zs + idx_zs + 4);
                    rZS2 = _mm256_load_pd( zs + idx_zs + 8);
                    rZS3 = _mm256_load_pd( zs + idx_zs + 12);
#ifdef AVX_FMA
		    rH0 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ0,rC), rZS0, rH0);
		    rH1 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ1,rC), rZS1, rH1);
		    rH2 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ2,rC), rZS2, rH2);
		    rH3 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ3,rC), rZS3, rH3);
#else
		    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
		    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		    rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2));
		    rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3));
#endif
		    _mm256_store_pd(H + idx,      rH0);
		    _mm256_store_pd(H + idx + 4,  rH1);
		    _mm256_store_pd(H + idx + 8 , rH2);
		    _mm256_store_pd(H + idx + 12, rH3);

		    idx += incrj + 16;
		    idx_zs += 16;
		}
		idx += incri;
	    }
	}
	else // H[idx0] is 16-aligned, preventing nice vectorization
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<16; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[16*n+i]*zy[16*n+j] );

                    rH0  = _mm256_loadu_pd( H+idx     );
                    rH1  = _mm256_loadu_pd( H+idx +  4);
                    rH2  = _mm256_loadu_pd( H+idx +  8);
                    rH3  = _mm256_loadu_pd( H+idx + 12);

                    rZS0 = _mm256_load_pd( zs + idx_zs     );
                    rZS1 = _mm256_load_pd( zs + idx_zs +  4);
                    rZS2 = _mm256_load_pd( zs + idx_zs +  8);                   
                    rZS3 = _mm256_load_pd( zs + idx_zs + 12);
#ifdef AVX_FMA
		    rH0 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ0,rC), rZS0, rH0);
		    rH1 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ1,rC), rZS1, rH1);
		    rH2 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ2,rC), rZS2, rH2);
		    rH3 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ3,rC), rZS3, rH3);
#else
		    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
		    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
		    rH2 = _mm256_add_pd(rH2,_mm256_mul_pd(_mm256_mul_pd(rZZ2,rC),rZS2));
		    rH3 = _mm256_add_pd(rH3,_mm256_mul_pd(_mm256_mul_pd(rZZ3,rC),rZS3));
#endif
		    _mm256_storeu_pd(H + idx,      rH0);
		    _mm256_storeu_pd(H + idx + 4,  rH1);
		    _mm256_storeu_pd(H + idx + 8,  rH2);
		    _mm256_storeu_pd(H + idx + 12, rH3);

		    idx += incrj + 16;
		    idx_zs += 16;
		}
		idx += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_AVX_P8(SE_FGG_work* work, const SE_state* st, 
			       const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;

    double qn;
    int idx, idx_zs, i, j;
    const int incrj = params->npdims[2]-8; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-8);// outer increment

    __m256d rZZ0, rZZ1; 
    __m256d rH0, rH1;
    __m256d rC, rZS0, rZS1;

    int i_end;
    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, 8);

    for(int n=0; n<N; n++)
    {
	idx = work->idx[n];
	idx_zs = 0;

	grid_thrd_slice(&grid_thrd_ws, &idx, &i, &i_end, &idx_zs);
	if (grid_thrd_ws.skip) 
	    continue;

	qn = st->q[n];
        rZZ0 = _mm256_load_pd(zz + n*8     );
        rZZ1 = _mm256_load_pd(zz + n*8 + 4 );

	if(idx%4 == 0) // H[idx0] is 32-aligned
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<8; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[8*n+i]*zy[8*n+j] );

                    rH0  = _mm256_load_pd( H+idx    );
                    rH1  = _mm256_load_pd( H+idx + 4);

                    rZS0 = _mm256_load_pd( zs + idx_zs);
                    rZS1 = _mm256_load_pd( zs + idx_zs + 4);
#ifdef AVX_FMA
		    rH0 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ0,rC), rZS0, rH0);
		    rH1 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ1,rC), rZS1, rH1);
#else
		    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
		    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
#endif
		    _mm256_store_pd(H + idx,      rH0);
		    _mm256_store_pd(H + idx + 4,  rH1);

		    idx += incrj + 8;
		    idx_zs += 8;
		}
		idx += incri;
	    }
	}
	else // H[idx0] is 16-aligned, preventing nice vectorization
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<8; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[8*n+i]*zy[8*n+j] );

                    rH0  = _mm256_loadu_pd( H+idx     );
                    rH1  = _mm256_loadu_pd( H+idx +  4);

                    rZS0 = _mm256_load_pd( zs + idx_zs     );
                    rZS1 = _mm256_load_pd( zs + idx_zs +   4);
#ifdef AVX_FMA
		    rH0 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ0,rC), rZS0, rH0);
		    rH1 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ1,rC), rZS1, rH1);
#else
		    rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
		    rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
#endif
		    _mm256_storeu_pd(H + idx,      rH0);
		    _mm256_storeu_pd(H + idx + 4,  rH1);

		    idx += incrj + 8;
		    idx_zs += 8;
		}
		idx += incri;
	    }
	}
    }
}

// -----------------------------------------------------------------------------
void SE_FGG_grid_split_AVX_u8(SE_FGG_work* work, const SE_state* st, 
			      const SE_FGG_params* params)
{
    // unpack parameters
    const int N=params->N;
    double* restrict H = work->H; // pointer to grid does NOT alias
    const double* restrict zs = work->zs;
    const double* restrict zx = work->zx;
    const double* restrict zy = work->zy;
    const double* restrict zz = work->zz;
    
    const int p = params->P;

    double qn;
    int idx0, idx_zs, idx_zz, i, j, k;
    const int incrj = params->npdims[2]-p; // middle increment
    const int incri = params->npdims[2]*(params->npdims[1]-p);// outer increment

    __m256d rH0, rZZ0, rZS0, rC;
    __m256d rH1, rZZ1, rZS1;

    int i_end;
    grid_thrd_ws_t grid_thrd_ws;
    grid_thrd_setup(&grid_thrd_ws, params->npdims, p);

    for(int n=0; n<N; n++)
    {
	idx0 = work->idx[n];
	idx_zs = 0;

	grid_thrd_slice(&grid_thrd_ws, &idx0, &i, &i_end, &idx_zs);
	if (grid_thrd_ws.skip) 
	    continue;

	qn = st->q[n];

	if(idx0%4 == 0) // H[idx0] is 32-aligned
	{
	    for(; i<i_end; i++)
	    {
		for(j = 0; j<p; j++)
		{
		    rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
		    idx_zz=p*n;

		    for(k = 0; k<p; k+=8)
		    {
			rH0  = _mm256_load_pd( H+idx0     );
			rH1  = _mm256_load_pd( H+idx0 + 4 );

			rZZ0 = _mm256_load_pd( zz + idx_zz     );
			rZZ1 = _mm256_load_pd( zz + idx_zz + 4 );

			rZS0 = _mm256_load_pd( zs + idx_zs    );
			rZS1 = _mm256_load_pd( zs + idx_zs + 4);
#ifdef AVX_FMA
			rH0 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ0,rC), rZS0, rH0);
			rH1 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ1,rC), rZS1, rH1);
#else
			rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
			rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
#endif
			_mm256_store_pd( H+idx0    , rH0 );
			_mm256_store_pd( H+idx0 + 4, rH1 );

			idx0  +=8;
			idx_zs+=8; 
			idx_zz+=8;
		    }
		    idx0 += incrj;
		}
		idx0 += incri;
	    }
	}
	else // H[idx0] is 16-aligned, preventing nice vectorization
	{
	    for(; i<i_end; i++)
	    {
	    	for(j = 0; j<p; j++)
	    	{
	    	    rC = _mm256_set1_pd( qn*zx[p*n+i]*zy[p*n+j] );
	    	    idx_zz=p*n;
		    
	    	    for(k = 0; k<p; k+=8)
	    	    {
	    		rH0  = _mm256_loadu_pd( H+idx0     );
	    		rH1  = _mm256_loadu_pd( H+idx0 + 4 );

	    		rZZ0 = _mm256_load_pd( zz + idx_zz     );
	    		rZZ1 = _mm256_load_pd( zz + idx_zz + 4 );

	    		rZS0 = _mm256_load_pd( zs + idx_zs    );
	    		rZS1 = _mm256_load_pd( zs + idx_zs + 4);
#ifdef AVX_FMA
			rH0 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ0,rC), rZS0, rH0);
			rH1 = _mm256_fmadd_pd(_mm256_mul_pd(rZZ1,rC), rZS1, rH1);
#else
			rH0 = _mm256_add_pd(rH0,_mm256_mul_pd(_mm256_mul_pd(rZZ0,rC),rZS0));
			rH1 = _mm256_add_pd(rH1,_mm256_mul_pd(_mm256_mul_pd(rZZ1,rC),rZS1));
#endif
	    		_mm256_storeu_pd( H+idx0    , rH0 );
	    		_mm256_storeu_pd( H+idx0 + 4, rH1 );

	    		idx0  +=8;
	    		idx_zs+=8;
	    		idx_zz+=8;
	    	    }
	    	    idx0 += incrj;
	    	}
	    	idx0 += incri;
	    }
	}
    }
}
#endif // AVX



// =============================================================================
// Reordering ==================================================================

// -----------------------------------------------------------------------------
static
int compare_idx_pair(const void* a, const void* b)
{
    int temp = ( * ((idx_reorder_t*) a) ).idx_on_grid - 
	( * ((idx_reorder_t*) b) ).idx_on_grid;
    if (temp > 0)
	return 1;
    else if (temp < 0)
	return -1;
    else
	return 0;
}

// -----------------------------------------------------------------------------
static
SE_state* SE_clone_state(const SE_state* s, const SE_FGG_params* params)
{
    const int N = params->N;

    SE_state* t = (SE_state*) SE_FGG_MALLOC( sizeof(SE_state) );
    t->x = (double*) SE_FGG_MALLOC(3*N*sizeof(double));
    t->q = (double*) SE_FGG_MALLOC(  N*sizeof(double));
    
    memcpy(t->x, s->x, 3*N*sizeof(double));
    memcpy(t->q, s->q,   N*sizeof(double));

    return t;
}

// -----------------------------------------------------------------------------
void SE_FGG_reorder_system(SE_state* s, 
			   const SE_FGG_work* work, 
			   const SE_FGG_params* params)
{
    const int N = params->N;

    /* pairing of grid-index and index in array */
    idx_reorder_t* order = (idx_reorder_t*) SE_FGG_MALLOC(N*sizeof(idx_reorder_t));
    
    /* fill index pairing */
    for(int i=0;i<N;i++)
    {
	order[i].idx_on_grid=work->idx[i];
	order[i].idx_in_array=i;
    }

    /* sort the temporary */
    qsort(order, N, sizeof(idx_reorder_t), compare_idx_pair);

    /* copy system */
    SE_state* s_copy = SE_clone_state(s, params);

    /* permute system */ 
    int j;
    for(int i=0; i<N; i++)
    {
	j = order[i].idx_in_array;

    	s->x[i]     = s_copy->x[j  ];
    	s->x[i+N]   = s_copy->x[j+N];
    	s->x[i+2*N] = s_copy->x[j+2*N];
    	s->q[i]     = s_copy->q[j];
    }

    /* deallocations */
    SE_FGG_FREE(order);
    SE_free_system(s_copy); // free contents of s_copy
    SE_FGG_FREE(s_copy);           // s_copy itself is malloc'd
}

double calc_energy(SE_state st, int N)
{
  int i;
  double energy=0;

#ifdef _OPENMP
#pragma omp parallel for private(i) reduction(+:energy)
#endif
  for(i=0;i<N;i++)
    energy += (st.phi[i])*(st.q[i]);

return energy/2.;
}

#ifdef FORCE
#include "fgg_force.c"
#endif
