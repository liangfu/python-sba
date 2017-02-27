/*
sbaprojs.c
Dennis Evangelista, 2013
copied from Lourakis' sba demo eucsbademo.c and from imageproj.c
(eucsbademo = Euclidian bundle adjustment demo using sba package)

This library provides projections and Jacobians for projecting 
from Lourakis' parameterization to points (e.g. from P=[a b] to X),
in both simple and expert mode. The return types and argtypes are
compatible with sba. 

I also included some defines here (which Python wrappings of this do
not see); and the definition of his globs_ struct which is one of
the input arguments to some of these projections. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include <sba.h> // This must be installed in /usr/local/include
#include "sbaprojs.h"


#define CLOCKS_PER_MSEC (CLOCKS_PER_SEC/1000.0)

#define MAXITER         100
#define MAXITER2        150


/* pointers to additional data, used for computed image projections and their jacobians */
struct globs_ globs;

/* unit quaternion from vector part */
#define _MK_QUAT_FRM_VEC(q, v){                                     \
  (q)[1]=(v)[0]; (q)[2]=(v)[1]; (q)[3]=(v)[2];                      \
  (q)[0]=sqrt(1.0 - (q)[1]*(q)[1] - (q)[2]*(q)[2]- (q)[3]*(q)[3]);  \
}

/*
 * multiplication of the two quaternions in q1 and q2 into p
 */
inline static void quatMult(double q1[FULLQUATSZ], double q2[FULLQUATSZ], double p[FULLQUATSZ])
{
  p[0]=q1[0]*q2[0]-q1[1]*q2[1]-q1[2]*q2[2]-q1[3]*q2[3];
  p[1]=q1[0]*q2[1]+q2[0]*q1[1]+q1[2]*q2[3]-q1[3]*q2[2];
  p[2]=q1[0]*q2[2]+q2[0]*q1[2]+q2[1]*q1[3]-q1[1]*q2[3];
  p[3]=q1[0]*q2[3]+q2[0]*q1[3]+q1[1]*q2[2]-q2[1]*q1[2];
}

/*
 * fast multiplication of the two quaternions in q1 and q2 into p
 * this is the second of the two schemes derived in pg. 8 of
 * T. D. Howell, J.-C. Lafon, The complexity of the quaternion product, TR 75-245, Cornell Univ., June 1975.
 *
 * total additions increase from 12 to 27 (28), but total multiplications decrease from 16 to 9 (12)
 */
inline static void quatMultFast(double q1[FULLQUATSZ], double q2[FULLQUATSZ], double p[FULLQUATSZ])
{
double t1, t2, t3, t4, t5, t6, t7, t8, t9;
//double t10, t11, t12;

  t1=(q1[0]+q1[1])*(q2[0]+q2[1]);
  t2=(q1[3]-q1[2])*(q2[2]-q2[3]);
  t3=(q1[1]-q1[0])*(q2[2]+q2[3]);
  t4=(q1[2]+q1[3])*(q2[1]-q2[0]);
  t5=(q1[1]+q1[3])*(q2[1]+q2[2]);
  t6=(q1[1]-q1[3])*(q2[1]-q2[2]);
  t7=(q1[0]+q1[2])*(q2[0]-q2[3]);
  t8=(q1[0]-q1[2])*(q2[0]+q2[3]);

#if 0
  t9 =t5+t6;
  t10=t7+t8;
  t11=t5-t6;
  t12=t7-t8;

  p[0]= t2 + 0.5*(-t9+t10);
  p[1]= t1 - 0.5*(t9+t10);
  p[2]=-t3 + 0.5*(t11+t12);
  p[3]=-t4 + 0.5*(t11-t12);
#endif

  /* following fragment it equivalent to the one above */
  t9=0.5*(t5-t6+t7+t8);
  p[0]= t2 + t9-t5;
  p[1]= t1 - t9-t6;
  p[2]=-t3 + t9-t8;
  p[3]=-t4 + t9-t7;
}


/* Routines to estimate the estimated measurement vector (i.e. "func") and
 * its sparse jacobian (i.e. "fjac") needed in BA. Code below makes use of the
 * routines calcImgProj() and calcImgProjJacXXX() which
 * compute the predicted projection & jacobian of a SINGLE 3D point (see imgproj.c).
 * In the terminology of TR-340, these routines compute Q and its jacobians A=dQ/da, B=dQ/db.
 * Notice also that what follows is two pairs of "func" and corresponding "fjac" routines.
 * The first is to be used in full (i.e. motion + structure) BA, the second in 
 * motion only BA.
 */

static const double zerorotquat[FULLQUATSZ]={1.0, 0.0, 0.0, 0.0};

/****************************************************************************************/
/* MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR VARYING CAMERA POSE AND 3D STRUCTURE */
/****************************************************************************************/

/*** MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE SIMPLE DRIVERS ***/

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given the parameter vectors aj and bi of camera j and point i, computes in xij the
 * predicted projection of point i on image j
 */
void img_projRTS(int j, int i, double *aj, double *bi, double *xij, void *adata)
{
  double *Kparms, *pr0;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  Kparms=gl->intrcalib;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

  calcImgProj(Kparms, pr0, aj, aj+3, bi, xij); // 3 is the quaternion's vector part length
}

/* Given the parameter vectors aj and bi of camera j and point i, computes in Aij, Bij the
 * jacobian of the predicted projection of point i on image j
 */
void img_projRTS_jac(int j, int i, double *aj, double *bi, double *Aij, double *Bij, void *adata)
{
  double *Kparms, *pr0;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  Kparms=gl->intrcalib;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

  calcImgProjJacRTS(Kparms, pr0, aj, aj+3, bi, (double (*)[6])Aij, (double (*)[3])Bij); // 3 is the quaternion's vector part length
}

/* BUNDLE ADJUSTMENT FOR CAMERA PARAMETERS ONLY */

/* Given the parameter vector aj of camera j, computes in xij the
 * predicted projection of point i on image j
 */
void img_projRT(int j, int i, double *aj, double *xij, void *adata)
{
  int pnp;

  double *Kparms, *pr0, *ptparams;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  pnp=gl->pnp;
  Kparms=gl->intrcalib;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
  ptparams=gl->ptparams;

  calcImgProj(Kparms, pr0, aj, aj+3, ptparams+i*pnp, xij); // 3 is the quaternion's vector part length
}

/* Given the parameter vector aj of camera j, computes in Aij
 * the jacobian of the predicted projection of point i on image j
 */
void img_projRT_jac(int j, int i, double *aj, double *Aij, void *adata)
{
  int pnp;

  double *Kparms, *ptparams, *pr0;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  pnp=gl->pnp;
  Kparms=gl->intrcalib;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
  ptparams=gl->ptparams;

  calcImgProjJacRT(Kparms, pr0, aj, aj+3, ptparams+i*pnp, (double (*)[6])Aij); // 3 is the quaternion's vector part length
}

/* BUNDLE ADJUSTMENT FOR STRUCTURE PARAMETERS ONLY */

/* Given the parameter vector bi of point i, computes in xij the
 * predicted projection of point i on image j
 */
void img_projS(int j, int i, double *bi, double *xij, void *adata)
{
  int cnp;

  double *Kparms, *camparams, *aj;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp;
  Kparms=gl->intrcalib;
  camparams=gl->camparams;
  aj=camparams+j*cnp;

  calcImgProjFullR(Kparms, aj, aj+3, bi, xij); // 3 is the quaternion's vector part length
  //calcImgProj(Kparms, (double *)zerorotquat, aj, aj+3, bi, xij); // 3 is the quaternion's vector part length
}

/* Given the parameter vector bi of point i, computes in Bij
 * the jacobian of the predicted projection of point i on image j
 */
void img_projS_jac(int j, int i, double *bi, double *Bij, void *adata)
{
  int cnp;

  double *Kparms, *camparams, *aj;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp;
  Kparms=gl->intrcalib;
  camparams=gl->camparams;
  aj=camparams+j*cnp;

  calcImgProjJacS(Kparms, (double *)zerorotquat, aj, aj+3, bi, (double (*)[3])Bij); // 3 is the quaternion's vector part length
}

/*** MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE EXPERT DRIVERS ***/

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsRTS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pa, *pb, *pqr, *pt, *ppt, *pmeas, *Kparms, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  Kparms=gl->intrcalib;

  //n=idxij->nr;
  m=idxij->nc;
  pa=p; pb=p+m*cnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=pa+j*cnp;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
    _MK_QUAT_FRM_VEC(lrot, pqr);
    quatMultFast(lrot, pr0, trot); // trot=lrot*pr0

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=pb + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcImgProjFullR(Kparms, trot, pt, ppt, pmeas); // evaluate Q in pmeas
      //calcImgProj(Kparms, pr0, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm, B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where A_ij=dx_ij/db_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 */
void img_projsRTS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pa, *pb, *pqr, *pt, *ppt, *pA, *pB, *Kparms, *pr0;
  //int n;
  int m, nnz, Asz, Bsz, ABsz;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  Kparms=gl->intrcalib;

  //n=idxij->nr;
  m=idxij->nc;
  pa=p; pb=p+m*cnp;
  Asz=mnp*cnp; Bsz=mnp*pnp; ABsz=Asz+Bsz;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=pa+j*cnp;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=pb + rcsubs[i]*pnp;
      pA=jac + idxij->val[rcidxs[i]]*ABsz; // set pA to point to A_ij
      pB=pA  + Asz; // set pB to point to B_ij

      calcImgProjJacRTS(Kparms, pr0, pqr, pt, ppt, (double (*)[6])pA, (double (*)[3])pB); // evaluate dQ/da, dQ/db in pA, pB
    }
  }
}

/* BUNDLE ADJUSTMENT FOR CAMERA PARAMETERS ONLY */

/* Given a parameter vector p made up of the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images.
 * The measurements are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T,
 * where hx_ij is the predicted projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsRT_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pmeas, *Kparms, *ptparams, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  Kparms=gl->intrcalib;
  ptparams=gl->ptparams;

  //n=idxij->nr;
  m=idxij->nc;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=p+j*cnp;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
    _MK_QUAT_FRM_VEC(lrot, pqr);
    quatMultFast(lrot, pr0, trot); // trot=lrot*pr0

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
	    ppt=ptparams + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcImgProjFullR(Kparms, trot, pt, ppt, pmeas); // evaluate Q in pmeas
      //calcImgProj(Kparms, pr0, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the parameters of m cameras, compute in jac
 * the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm),
 * where A_ij=dx_ij/db_j (see HZ).
 * Notice that depending on idxij, some of the A_ij might be missing
 *
 */
void img_projsRT_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pA, *Kparms, *ptparams, *pr0;
  //int n;
  int m, nnz, Asz;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  Kparms=gl->intrcalib;
  ptparams=gl->ptparams;

  //n=idxij->nr;
  m=idxij->nc;
  Asz=mnp*cnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=p+j*cnp;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=ptparams + rcsubs[i]*pnp;
      pA=jac + idxij->val[rcidxs[i]]*Asz; // set pA to point to A_ij

      calcImgProjJacRT(Kparms, pr0, pqr, pt, ppt, (double (*)[6])pA); // evaluate dQ/da in pA
    }
  }
}

/* BUNDLE ADJUSTMENT FOR STRUCTURE PARAMETERS ONLY */

/* Given a parameter vector p made up of the 3D coordinates of n points, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pmeas, *Kparms, *camparams, trot[FULLQUATSZ];
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  Kparms=gl->intrcalib;
  camparams=gl->camparams;

  //n=idxij->nr;
  m=idxij->nc;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=camparams+j*cnp;
    pt=pqr+3; // quaternion vector part has 3 elements
    _MK_QUAT_FRM_VEC(trot, pqr);

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=p + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcImgProjFullR(Kparms, trot, pt, ppt, pmeas); // evaluate Q in pmeas
      //calcImgProj(Kparms, (double *)zerorotquat, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the 3D coordinates of n points, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the B_ij might be missing
 *
 */
void img_projsS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pB, *Kparms, *camparams;
  //int n;
  int m, nnz, Bsz;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  Kparms=gl->intrcalib;
  camparams=gl->camparams;

  //n=idxij->nr;
  m=idxij->nc;
  Bsz=mnp*pnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pqr=camparams+j*cnp;
    pt=pqr+3; // quaternion vector part has 3 elements

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=p + rcsubs[i]*pnp;
      pB=jac + idxij->val[rcidxs[i]]*Bsz; // set pB to point to B_ij

      calcImgProjJacS(Kparms, (double *)zerorotquat, pqr, pt, ppt, (double (*)[3])pB); // evaluate dQ/da in pB
    }
  }
}

/****************************************************************************************************/
/* MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR VARYING CAMERA INTRINSICS, POSE AND 3D STRUCTURE */
/****************************************************************************************************/

/*** MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE SIMPLE DRIVERS ***/

/* A note about the computation of Jacobians below:
 *
 * When performing BA that includes the camera intrinsics, it would be
 * very desirable to allow for certain parameters such as skew, aspect
 * ratio and principal point to be fixed. The straighforward way to
 * implement this would be to code a separate version of the Jacobian
 * computation routines for each subset of non-fixed parameters. Here,
 * this is bypassed by developing only one set of Jacobian computation
 * routines which estimate the former for all 5 intrinsics and then set
 * the columns corresponding to fixed parameters to zero.
 */

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given the parameter vectors aj and bi of camera j and point i, computes in xij the
 * predicted projection of point i on image j
 */
void img_projKRTS(int j, int i, double *aj, double *bi, double *xij, void *adata)
{
  double *pr0;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

  calcImgProj(aj, pr0, aj+5, aj+5+3, bi, xij); // 5 for the calibration + 3 for the quaternion's vector part
}

/* Given the parameter vectors aj and bi of camera j and point i, computes in Aij, Bij the
 * jacobian of the predicted projection of point i on image j
 */
void img_projKRTS_jac(int j, int i, double *aj, double *bi, double *Aij, double *Bij, void *adata)
{
struct globs_ *gl;
double *pr0;
int ncK;

  gl=(struct globs_ *)adata;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
  calcImgProjJacKRTS(aj, pr0, aj+5, aj+5+3, bi, (double (*)[5+6])Aij, (double (*)[3])Bij); // 5 for the calibration + 3 for the quaternion's vector part

  /* clear the columns of the Jacobian corresponding to fixed calibration parameters */
  gl=(struct globs_ *)adata;
  ncK=gl->nccalib;
  if(ncK){
    int cnp, mnp, j0;

    cnp=gl->cnp;
    mnp=gl->mnp;
    j0=5-ncK;

    for(i=0; i<mnp; ++i, Aij+=cnp)
      for(j=j0; j<5; ++j)
        Aij[j]=0.0; // Aij[i*cnp+j]=0.0;
  }
}

/* BUNDLE ADJUSTMENT FOR CAMERA PARAMETERS ONLY */

/* Given the parameter vector aj of camera j, computes in xij the
 * predicted projection of point i on image j
 */
void img_projKRT(int j, int i, double *aj, double *xij, void *adata)
{
  int pnp;

  double *ptparams, *pr0;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  pnp=gl->pnp;
  ptparams=gl->ptparams;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

  calcImgProj(aj, pr0, aj+5, aj+5+3, ptparams+i*pnp, xij); // 5 for the calibration + 3 for the quaternion's vector part
}

/* Given the parameter vector aj of camera j, computes in Aij
 * the jacobian of the predicted projection of point i on image j
 */
void img_projKRT_jac(int j, int i, double *aj, double *Aij, void *adata)
{
struct globs_ *gl;
double *ptparams, *pr0;
int pnp, ncK;
  
  gl=(struct globs_ *)adata;
  pnp=gl->pnp;
  ptparams=gl->ptparams;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

  calcImgProjJacKRT(aj, pr0, aj+5, aj+5+3, ptparams+i*pnp, (double (*)[5+6])Aij); // 5 for the calibration + 3 for the quaternion's vector part

  /* clear the columns of the Jacobian corresponding to fixed calibration parameters */
  ncK=gl->nccalib;
  if(ncK){
    int cnp, mnp, j0;

    cnp=gl->cnp;
    mnp=gl->mnp;
    j0=5-ncK;

    for(i=0; i<mnp; ++i, Aij+=cnp)
      for(j=j0; j<5; ++j)
        Aij[j]=0.0; // Aij[i*cnp+j]=0.0;
  }
}

/* BUNDLE ADJUSTMENT FOR STRUCTURE PARAMETERS ONLY */

/* Given the parameter vector bi of point i, computes in xij the
 * predicted projection of point i on image j
 */
void img_projKS(int j, int i, double *bi, double *xij, void *adata)
{
  int cnp;

  double *camparams, *aj;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp;
  camparams=gl->camparams;
  aj=camparams+j*cnp;

  calcImgProjFullR(aj, aj+5, aj+5+3, bi, xij); // 5 for the calibration + 3 for the quaternion's vector part
  //calcImgProj(aj, (double *)zerorotquat, aj+5, aj+5+3, bi, xij); // 5 for the calibration + 3 for the quaternion's vector part
}

/* Given the parameter vector bi of point i, computes in Bij
 * the jacobian of the predicted projection of point i on image j
 */
void img_projKS_jac(int j, int i, double *bi, double *Bij, void *adata)
{
  int cnp;

  double *camparams, *aj;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp;
  camparams=gl->camparams;
  aj=camparams+j*cnp;

  calcImgProjJacS(aj, (double *)zerorotquat, aj+5, aj+5+3, bi, (double (*)[3])Bij); // 5 for the calibration + 3 for the quaternion's vector part
}

/*** MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE EXPERT DRIVERS ***/

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsKRTS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pa, *pb, *pqr, *pt, *ppt, *pmeas, *pcalib, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
  //int n;
  int m, nnz;
  struct globs_ *gl;

  //fprintf(stderr,"img_projsKRTS_x\n");

  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  //fprintf(stderr,"img_projsKRTS_x unpacked glob\n");

  //n=idxij->nr;
  m=idxij->nc;
  pa=p; pb=p+m*cnp;
  //fprintf(stderr,"img_projsKRTS_x diddled idxij\n");

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=pa+j*cnp;
    //fprintf(stderr,"img_projsKRTS_x pcalib\n");
    pqr=pcalib+5;
    //fprintf(stderr,"img_projsKRTS_x pqr\n");
    pt=pqr+3; // quaternion vector part has 3 elements
    //fprintf(stderr,"img_projsKRTS_x pt\n");
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
    //fprintf(stderr,"pr0 %f %f %f %f\n",pr0[0],pr0[1],pr0[2],pr0[3]);
    _MK_QUAT_FRM_VEC(lrot, pqr);
    //fprintf(stderr,"made quats from vec imag part\n");
    quatMultFast(lrot, pr0, trot); // trot=lrot*pr0
    //fprintf(stderr,"multipleid some quats\n");

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=pb + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcImgProjFullR(pcalib, trot, pt, ppt, pmeas); // evaluate Q in pmeas
      //calcImgProj(pcalib, pr0, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm, B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where A_ij=dx_ij/db_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 */
void img_projsKRTS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j, ii, jj;
  int cnp, pnp, mnp, ncK;
  double *pa, *pb, *pqr, *pt, *ppt, *pA, *pB, *pcalib, *pr0;
  //int n;
  int m, nnz, Asz, Bsz, ABsz;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  ncK=gl->nccalib;

  //n=idxij->nr;
  m=idxij->nc;
  pa=p; pb=p+m*cnp;
  Asz=mnp*cnp; Bsz=mnp*pnp; ABsz=Asz+Bsz;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=pa+j*cnp;
    pqr=pcalib+5;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=pb + rcsubs[i]*pnp;
      pA=jac + idxij->val[rcidxs[i]]*ABsz; // set pA to point to A_ij
      pB=pA  + Asz; // set pB to point to B_ij

      calcImgProjJacKRTS(pcalib, pr0, pqr, pt, ppt, (double (*)[5+6])pA, (double (*)[3])pB); // evaluate dQ/da, dQ/db in pA, pB

      /* clear the columns of the Jacobian corresponding to fixed calibration parameters */
      if(ncK){
        int jj0=5-ncK;

        for(ii=0; ii<mnp; ++ii, pA+=cnp)
          for(jj=jj0; jj<5; ++jj)
            pA[jj]=0.0; // pA[ii*cnp+jj]=0.0;
      }
    }
  }
}

/* BUNDLE ADJUSTMENT FOR CAMERA PARAMETERS ONLY */

/* Given a parameter vector p made up of the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images.
 * The measurements are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T,
 * where hx_ij is the predicted projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsKRT_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pmeas, *pcalib, *ptparams, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  ptparams=gl->ptparams;

  //n=idxij->nr;
  m=idxij->nc;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=p+j*cnp;
    pqr=pcalib+5;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
    _MK_QUAT_FRM_VEC(lrot, pqr);
    quatMultFast(lrot, pr0, trot); // trot=lrot*pr0

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
	    ppt=ptparams + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcImgProjFullR(pcalib, trot, pt, ppt, pmeas); // evaluate Q in pmeas
      //calcImgProj(pcalib, pr0, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the parameters of m cameras, compute in jac
 * the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm),
 * where A_ij=dx_ij/db_j (see HZ).
 * Notice that depending on idxij, some of the A_ij might be missing
 *
 */
void img_projsKRT_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j, ii, jj;
  int cnp, pnp, mnp, ncK;
  double *pqr, *pt, *ppt, *pA, *pcalib, *ptparams, *pr0;
  //int n;
  int m, nnz, Asz;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  ncK=gl->nccalib;
  ptparams=gl->ptparams;

  //n=idxij->nr;
  m=idxij->nc;
  Asz=mnp*cnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=p+j*cnp;
    pqr=pcalib+5;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=ptparams + rcsubs[i]*pnp;
      pA=jac + idxij->val[rcidxs[i]]*Asz; // set pA to point to A_ij

      calcImgProjJacKRT(pcalib, pr0, pqr, pt, ppt, (double (*)[5+6])pA); // evaluate dQ/da in pA

      /* clear the columns of the Jacobian corresponding to fixed calibration parameters */
      if(ncK){
        int jj0;

        jj0=5-ncK;
        for(ii=0; ii<mnp; ++ii, pA+=cnp)
          for(jj=jj0; jj<5; ++jj)
            pA[jj]=0.0; // pA[ii*cnp+jj]=0.0;
      }
    }
  }
}

/* BUNDLE ADJUSTMENT FOR STRUCTURE PARAMETERS ONLY */

/* Given a parameter vector p made up of the 3D coordinates of n points, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsKS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pmeas, *pcalib, *camparams, trot[FULLQUATSZ];
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  camparams=gl->camparams;

  //n=idxij->nr;
  m=idxij->nc;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=camparams+j*cnp;
    pqr=pcalib+5;
    pt=pqr+3; // quaternion vector part has 3 elements
    _MK_QUAT_FRM_VEC(trot, pqr);

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=p + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcImgProjFullR(pcalib, trot, pt, ppt, pmeas); // evaluate Q in pmeas
      //calcImgProj(pcalib, (double *)zerorotquat, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the 3D coordinates of n points, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the B_ij might be missing
 *
 */
void img_projsKS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pB, *pcalib, *camparams;
  //int n;
  int m, nnz, Bsz;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  camparams=gl->camparams;

  //n=idxij->nr;
  m=idxij->nc;
  Bsz=mnp*pnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=camparams+j*cnp;
    pqr=pcalib+5;
    pt=pqr+3; // quaternion vector part has 3 elements

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=p + rcsubs[i]*pnp;
      pB=jac + idxij->val[rcidxs[i]]*Bsz; // set pB to point to B_ij

      calcImgProjJacS(pcalib, (double *)zerorotquat, pqr, pt, ppt, (double (*)[3])pB); // evaluate dQ/da in pB
    }
  }
}

/****************************************************************************************************************/
/* MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR VARYING CAMERA INTRINSICS, DISTORTION, POSE AND 3D STRUCTURE */
/****************************************************************************************************************/

/*** MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE SIMPLE DRIVERS ***/

/* A note about the computation of Jacobians below:
 *
 * When performing BA that includes the camera intrinsics & distortion, it would be
 * very desirable to allow for certain parameters such as skew, aspect ratio and principal
 * point (also high order radial distortion, tangential distortion), to be fixed. The
 * straighforward way to implement this would be to code a separate version of the
 * Jacobian computation routines for each subset of non-fixed parameters. Here,
 * this is bypassed by developing only one set of Jacobian computation
 * routines which estimate the former for all 5 intrinsics and all 5 distortion
 * coefficients and then set the columns corresponding to fixed parameters to zero.
 */

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given the parameter vectors aj and bi of camera j and point i, computes in xij the
 * predicted projection of point i on image j
 */
void img_projKDRTS(int j, int i, double *aj, double *bi, double *xij, void *adata)
{
  double *pr0;
  struct globs_ *gl;
  //fprintf(stderr,"Got to img_projKDRTS\n");

  gl=(struct globs_ *)adata;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

  calcDistImgProj(aj, aj+5, pr0, aj+5+5, aj+5+5+3, bi, xij); // 5 for the calibration + 5 for the distortion + 3 for the quaternion's vector part
}

/* Given the parameter vectors aj and bi of camera j and point i, computes in Aij, Bij the
 * jacobian of the predicted projection of point i on image j
 */
void img_projKDRTS_jac(int j, int i, double *aj, double *bi, double *Aij, double *Bij, void *adata)
{
struct globs_ *gl;
double *pA, *pr0;
int nc;

  gl=(struct globs_ *)adata;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
  calcDistImgProjJacKDRTS(aj, aj+5, pr0, aj+5+5, aj+5+5+3, bi, (double (*)[5+5+6])Aij, (double (*)[3])Bij); // 5 for the calibration + 5 for the distortion + 3 for the quaternion's vector part

  /* clear the columns of the Jacobian corresponding to fixed calibration parameters */
  gl=(struct globs_ *)adata;
  nc=gl->nccalib;
  if(nc){
    int cnp, mnp, j0;

    pA=Aij;
    cnp=gl->cnp;
    mnp=gl->mnp;
    j0=5-nc;

    for(i=0; i<mnp; ++i, pA+=cnp)
      for(j=j0; j<5; ++j)
        pA[j]=0.0; // pA[i*cnp+j]=0.0;
  }

  /* clear the columns of the Jacobian corresponding to fixed distortion parameters */
  nc=gl->ncdist;
  if(nc){
    int cnp, mnp, j0;

    pA=Aij;
    cnp=gl->cnp;
    mnp=gl->mnp;
    j0=5-nc;

    for(i=0; i<mnp; ++i, pA+=cnp)
      for(j=j0; j<5; ++j)
        pA[5+j]=0.0; // pA[i*cnp+5+j]=0.0;
  }
}

/* BUNDLE ADJUSTMENT FOR CAMERA PARAMETERS ONLY */

/* Given the parameter vector aj of camera j, computes in xij the
 * predicted projection of point i on image j
 */
void img_projKDRT(int j, int i, double *aj, double *xij, void *adata)
{
  int pnp;

  double *ptparams, *pr0;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  pnp=gl->pnp;
  ptparams=gl->ptparams;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

  calcDistImgProj(aj, aj+5, pr0, aj+5+5, aj+5+5+3, ptparams+i*pnp, xij); // 5 for the calibration + 5 for the distortion + 3 for the quaternion's vector part
}

/* Given the parameter vector aj of camera j, computes in Aij
 * the jacobian of the predicted projection of point i on image j
 */
void img_projKDRT_jac(int j, int i, double *aj, double *Aij, void *adata)
{
struct globs_ *gl;
double *pA, *ptparams, *pr0;
int pnp, nc;
  
  gl=(struct globs_ *)adata;
  pnp=gl->pnp;
  ptparams=gl->ptparams;
  pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

  calcDistImgProjJacKDRT(aj, aj+5, pr0, aj+5+5, aj+5+5+3, ptparams+i*pnp, (double (*)[5+5+6])Aij); // 5 for the calibration + 5 for the distortion + 3 for the quaternion's vector part

  /* clear the columns of the Jacobian corresponding to fixed calibration parameters */
  nc=gl->nccalib;
  if(nc){
    int cnp, mnp, j0;

    pA=Aij;
    cnp=gl->cnp;
    mnp=gl->mnp;
    j0=5-nc;

    for(i=0; i<mnp; ++i, pA+=cnp)
      for(j=j0; j<5; ++j)
        pA[j]=0.0; // pA[i*cnp+j]=0.0;
  }
  nc=gl->ncdist;
  if(nc){
    int cnp, mnp, j0;

    pA=Aij;
    cnp=gl->cnp;
    mnp=gl->mnp;
    j0=5-nc;

    for(i=0; i<mnp; ++i, pA+=cnp)
      for(j=j0; j<5; ++j)
        pA[5+j]=0.0; // pA[i*cnp+5+j]=0.0;
  }
}

/* BUNDLE ADJUSTMENT FOR STRUCTURE PARAMETERS ONLY */

/* Given the parameter vector bi of point i, computes in xij the
 * predicted projection of point i on image j
 */
void img_projKDS(int j, int i, double *bi, double *xij, void *adata)
{
  int cnp;

  double *camparams, *aj;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp;
  camparams=gl->camparams;
  aj=camparams+j*cnp;

  calcDistImgProjFullR(aj, aj+5, aj+5+5, aj+5+5+3, bi, xij); // 5 for the calibration + 5 for the distortion + 3 for the quaternion's vector part
  //calcDistImgProj(aj, aj+5, (double *)zerorotquat, aj+5+5, aj+5+5+3, bi, xij); // 5 for the calibration + 5 for the distortion + 3 for the quaternion's vector part
}

/* Given the parameter vector bi of point i, computes in Bij the
 * jacobian of the predicted projection of point i on image j
 */
void img_projKDS_jac(int j, int i, double *bi, double *Bij, void *adata)
{
  int cnp;

  double *camparams, *aj;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp;
  camparams=gl->camparams;
  aj=camparams+j*cnp;

  calcDistImgProjJacS(aj, aj+5, (double *)zerorotquat, aj+5+5, aj+5+5+3, bi, (double (*)[3])Bij); // 5 for the calibration + 5 for the distortion + 3 for the quaternion's vector part
}


/*** MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE EXPERT DRIVERS ***/

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsKDRTS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pa, *pb, *pqr, *pt, *ppt, *pmeas, *pcalib, *pdist, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
  //int n;
  int m, nnz;
  struct globs_ *gl;

  //fprintf(stderr,"Got to img_projsKDRTS_x\n");

  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  //fprintf(stderr,"unpacked globs OK\n");

  //n=idxij->nr;
  m=idxij->nc;
  pa=p; pb=p+m*cnp;

  //fprintf(stderr,"just before for loop\n");
  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    //fprintf(stderr,"in the loop\n");
    pcalib=pa+j*cnp;
    pdist=pcalib+5;
    //fprintf(stderr,"before dist\n");
    pqr=pdist+5;
    //fprintf(stderr,"after dist\n");
    pt=pqr+3; // quaternion vector part has 3 elements
    //fprintf(stderr,"after quat thing\n");
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
    //fprintf(stderr,"got rot0params from glob\n");
    _MK_QUAT_FRM_VEC(lrot, pqr);
    //fprintf(stderr,"made quats from vec imag part\n");
    quatMultFast(lrot, pr0, trot); // trot=lrot*pr0
    //fprintf(stderr,"multiplied some quaternions\n");
  
    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */
    //fprintf(stderr,"messed with crsm\n");

    for(i=0; i<nnz; ++i){
      ppt=pb + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      //fprintf(stderr,"about to calculate dist img proj ful r\n");
      calcDistImgProjFullR(pcalib, pdist, trot, pt, ppt, pmeas); // evaluate Q in pmeas
      //calcDistImgProj(pcalib, pdist, pr0, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm, B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where A_ij=dx_ij/db_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 */
void img_projsKDRTS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j, ii, jj;
  int cnp, pnp, mnp, ncK, ncD;
  double *pa, *pb, *pqr, *pt, *ppt, *pA, *pB, *ptr, *pcalib, *pdist, *pr0;
  //int n;
  int m, nnz, Asz, Bsz, ABsz;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  ncK=gl->nccalib;
  ncD=gl->ncdist;

  //n=idxij->nr;
  m=idxij->nc;
  pa=p; pb=p+m*cnp;
  Asz=mnp*cnp; Bsz=mnp*pnp; ABsz=Asz+Bsz;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=pa+j*cnp;
    pdist=pcalib+5;
    pqr=pdist+5;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=pb + rcsubs[i]*pnp;
      pA=jac + idxij->val[rcidxs[i]]*ABsz; // set pA to point to A_ij
      pB=pA  + Asz; // set pB to point to B_ij

      calcDistImgProjJacKDRTS(pcalib, pdist, pr0, pqr, pt, ppt, (double (*)[5+5+6])pA, (double (*)[3])pB); // evaluate dQ/da, dQ/db in pA, pB

      /* clear the columns of the Jacobian corresponding to fixed calibration parameters */
      if(ncK){
        int jj0=5-ncK;

        ptr=pA;
        for(ii=0; ii<mnp; ++ii, ptr+=cnp)
          for(jj=jj0; jj<5; ++jj)
            ptr[jj]=0.0; // ptr[ii*cnp+jj]=0.0;
      }

      /* clear the columns of the Jacobian corresponding to fixed distortion parameters */
      if(ncD){
        int jj0=5-ncD;

        ptr=pA;
        for(ii=0; ii<mnp; ++ii, ptr+=cnp)
          for(jj=jj0; jj<5; ++jj)
            ptr[5+jj]=0.0; // ptr[ii*cnp+5+jj]=0.0;
      }
    }
  }
}

/* BUNDLE ADJUSTMENT FOR CAMERA PARAMETERS ONLY */

/* Given a parameter vector p made up of the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images.
 * The measurements are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T,
 * where hx_ij is the predicted projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsKDRT_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pmeas, *pcalib, *pdist, *ptparams, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  ptparams=gl->ptparams;

  //n=idxij->nr;
  m=idxij->nc;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=p+j*cnp;
    pdist=pcalib+5;
    pqr=pdist+5;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
    _MK_QUAT_FRM_VEC(lrot, pqr);
    quatMultFast(lrot, pr0, trot); // trot=lrot*pr0

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
	    ppt=ptparams + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcDistImgProjFullR(pcalib, pdist, trot, pt, ppt, pmeas); // evaluate Q in pmeas
      //calcDistImgProj(pcalib, pdist, pr0, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the parameters of m cameras, compute in jac
 * the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm),
 * where A_ij=dx_ij/db_j (see HZ).
 * Notice that depending on idxij, some of the A_ij might be missing
 *
 */
void img_projsKDRT_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j, ii, jj;
  int cnp, pnp, mnp, ncK, ncD;
  double *pqr, *pt, *ppt, *pA, *ptr, *pcalib, *pdist, *ptparams, *pr0;
  //int n;
  int m, nnz, Asz;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  ncK=gl->nccalib;
  ncD=gl->ncdist;
  ptparams=gl->ptparams;

  //n=idxij->nr;
  m=idxij->nc;
  Asz=mnp*cnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=p+j*cnp;
    pdist=pcalib+5;
    pqr=pdist+5;
    pt=pqr+3; // quaternion vector part has 3 elements
    pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=ptparams + rcsubs[i]*pnp;
      pA=jac + idxij->val[rcidxs[i]]*Asz; // set pA to point to A_ij

      calcDistImgProjJacKDRT(pcalib, pdist, pr0, pqr, pt, ppt, (double (*)[5+5+6])pA); // evaluate dQ/da in pA

      /* clear the columns of the Jacobian corresponding to fixed calibration parameters */
      if(ncK){
        int jj0;

        ptr=pA;
        jj0=5-ncK;
        for(ii=0; ii<mnp; ++ii, ptr+=cnp)
          for(jj=jj0; jj<5; ++jj)
            ptr[jj]=0.0; // ptr[ii*cnp+jj]=0.0;
      }

      /* clear the columns of the Jacobian corresponding to fixed distortion parameters */
      if(ncD){
        int jj0;

        ptr=pA;
        jj0=5-ncD;
        for(ii=0; ii<mnp; ++ii, ptr+=cnp)
          for(jj=jj0; jj<5; ++jj)
            ptr[5+jj]=0.0; // ptr[ii*cnp+5+jj]=0.0;
      }
    }
  }
}

/* BUNDLE ADJUSTMENT FOR STRUCTURE PARAMETERS ONLY */

/* Given a parameter vector p made up of the 3D coordinates of n points, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsKDS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pmeas, *pcalib, *pdist, *camparams, trot[FULLQUATSZ];
  //int n;
  int m, nnz;
  struct globs_ *gl;

  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  camparams=gl->camparams;

  //n=idxij->nr;
  m=idxij->nc;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=camparams+j*cnp;
    pdist=pcalib+5;
    pqr=pdist+5;
    pt=pqr+3; // quaternion vector part has 3 elements
    _MK_QUAT_FRM_VEC(trot, pqr);

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=p + rcsubs[i]*pnp;
      pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

      calcDistImgProjFullR(pcalib, pdist, trot, pt, ppt, pmeas); // evaluate Q in pmeas
      //calcDistImgProj(pcalib, pdist, (double *)zerorotquat, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
    }
  }
}

/* Given a parameter vector p made up of the 3D coordinates of n points, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the B_ij might be missing
 *
 */
void img_projsKDS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
  register int i, j;
  int cnp, pnp, mnp;
  double *pqr, *pt, *ppt, *pB, *pcalib, *pdist, *camparams;
  //int n;
  int m, nnz, Bsz;
  struct globs_ *gl;
  
  gl=(struct globs_ *)adata;
  cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
  camparams=gl->camparams;

  //n=idxij->nr;
  m=idxij->nc;
  Bsz=mnp*pnp;

  for(j=0; j<m; ++j){
    /* j-th camera parameters */
    pcalib=camparams+j*cnp;
    pdist=pcalib+5;
    pqr=pdist+5;
    pt=pqr+3; // quaternion vector part has 3 elements

    nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

    for(i=0; i<nnz; ++i){
      ppt=p + rcsubs[i]*pnp;
      pB=jac + idxij->val[rcidxs[i]]*Bsz; // set pB to point to B_ij

      calcDistImgProjJacS(pcalib, pdist, (double *)zerorotquat, pqr, pt, ppt, (double (*)[3])pB); // dQ/db in pB
    }
  }
}


