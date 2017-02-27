/*
sbaprojs.h
Dennis Evangelista, 2013
copied from Lourakis' sba demo eucsbademo.c and from imageproj.c

This library provides projections and Jacobians for projecting 
from Lourakis' parameterization to points (e.g. from P=[a b] to X),
in both simple and expert mode. The return types and argtypes are
compatible with sba. 

I also included some defines here (which Python wrappings of this do
not see); and the definition of his globs_ struct which is one of
the input arguments to some of these projections. 
*/

#ifndef _SBA_PROJS_H_
#define _SBA_PROJS_H_

#define SBA_MAX_REPROJ_ERROR    4.0 // max motion only reprojection error

#define BA_NONE                 -1
#define BA_MOTSTRUCT            0
#define BA_MOT                  1
#define BA_STRUCT               2
#define BA_MOT_MOTSTRUCT        3

#define FULLQUATSZ 4

#include <sba.h>
// is this OK to do?  need it so use of struct Crsm here in sbaprojs.h doesn't confuse things

/* in imgproj.c */
extern void calcImgProj(double a[5], double qr0[4], double v[3], double t[3], double M[3], double n[2]);
extern void calcImgProjFullR(double a[5], double qr0[4], double t[3], double M[3], double n[2]);
extern void calcImgProjJacKRTS(double a[5], double qr0[4], double v[3], double t[3], double M[3], double jacmKRT[2][11], double jacmS[2][3]);
extern void calcImgProjJacKRT(double a[5], double qr0[4], double v[3], double t[3], double M[3], double jacmKRT[2][11]);
extern void calcImgProjJacS(double a[5], double qr0[4], double v[3], double t[3], double M[3], double jacmS[2][3]);
extern void calcImgProjJacRTS(double a[5], double qr0[4], double v[3], double t[3], double M[3], double jacmRT[2][6], double jacmS[2][3]);
extern void calcImgProjJacRT(double a[5], double qr0[4], double v[3], double t[3], double M[3], double jacmRT[2][6]);
extern void calcDistImgProj(double a[5], double kc[5], double qr0[4], double v[3], double t[3], double M[3], double n[2]);
extern void calcDistImgProjFullR(double a[5], double kc[5], double qr0[4], double t[3], double M[3], double n[2]);
extern void calcDistImgProjJacKDRTS(double a[5], double kc[5], double qr0[4], double v[3], double t[3], double M[3], double jacmKDRT[2][16], double jacmS[2][3]);
extern void calcDistImgProjJacKDRT(double a[5], double kc[5], double qr0[4], double v[3], double t[3], double M[3], double jacmKDRT[2][16]);
extern void calcDistImgProjJacS(double a[5], double kc[5], double qr0[4], double v[3], double t[3], double M[3], double jacmS[2][3]);




/* Projection functions and their Jacobians in form with correct 
arguments for use with sba */

// simple mode, full bundle adjustment
extern void img_projRTS(int j, int i, double *aj, double *bi, double *xij, void *adata);
extern void img_projRTS_jac(int j, int i, double *aj, double *bi, double *Aij, double *xij, void *adata);

// simple mode, bundle adjustment for camera motion only
extern void img_projRT(int j, int i, double *aj, double *xij, void *adata);
extern void img_projRT_jac(int j, int i, double *aj, double *Aij, void *adata);

// simple mode, bundle adjustment for structure parameters only
extern void img_projS(int j, int i, double *bi, double *xij, void *adata);
extern void img_projS_jac(int j, int i, double *bi, double *Bij, void *adata);

// expert mode, full bundle adjustment
extern void img_projsRTS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
extern void img_projsRTS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);

// expert mode, camera motion only
extern void img_projsRT_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
extern void img_projsRT_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);

// expert mode, structure parameters only
extern void img_projsS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
extern void img_projsS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);




// simple mode, for camera with no distortion? 
extern void img_projKRTS(int j, int i, double *aj, double *bi, double *xij, void *adata);
extern void img_projKRTS_jac(int j, int i, double *aj, double *bi, double *Aij, double *Bij, void *adata);
extern void img_projKRT(int j, int i, double *aj, double *xij, void *adata);
extern void img_projKRT_jac(int j, int i, double *aj, double *Aij, void *adata);
extern void img_projKS(int j, int i, double *bi, double *xij, void *adata);
extern void img_projKS_jac(int j, int i, double *bi, double *Bij, void *adata);

// expert mode for camera with no distortion???
extern void img_projsKRTS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
extern void img_projsKRTS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);
extern void img_projsKRT_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
extern void img_projsKRT_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);
extern void img_projsKS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
extern void img_projsKS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);



// simple mode for cameras, structure and motion with distortion
extern void img_projKDRTS(int j, int i, double *aj, double *bi, double *xij, void *adata);
extern void img_projKDRTS_jac(int j, int i, double *aj, double *bi, double *Aij, double *Bij, void *adata);
extern void img_projKDRT(int j, int i, double *aj, double *xij, void *adata);
extern void img_projKDRT_jac(int j, int i, double *aj, double *Aij, void *adata);
extern void img_projKDS(int j, int i, double *bi, double *xij, void *adata);
extern void img_projKDS_jac(int j, int i, double *bi, double *Bij, void *adata);

// expert mode for cameras, structure and motion with distortion
extern void img_projsKDRTS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
extern void img_projsKDRTS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);
extern void img_projsKDRT_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
extern void img_projsKDRT_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);
extern void img_projsKDS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
extern void img_projsKDS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);



/* pointers to additional data, used for computed image projections and their jacobians */
struct globs_{
	double *rot0params; /* initial rotation parameters, combined with a local rotation parameterization */
	double *intrcalib; /* the 5 intrinsic calibration parameters in the order [fu, u0, v0, ar, skew],
                      * where ar is the aspect ratio fv/fu.
                      * Used only when calibration is fixed for all cameras;
                      * otherwise, it is null and the intrinsic parameters are
                      * included in the set of motion parameters for each camera
                      */
  int nccalib; /* number of calibration parameters that must be kept constant.
                * 0: all parameters are free 
                * 1: skew is fixed to its initial value, all other parameters vary (i.e. fu, u0, v0, ar) 
                * 2: skew and aspect ratio are fixed to their initial values, all other parameters vary (i.e. fu, u0, v0)
                * 3: meaningless
                * 4: skew, aspect ratio and principal point are fixed to their initial values, only the focal length varies (i.e. fu)
                * 5: all intrinsics are kept fixed to their initial values
                * >5: meaningless
                * Used only when calibration varies among cameras
                */

  int ncdist; /* number of distortion parameters in Bouguet's model that must be kept constant.
               * 0: all parameters are free 
               * 1: 6th order radial distortion term (kc[4]) is fixed
               * 2: 6th order radial distortion and one of the tangential distortion terms (kc[3]) are fixed
               * 3: 6th order radial distortion and both tangential distortion terms (kc[3], kc[2]) are fixed [i.e., only 2nd & 4th order radial dist.]
               * 4: 4th & 6th order radial distortion terms and both tangential distortion ones are fixed [i.e., only 2nd order radial dist.]
               * 5: all distortion parameters are kept fixed to their initial values
               * >5: meaningless
               * Used only when calibration varies among cameras and distortion is to be estimated
               */
  int cnp, pnp, mnp; /* dimensions */

	double *ptparams; /* needed only when bundle adjusting for camera parameters only */
	double *camparams; /* needed only when bundle adjusting for structure parameters only */
};







#endif /* _SBA_PROJS_H_ */
