#ifndef _CONSTANTS_H
#define _CONSTANTS_H 

/* unused grid parameters */
#define NRADII 100
#define NAZ    360

#define NZ     20
#define RADMAX 100*AU 
#define RADMIN 0.001*AU
#define ZMAX   10*AU
#define ZMIN   0.0001*AU
#define SUBPIXREGION 20 /*AU*/

/* physical constants */
#define PI 4*atan(1)
#define AU 1.49597871E13
#define DTOR PI/180.
#define NU 345.0E9
#define PLANCK 6.6260755E-27
#define CC    2.9979245810E10
#define GG    6.67259E-8
#define MP   1.6726219E-24
#define MSUN  1.989E33
#define MSTAR 1.0*MSUN
#define AD_GAMMA 1.400
#define KB   1.380658E-16
#define TOJY 1E23
#define DUST_SUB_T 1500.

/* image stuff */
#define IMSIZE 2.0 /*arcsec*/
#define NX     2048
#define DIST 145*AU


/* model starts */
#define SIG0 1E-5
#define RC   1E1*AU
#define GRAD 1
#define INC  45.0
#define PA   100.0
#define T10 45.0
#define QQ  0.5 

/* debug shiz */
#define DEBUG 1

#define NPARAMS 6
#endif
