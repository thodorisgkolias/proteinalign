/* sys/gsl_sys.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS
     
     double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS
     
#endif /* __GSL_SYS_H__ */
     
     
     /* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__
     
#include <limits.h>
#include <float.h>
     
     /* magic constants; mostly for the benefit of the implementation */
     
     /* -*-MACHINE CONSTANTS-*-
     *
     * PLATFORM: Whiz-O-Matic 9000
     * FP_PLATFORM: IEEE-Virtual
     * HOSTNAME: nnn.lanl.gov
     * DATE: Fri Nov 20 17:53:26 MST 1998
     */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)
     
#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)
     
#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02
     
#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)
     
#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)
     
#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01
     
#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)
     
     /* !MACHINE CONSTANTS! */
     
     
     /* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON
     
     
     
     /* Here are the constants related to or derived from
     * machine constants. These are not to be confused with
     * the constants that define various precision levels
     * for the precision/error system.
     *
     * This information is determined at configure time
     * and is platform dependent. Edit at your own risk.
     *
     * PLATFORM: WHIZ-O-MATIC
     * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
     * CONFIG-HOST: nnn.lanl.gov
     */
     
     /* machine precision constants */
     /* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)
     
     
#endif /* __GSL_MACHINE_H__ */
     
     /* gsl_types.h
     * 
     * Copyright (C) 2001 Brian Gough
     * 
     * This program is free software; you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation; either version 2 of the License, or (at
     * your option) any later version.
     * 
     * This program is distributed in the hope that it will be useful, but
     * WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     * General Public License for more details.
     * 
     * You should have received a copy of the GNU General Public License
     * along with this program; if not, write to the Free Software
     * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
     */
     
#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__
     
#ifndef GSL_VAR
     
#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif
     
#endif
     
#endif /* __GSL_TYPES_H__ */
     
     /* gsl_precision.h
* 
* Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
* 
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
     
     /* Author:  B. Gough and G. Jungman */
     
#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
     //#include <gsl/gsl_types.h>
     
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
     
     __BEGIN_DECLS
          
          
          /* A type for the precision indicator.
          * This is mainly for pedagogy.
          */
          typedef  unsigned int  gsl_prec_t;
     
     
     /* The number of precision types.
     * Remember that precision-mode
     * can index an array.
     */
#define _GSL_PREC_T_NUM 3
     
     
     /* Arrays containing derived
* precision constants for the
* different precision levels.
*/
     GSL_VAR const double gsl_prec_eps[];
     GSL_VAR const double gsl_prec_sqrt_eps[];
     GSL_VAR const double gsl_prec_root3_eps[];
     GSL_VAR const double gsl_prec_root4_eps[];
     GSL_VAR const double gsl_prec_root5_eps[];
     GSL_VAR const double gsl_prec_root6_eps[];
     
     
     __END_DECLS
          
#endif /* __GSL_PRECISION_H__ */
          
          
          /* gsl_nan.h
* 
* Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
* 
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
          
#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__
          
#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif
          
#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif
          
#define GSL_POSZERO (+0)
#define GSL_NEGZERO (-0)
          
#endif /* __GSL_NAN_H__ */
          
          /* gsl_pow_int.h
* 
* Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman, Brian Gough
* 
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
          
#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
          
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
          
          __BEGIN_DECLS
               
#ifdef HAVE_INLINE
               extern inline double gsl_pow_2(const double x);
               extern inline double gsl_pow_3(const double x);
               extern inline double gsl_pow_4(const double x);
               extern inline double gsl_pow_5(const double x);
               extern inline double gsl_pow_6(const double x);
               extern inline double gsl_pow_7(const double x);
               extern inline double gsl_pow_8(const double x);
               extern inline double gsl_pow_9(const double x);
               
               extern inline double gsl_pow_2(const double x) { return x*x;   }
               extern inline double gsl_pow_3(const double x) { return x*x*x; }
               extern inline double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
               extern inline double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
               extern inline double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
               extern inline double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
               extern inline double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
               extern inline double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#else
               double gsl_pow_2(const double x);
               double gsl_pow_3(const double x);
               double gsl_pow_4(const double x);
               double gsl_pow_5(const double x);
               double gsl_pow_6(const double x);
               double gsl_pow_7(const double x);
               double gsl_pow_8(const double x);
               double gsl_pow_9(const double x);
#endif
               
               double gsl_pow_int(double x, int n);
               
               __END_DECLS
                    
#endif /* __GSL_POW_INT_H__ */
                    
                    /* gsl_math.h
* 
* Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman, Brian Gough
* 
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
                    
#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
                    //#include <gsl/gsl_sys.h>
                    //#include <gsl/gsl_machine.h>
                    //#include <gsl/gsl_precision.h>
                    //#include <gsl/gsl_nan.h>
                    //#include <gsl/gsl_pow_int.h>
                    
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif
                    
#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif
                    
#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif
                    
#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif
                    
#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif
                    
                    
#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif
                    
#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif
                    
#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif
                    
#ifndef M_PI_4
#define M_PI_4     0.78539816339744830966156608458      /* pi/4 */
#endif
                    
#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif
                    
#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif
                    
#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif
                    
#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif
                    
#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif
                    
#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif
                    
#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif
                    
#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif
                    
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
                    
                    __BEGIN_DECLS
                         
                         /* other needlessly compulsive abstractions */
                         
#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)
                         
                         /* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))
                         
                         /* Define MAX and MIN macros/functions if they don't exist. */
                         
                         /* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))
                         
                         /* function versions of the above, in case they are needed */
                         double gsl_max (double a, double b);
                         double gsl_min (double a, double b);
                         
                         /* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE
                         
                         extern inline int GSL_MAX_INT (int a, int b);
                         extern inline int GSL_MIN_INT (int a, int b);
                         extern inline double GSL_MAX_DBL (double a, double b);
                         extern inline double GSL_MIN_DBL (double a, double b);
                         extern inline long double GSL_MAX_LDBL (long double a, long double b);
                         extern inline long double GSL_MIN_LDBL (long double a, long double b);
                         
                         extern inline int
                              GSL_MAX_INT (int a, int b)
                              {
                                   return GSL_MAX (a, b);
                              }
                         
                         extern inline int
                              GSL_MIN_INT (int a, int b)
                              {
                                   return GSL_MIN (a, b);
                              }
                         
                         extern inline double
                              GSL_MAX_DBL (double a, double b)
                              {
                                   return GSL_MAX (a, b);
                              }
                         
                         extern inline double
                              GSL_MIN_DBL (double a, double b)
                              {
                                   return GSL_MIN (a, b);
                              }
                         
                         extern inline long double
                              GSL_MAX_LDBL (long double a, long double b)
                              {
                                   return GSL_MAX (a, b);
                              }
                         
                         extern inline long double
                              GSL_MIN_LDBL (long double a, long double b)
                              {
                                   return GSL_MIN (a, b);
                              }
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */
                         
                         /* Definition of an arbitrary function with parameters */
                         
                         struct gsl_function_struct 
                         {
                              double (* function) (double x, void * params);
                              void * params;
                         };
                         
                         typedef struct gsl_function_struct gsl_function ;
                         
#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)
                         
                         /* Definition of an arbitrary function returning two values, r1, r2 */
                         
                         struct gsl_function_fdf_struct 
                         {
                              double (* f) (double x, void * params);
                              double (* df) (double x, void * params);
                              void (* fdf) (double x, void * params, double * f, double * df);
                              void * params;
                         };
                         
                         typedef struct gsl_function_fdf_struct gsl_function_fdf ;
                         
#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))
                         
                         
                         /* Definition of an arbitrary vector-valued function with parameters */
                         
                         struct gsl_function_vec_struct 
                         {
                              int (* function) (double x, double y[], void * params);
                              void * params;
                         };
                         
                         typedef struct gsl_function_vec_struct gsl_function_vec ;
                         
#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)
                         
                         __END_DECLS
                              
#endif /* __GSL_MATH_H__ */
                              

                              /* err/gsl_errno.h
                               * 
                               * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
                               * 
                               * This program is free software; you can redistribute it and/or modify
                               * it under the terms of the GNU General Public License as published by
                               * the Free Software Foundation; either version 2 of the License, or (at
                               * your option) any later version.
                               * 
                               * This program is distributed in the hope that it will be useful, but
                               * WITHOUT ANY WARRANTY; without even the implied warranty of
                               * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
                               * General Public License for more details.
                               * 
                               * You should have received a copy of the GNU General Public License
                               * along with this program; if not, write to the Free Software
                               * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
                               */
                              
#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__
                              
#include <stdio.h>
#include <errno.h>
                              //#include <gsl/gsl_types.h>
                              
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
                              
                              __BEGIN_DECLS
                                   
                                   enum { 
                                   GSL_SUCCESS  = 0, 
                                        GSL_FAILURE  = -1,
                                        GSL_CONTINUE = -2,  /* iteration has not converged */
                              GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
                              GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
                              GSL_EFAULT   = 3,   /* invalid pointer */
                              GSL_EINVAL   = 4,   /* invalid argument supplied by user */
                              GSL_EFAILED  = 5,   /* generic failure */
                              GSL_EFACTOR  = 6,   /* factorization failed */
                              GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
                              GSL_ENOMEM   = 8,   /* malloc failed */
                              GSL_EBADFUNC = 9,   /* problem with user-supplied function */
                              GSL_ERUNAWAY = 10,  /* iterative process is out of control */
                              GSL_EMAXITER = 11,  /* exceeded max number of iterations */
                              GSL_EZERODIV = 12,  /* tried to divide by zero */
                              GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
                              GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
                              GSL_EUNDRFLW = 15,  /* underflow */
                              GSL_EOVRFLW  = 16,  /* overflow  */
                              GSL_ELOSS    = 17,  /* loss of accuracy */
                              GSL_EROUND   = 18,  /* failed because of roundoff error */
                              GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
                              GSL_ENOTSQR  = 20,  /* matrix not square */
                              GSL_ESING    = 21,  /* apparent singularity detected */
                              GSL_EDIVERGE = 22,  /* integral or series is divergent */
                              GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
                              GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
                              GSL_ECACHE   = 25,  /* cache limit exceeded */
                              GSL_ETABLE   = 26,  /* table limit exceeded */
                              GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
                              GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
                              GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
                              GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
                              GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
                              GSL_EOF      = 32   /* end of file */
                              } ;
                              
                              void gsl_error (const char * reason, const char * file, int line,
                                              int gsl_errno);
                              
                              void gsl_stream_printf (const char *label, const char *file,
                                                      int line, const char *reason);
                              
                              const char * gsl_strerror (const int gsl_errno);
                              
                              typedef void gsl_error_handler_t (const char * reason, const char * file,
                                                                int line, int gsl_errno);
                              
                              typedef void gsl_stream_handler_t (const char * label, const char * file,
                                                                 int line, const char * reason);
                              
                              gsl_error_handler_t * 
                                   gsl_set_error_handler (gsl_error_handler_t * new_handler);
                              
                              gsl_error_handler_t *
                                   gsl_set_error_handler_off (void);
                              
                              gsl_stream_handler_t * 
                                   gsl_set_stream_handler (gsl_stream_handler_t * new_handler);
                              
                              FILE * gsl_set_stream (FILE * new_stream);
                              
                              /* GSL_ERROR: call the error handler, and return the error code */
                              
#define GSL_ERROR(reason, gsl_errno)                                                   \
                              do {                                                     \
                                   gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
                                   return gsl_errno ;                                  \
                              } while (0)
                         
                         /* GSL_ERROR_VAL: call the error handler, and return the given value */
                         
#define GSL_ERROR_VAL(reason, gsl_errno, value)                                   \
                         do {                                                     \
                              gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
                              return value ;                                      \
                         } while (0)
                    
                    /* GSL_ERROR_VOID: call the error handler, and then return
                              (for void functions which still need to generate an error) */
                    
#define GSL_ERROR_VOID(reason, gsl_errno)                                    \
                    do {                                                     \
                         gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
                         return ;                                            \
                    } while (0)
               
               /* GSL_ERROR_NULL suitable for out-of-memory conditions */
               
#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)
               
               /* Sometimes you have several status results returned from
               * function calls and you want to combine them in some sensible
               * way. You cannot produce a "total" status condition, but you can
               * pick one from a set of conditions based on an implied hierarchy.
               *
               * In other words:
               *    you have: status_a, status_b, ...
               *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
               *
               * In this example you consider status_a to be more important and
               * it is checked first, followed by the others in the order specified.
               *
               * Here are some dumb macros to do this.
               */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))
               
#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)
               
               __END_DECLS
                    
#endif /* __GSL_ERRNO_H__ */
                    
                    
                    /* err/gsl_message.h
* 
* Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
* 
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
                    
#ifndef __GSL_MESSAGE_H__
#define __GSL_MESSAGE_H__
                    //#include <gsl/gsl_types.h>
                    
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif
                    
                    __BEGIN_DECLS
                         
                         /* Provide a general messaging service for client use.  Messages can
                         * be selectively turned off at compile time by defining an
                         * appropriate message mask. Client code which uses the GSL_MESSAGE()
                         * macro must provide a mask which is or'ed with the GSL_MESSAGE_MASK.
                         *
                         * The messaging service can be completely turned off
                         * by defining GSL_MESSAGING_OFF.  */
                         
                         void gsl_message(const char * message, const char * file, int line,
                                          unsigned int mask);
                    
#ifndef GSL_MESSAGE_MASK
#define GSL_MESSAGE_MASK 0xffffffffu /* default all messages allowed */
#endif
                    
                    GSL_VAR unsigned int gsl_message_mask ;
                    
                    /* Provide some symolic masks for client ease of use. */
                    
                    enum {
                         GSL_MESSAGE_MASK_A = 1,
                         GSL_MESSAGE_MASK_B = 2,
                         GSL_MESSAGE_MASK_C = 4,
                         GSL_MESSAGE_MASK_D = 8,
                         GSL_MESSAGE_MASK_E = 16,
                         GSL_MESSAGE_MASK_F = 32,
                         GSL_MESSAGE_MASK_G = 64,
                         GSL_MESSAGE_MASK_H = 128
                    } ;
                    
#ifdef GSL_MESSAGING_OFF        /* throw away messages */ 
#define GSL_MESSAGE(message, mask) do { } while(0)
#else                           /* output all messages */
#define GSL_MESSAGE(message, mask)                                              \
                    do {                                                        \
                         if (mask & GSL_MESSAGE_MASK)                           \
                              gsl_message (message, __FILE__, __LINE__, mask) ; \
                    } while (0)
#endif
               
               __END_DECLS
                         
#endif /* __GSL_MESSAGE_H__ */
                         
                         /* err/stream.h
* 
* Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
* 
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
                         
                         //#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
                         
                         //#include <gsl/gsl_errno.h>
                         //#include <gsl/gsl_message.h>
                         
                         FILE * gsl_stream = NULL ;
                         gsl_stream_handler_t * gsl_stream_handler = NULL;
                         
                         void
                              gsl_stream_printf (const char *label, const char *file, int line, 
                                                 const char *reason)
                              {
                                   if (gsl_stream == NULL)
                                   {
                                        //    gsl_stream = stderr;
                                   }
                                   if (gsl_stream_handler)
                                   {
                                        (*gsl_stream_handler) (label, file, line, reason);
                                        return;
                                   }
                                   fprintf (gsl_stream, "gsl: %s:%d: %s: %s\n", file, line, label, reason);
                                   
                              }
                         
                         gsl_stream_handler_t *
                              gsl_set_stream_handler (gsl_stream_handler_t * new_handler)
                              {
                                   gsl_stream_handler_t * previous_handler = gsl_stream_handler;
                                   gsl_stream_handler = new_handler;
                                   return previous_handler;
                              }
                         
                         FILE *
                              gsl_set_stream (FILE * new_stream)
                              {
                                   FILE * previous_stream;
                                   if (gsl_stream == NULL) {
                                        //   gsl_stream = stderr;
                                   }
                                   previous_stream = gsl_stream;
                                   gsl_stream = new_stream;
                                   return previous_stream;
                              }
                         
                         
                         
                         /* err/error.c
                         * 
                         * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
                         * 
                         * This program is free software; you can redistribute it and/or modify
                         * it under the terms of the GNU General Public License as published by
                         * the Free Software Foundation; either version 2 of the License, or (at
                         * your option) any later version.
                         * 
                         * This program is distributed in the hope that it will be useful, but
                         * WITHOUT ANY WARRANTY; without even the implied warranty of
                         * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
                         * General Public License for more details.
                         * 
                         * You should have received a copy of the GNU General Public License
                         * along with this program; if not, write to the Free Software
                         * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
                         */
                         
                         //#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
                         
                         //#include <gsl/gsl_errno.h>
                         //#include <gsl/gsl_message.h>
                         
                         gsl_error_handler_t * gsl_error_handler = NULL;
                         
                         static void no_error_handler (const char *reason, const char *file, int line, int gsl_errno);
                         
                         void
                              gsl_error (const char * reason, const char * file, int line, int gsl_errno)
                              {
                                   if (gsl_error_handler) 
                                   {
                                        (*gsl_error_handler) (reason, file, line, gsl_errno);
                                        return ;
                                   }
                                   
                                   gsl_stream_printf ("ERROR", file, line, reason);
                                   
                                   //  fflush (stdout);
                                   // fprintf (stderr, "Default GSL error handler invoked.\n");
                                   // fflush (stderr);
                                   
                                   // abort ();
                              }
                         
                         gsl_error_handler_t *
                              gsl_set_error_handler (gsl_error_handler_t * new_handler)
                              {
                                   gsl_error_handler_t * previous_handler = gsl_error_handler;
                                   gsl_error_handler = new_handler;
                                   return previous_handler;
                              }
                         
                         
                         gsl_error_handler_t *
                              gsl_set_error_handler_off (void)
                              {
                                   gsl_error_handler_t * previous_handler = gsl_error_handler;
                                   gsl_error_handler = no_error_handler;
                                   return previous_handler;
                              }
                         
                         static void
no_error_handler (const char *reason, const char *file, int line, int gsl_errno)
                              {
                                   /* do nothing */
                                   reason = 0;
                                   file = 0;
                                   line = 0;
                                   gsl_errno = 0;
                                   return;
                              }
                         
                         
                         

/* integration/gsl_integration.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef __GSL_INTEGRATION_H__
#define __GSL_INTEGRATION_H__
#include <stdlib.h>
//#include "func.h"
//#include "error.h"
//#include "qag.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* Workspace for adaptive integrators */

typedef struct
  {
    size_t limit;
    size_t size;
    size_t nrmax;
    size_t i;
    size_t maximum_level;
    double *alist;
    double *blist;
    double *rlist;
    double *elist;
    size_t *order;
    size_t *level;
  }
gsl_integration_workspace;

gsl_integration_workspace *
  gsl_integration_workspace_alloc (const size_t n);

void
  gsl_integration_workspace_free (gsl_integration_workspace * w);


/* Workspace for QAWS integrator */

typedef struct
{
  double alpha;
  double beta;
  int mu;
  int nu;
  double ri[25];
  double rj[25];
  double rg[25];
  double rh[25];
}
gsl_integration_qaws_table;

gsl_integration_qaws_table * 
gsl_integration_qaws_table_alloc (double alpha, double beta, int mu, int nu);

int
gsl_integration_qaws_table_set (gsl_integration_qaws_table * t,
                                double alpha, double beta, int mu, int nu);

void
gsl_integration_qaws_table_free (gsl_integration_qaws_table * t);

/* Workspace for QAWO integrator */

enum gsl_integration_qawo_enum { GSL_INTEG_COSINE, GSL_INTEG_SINE };

typedef struct
{
  size_t n;
  double omega;
  double L;
  double par;
  enum gsl_integration_qawo_enum sine;
  double *chebmo;
}
gsl_integration_qawo_table;

gsl_integration_qawo_table * 
gsl_integration_qawo_table_alloc (double omega, double L, 
                                  enum gsl_integration_qawo_enum sine,
                                  size_t n);

int
gsl_integration_qawo_table_set (gsl_integration_qawo_table * t,
                                double omega, double L,
                                enum gsl_integration_qawo_enum sine);

int
gsl_integration_qawo_table_set_length (gsl_integration_qawo_table * t,
                                       double L);

void
gsl_integration_qawo_table_free (gsl_integration_qawo_table * t);


/* Definition of an integration rule */

typedef void gsl_integration_rule (const gsl_function * f,
                                   double a, double b,
                                   double *result, double *abserr,
                                   double *defabs, double *resabs);

void gsl_integration_qk15 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk21 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk31 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk41 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk51 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qk61 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

void gsl_integration_qcheb (gsl_function * f, double a, double b, 
                            double *cheb12, double *cheb24);

/* The low-level integration rules in QUADPACK are identified by small
   integers (1-6). We'll use symbolic constants to refer to them.  */

enum
  {
    GSL_INTEG_GAUSS15 = 1,      /* 15 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS21 = 2,      /* 21 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS31 = 3,      /* 31 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS41 = 4,      /* 41 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS51 = 5,      /* 51 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS61 = 6       /* 61 point Gauss-Kronrod rule */
  };

void 
gsl_integration_qk (const int n, const double xgk[], 
                    const double wg[], const double wgk[],
                    double fv1[], double fv2[],
                    const gsl_function *f, double a, double b,
                    double * result, double * abserr, 
                    double * resabs, double * resasc);


int gsl_integration_qng (const gsl_function * f,
                         double a, double b,
                         double epsabs, double epsrel,
                         double *result, double *abserr,
                         size_t * neval);

int gsl_integration_qag (const gsl_function * f,
                         double a, double b,
                         double epsabs, double epsrel, size_t limit,
                         int key,
                         gsl_integration_workspace * workspace,
                         double *result, double *abserr);

int gsl_integration_qagi (gsl_function * f,
                          double epsabs, double epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          double *result, double *abserr);

int gsl_integration_qagiu (gsl_function * f,
                           double a,
                           double epsabs, double epsrel, size_t limit,
                           gsl_integration_workspace * workspace,
                           double *result, double *abserr);

int gsl_integration_qagil (gsl_function * f,
                           double b,
                           double epsabs, double epsrel, size_t limit,
                           gsl_integration_workspace * workspace,
                           double *result, double *abserr);


int gsl_integration_qags (const gsl_function * f,
                          double a, double b,
                          double epsabs, double epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          double *result, double *abserr);

int gsl_integration_qagp (const gsl_function * f,
                          double *pts, size_t npts,
                          double epsabs, double epsrel, size_t limit,
                          gsl_integration_workspace * workspace,
                          double *result, double *abserr);

int gsl_integration_qawc (gsl_function *f,
                          const double a, const double b, const double c,
                          const double epsabs, const double epsrel, const size_t limit,
                          gsl_integration_workspace * workspace,
                          double * result, double * abserr);

int gsl_integration_qaws (gsl_function * f,
                          const double a, const double b,
                          gsl_integration_qaws_table * t,
                          const double epsabs, const double epsrel,
                          const size_t limit,
                          gsl_integration_workspace * workspace,
                          double *result, double *abserr);

int gsl_integration_qawo (gsl_function * f,
                          const double a,
                          const double epsabs, const double epsrel,
                          const size_t limit,
                          gsl_integration_workspace * workspace,
                          gsl_integration_qawo_table * wf,
                          double *result, double *abserr);

int gsl_integration_qawf (gsl_function * f,
                          const double a,
                          const double epsabs,
                          const size_t limit,
                          gsl_integration_workspace * workspace,
                          gsl_integration_workspace * cycle_workspace,
                          gsl_integration_qawo_table * wf,
                          double *result, double *abserr);

__END_DECLS

#endif /* __GSL_INTEGRATION_H__ */




/* integration/workspace.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//#include <config.h>
#include <stdlib.h>
//#include <gsl/gsl_integration.h>
//#include <gsl/gsl_errno.h>

gsl_integration_workspace *
gsl_integration_workspace_alloc (const size_t n) 
{
  gsl_integration_workspace * w ;
  
  if (n == 0)
    {
      GSL_ERROR_VAL ("workspace length n must be positive integer",
                        GSL_EDOM, 0);
    }

  w = (gsl_integration_workspace *) 
    malloc (sizeof (gsl_integration_workspace));

  if (w == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for workspace struct",
                        GSL_ENOMEM, 0);
    }

  w->alist = (double *) malloc (n * sizeof (double));

  if (w->alist == 0)
    {
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for alist ranges",
                        GSL_ENOMEM, 0);
    }

  w->blist = (double *) malloc (n * sizeof (double));

  if (w->blist == 0)
    {
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for blist ranges",
                        GSL_ENOMEM, 0);
    }

  w->rlist = (double *) malloc (n * sizeof (double));

  if (w->rlist == 0)
    {
      free (w->blist);
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for rlist ranges",
                        GSL_ENOMEM, 0);
    }


  w->elist = (double *) malloc (n * sizeof (double));

  if (w->elist == 0)
    {
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for elist ranges",
                        GSL_ENOMEM, 0);
    }

  w->order = (size_t *) malloc (n * sizeof (size_t));

  if (w->order == 0)
    {
      free (w->elist);
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for order ranges",
                        GSL_ENOMEM, 0);
    }

  w->level = (size_t *) malloc (n * sizeof (size_t));

  if (w->level == 0)
    {
      free (w->order);
      free (w->elist);
      free (w->rlist);
      free (w->blist);
      free (w->alist);
      free (w);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for order ranges",
                        GSL_ENOMEM, 0);
    }

  w->size = 0 ;
  w->limit = n ;
  w->maximum_level = 0 ;
  
  return w ;
}

void
gsl_integration_workspace_free (gsl_integration_workspace * w)
{
  free (w->level) ;
  free (w->order) ;
  free (w->elist) ;
  free (w->rlist) ;
  free (w->blist) ;
  free (w->alist) ;
  free (w) ;
}

/*
size_t 
gsl_integration_workspace_limit (gsl_integration_workspace * w) 
{
  return w->limit ;
}


size_t 
gsl_integration_workspace_size (gsl_integration_workspace * w) 
{
  return w->size ;
}
*/

/* END__WORKSPACE_H__ */


/* sys/coerce.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//#include <config.h>
//#include <math.h>

double gsl_coerce_double (const double x);

double 
gsl_coerce_double (const double x)
{
  volatile double y;
  y = x;
  return y;
}

float gsl_coerce_float (const float x);

float 
gsl_coerce_float (const float x)
{
  volatile float y;
  y = x;
  return y;
}

/* The following function is not needed, but is included for completeness */

long double gsl_coerce_long_double (const long double x);

long double 
gsl_coerce_long_double (const long double x)
{
  volatile long double y;
  y = x;
  return y;
}

/* END_COERCE_H */

/* integration/err.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//#include <config.h>
#include <math.h>
#include <float.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_integration.h>

static double rescale_error (double err, const double result_abs, const double result_asc) ;

static double
rescale_error (double err, const double result_abs, const double result_asc)
{
  err = fabs(err) ;

  if (result_asc != 0 && err != 0)
      {
        double scale = pow((200 * err / result_asc), 1.5) ;
        
        if (scale < 1)
          {
            err = result_asc * scale ;
          }
        else 
          {
            err = result_asc ;
          }
      }
  if (result_abs > GSL_DBL_MIN / (50 * GSL_DBL_EPSILON))
    {
      double min_err = 50 * GSL_DBL_EPSILON * result_abs ;

      if (min_err > err) 
        {
          err = min_err ;
        }
    }
  
  return err ;
}

/* END__ERR_H__ */

/* integration/qag.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//#include <config.h>
#include <stdlib.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_integration.h>
/* integration/initialise.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

static inline
     void initialise (gsl_integration_workspace * workspace, double a, double b);

static inline
     void initialise (gsl_integration_workspace * workspace, double a, double b)
     {
          workspace->size = 0;
          workspace->nrmax = 0;
          workspace->i = 0;
          workspace->alist[0] = a;
          workspace->blist[0] = b;
          workspace->rlist[0] = 0.0;
          workspace->elist[0] = 0.0;
          workspace->order[0] = 0;
          workspace->level[0] = 0;
          
          workspace->maximum_level = 0;
     }

//#include "initialise.h"

static inline
     void set_initial_result (gsl_integration_workspace * workspace, 
                              double result, double error);

static inline
     void set_initial_result (gsl_integration_workspace * workspace, 
                              double result, double error)
     {
          workspace->size = 1;
          workspace->rlist[0] = result;
          workspace->elist[0] = error;
     }

//#include "set_initial.h"

/* integration/qpsrt.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

static inline void 
     qpsrt (gsl_integration_workspace * workspace);

static inline
     void qpsrt (gsl_integration_workspace * workspace)
     {
          const size_t last = workspace->size - 1;
          const size_t limit = workspace->limit;
          
          double * elist = workspace->elist;
          size_t * order = workspace->order;
          
          double errmax ;
          double errmin ;
          int i, k, top;
          
          size_t i_nrmax = workspace->nrmax;
          size_t i_maxerr = order[i_nrmax] ;
          
          /* Check whether the list contains more than two error estimates */
          
          if (last < 2) 
          {
               order[0] = 0 ;
               order[1] = 1 ;
               workspace->i = i_maxerr ;
               return ;
          }
          
          errmax = elist[i_maxerr] ;
          
          /* This part of the routine is only executed if, due to a difficult
          integrand, subdivision increased the error estimate. In the normal
          case the insert procedure should start after the nrmax-th largest
          error estimate. */
          
          while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]]) 
          {
               order[i_nrmax] = order[i_nrmax - 1] ;
               i_nrmax-- ;
          } 
          
          /* Compute the number of elements in the list to be maintained in
          descending order. This number depends on the number of
          subdivisions still allowed. */
          
          if(last < (limit/2 + 2)) 
          {
               top = last ;
          }
          else
          {
               top = limit - last + 1;
          }
          
          /* Insert errmax by traversing the list top-down, starting
          comparison from the element elist(order(i_nrmax+1)). */
          
          i = i_nrmax + 1 ;
          
          /* The order of the tests in the following line is important to
          prevent a segmentation fault */
          
          while (i < top && errmax < elist[order[i]])
          {
               order[i-1] = order[i] ;
               i++ ;
          }
          
          order[i-1] = i_maxerr ;
          
          /* Insert errmin by traversing the list bottom-up */
          
          errmin = elist[last] ;
          
          k = top - 1 ;
          
          while (k > i - 2 && errmin >= elist[order[k]])
          {
               order[k+1] = order[k] ;
               k-- ;
          }
          
          order[k+1] = last ;
          
          /* Set i_max and e_max */
          
          i_maxerr = order[i_nrmax] ;
          
          workspace->i = i_maxerr ;
          workspace->nrmax = i_nrmax ;
     }



//#include "qpsrt.h"

/* integration/util.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

static inline
     void update (gsl_integration_workspace * workspace,
                  double a1, double b1, double area1, double error1,
                  double a2, double b2, double area2, double error2);

static inline void
     retrieve (const gsl_integration_workspace * workspace, 
               double * a, double * b, double * r, double * e);



static inline
     void update (gsl_integration_workspace * workspace,
                  double a1, double b1, double area1, double error1,
                  double a2, double b2, double area2, double error2)
     {
          double * alist = workspace->alist ;
          double * blist = workspace->blist ;
          double * rlist = workspace->rlist ;
          double * elist = workspace->elist ;
          size_t * level = workspace->level ;
          
          const size_t i_max = workspace->i ;
          const size_t i_new = workspace->size ;
          
          const size_t new_level = workspace->level[i_max] + 1;
          
          /* append the newly-created intervals to the list */
          
          if (error2 > error1)
          {
               alist[i_max] = a2;        /* blist[maxerr] is already == b2 */
          rlist[i_max] = area2;
          elist[i_max] = error2;
          level[i_max] = new_level;
          
          alist[i_new] = a1;
          blist[i_new] = b1;
          rlist[i_new] = area1;
          elist[i_new] = error1;
          level[i_new] = new_level;
          }
          else
          {
               blist[i_max] = b1;        /* alist[maxerr] is already == a1 */
          rlist[i_max] = area1;
          elist[i_max] = error1;
          level[i_max] = new_level;
          
          alist[i_new] = a2;
          blist[i_new] = b2;
          rlist[i_new] = area2;
          elist[i_new] = error2;
          level[i_new] = new_level;
          }
          
          workspace->size++;
          
          if (new_level > workspace->maximum_level)
          {
               workspace->maximum_level = new_level;
          }
          
          qpsrt (workspace) ;
     }

static inline void
     retrieve (const gsl_integration_workspace * workspace, 
               double * a, double * b, double * r, double * e)
     {
          const size_t i = workspace->i;
          double * alist = workspace->alist;
          double * blist = workspace->blist;
          double * rlist = workspace->rlist;
          double * elist = workspace->elist;
          
          *a = alist[i] ;
          *b = blist[i] ;
          *r = rlist[i] ;
          *e = elist[i] ;
     }

static inline double
     sum_results (const gsl_integration_workspace * workspace);

static inline double
     sum_results (const gsl_integration_workspace * workspace)
     {
          const double * const rlist = workspace->rlist ;
          const size_t n = workspace->size;
          
          size_t k;
          double result_sum = 0;
          
          for (k = 0; k < n; k++)
          {
               result_sum += rlist[k];
          }
          
          return result_sum;
     }

static inline int
     subinterval_too_small (double a1, double a2, double b2);

static inline int
     subinterval_too_small (double a1, double a2, double b2)
     {
          const double e = GSL_DBL_EPSILON;
          const double u = GSL_DBL_MIN;
          
          double tmp = (1 + 100 * e) * (fabs (a2) + 1000 * u);
          
          int status = fabs (a1) <= tmp && fabs (b2) <= tmp;
          
          return status;
     }


//#include "util.h"

static int
qag (const gsl_function *f,
     const double a, const double b,
     const double epsabs, const double epsrel,
     const size_t limit,
     gsl_integration_workspace * workspace,
     double * result, double * abserr,
     gsl_integration_rule * q) ;

int
gsl_integration_qag (const gsl_function *f,
                     double a, double b,
                     double epsabs, double epsrel, size_t limit,
                     int key,
                     gsl_integration_workspace * workspace,
                     double * result, double * abserr)
{
  int status ;
  gsl_integration_rule * integration_rule = gsl_integration_qk61 ;

 


  status = qag (f, a, b, epsabs, epsrel, limit,
                workspace, 
                result, abserr, 
                integration_rule) ;
  
  return status ;
}

static int
qag (const gsl_function * f,
     const double a, const double b,
     const double epsabs, const double epsrel,
     const size_t limit,
     gsl_integration_workspace * workspace,
     double *result, double *abserr,
     gsl_integration_rule * q)
{
  double area, errsum;
  double result0, abserr0, resabs0, resasc0;
  double tolerance;
  size_t iteration = 0;
  int roundoff_type1 = 0, roundoff_type2 = 0, error_type = 0;

  double round_off;     

  /* Initialize results */

  initialise (workspace, a, b);

  *result = 0;
  *abserr = 0;

  if (limit > workspace->limit)
    {
      GSL_ERROR ("iteration limit exceeds available workspace", GSL_EINVAL) ;
    }

  if (epsabs <= 0 && (epsrel < 50 * GSL_DBL_EPSILON || epsrel < 0.5e-28))
    {
      GSL_ERROR ("tolerance cannot be acheived with given epsabs and epsrel",
                 GSL_EBADTOL);
    }

  /* perform the first integration */

  q (f, a, b, &result0, &abserr0, &resabs0, &resasc0);

  set_initial_result (workspace, result0, abserr0);

  /* Test on accuracy */

  tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (result0));

  /* need IEEE rounding here to match original quadpack behavior */

  //round_off = GSL_COERCE_DBL (50 * GSL_DBL_EPSILON * resabs0);
  /* Cannot use GSL_COERSE_DBL */

  //round_off = 50 * GSL_DBL_EPSILON * resabs0;

  round_off = gsl_coerce_double(50 * GSL_DBL_EPSILON * resabs0);

  if (abserr0 <= round_off && abserr0 > tolerance)
    {
      *result = result0;
      *abserr = abserr0;

      GSL_ERROR ("cannot reach tolerance because of roundoff error "
                 "on first attempt", GSL_EROUND);
    }
  else if ((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0.0)
    {
      *result = result0;
      *abserr = abserr0;

      return GSL_SUCCESS;
    }
  else if (limit == 1)
    {
      *result = result0;
      *abserr = abserr0;

      GSL_ERROR ("a maximum of one iteration was insufficient", GSL_EMAXITER);
    }

  area = result0;
  errsum = abserr0;

  iteration = 1;

  do
    {
      double a1, b1, a2, b2;
      double a_i, b_i, r_i, e_i;
      double area1 = 0, area2 = 0, area12 = 0;
      double error1 = 0, error2 = 0, error12 = 0;
      double resasc1, resasc2;
      double resabs1, resabs2;

      /* Bisect the subinterval with the largest error estimate */

      retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

      a1 = a_i; 
      b1 = 0.5 * (a_i + b_i);
      a2 = b1;
      b2 = b_i;

      q (f, a1, b1, &area1, &error1, &resabs1, &resasc1);
      q (f, a2, b2, &area2, &error2, &resabs2, &resasc2);

      area12 = area1 + area2;
      error12 = error1 + error2;

      errsum += (error12 - e_i);
      area += area12 - r_i;

      if (resasc1 != error1 && resasc2 != error2)
        {
          double delta = r_i - area12;

          if (fabs (delta) <= 1.0e-5 * fabs (area12) && error12 >= 0.99 * e_i)
            {
              roundoff_type1++;
            }
          if (iteration >= 10 && error12 > e_i)
            {
              roundoff_type2++;
            }
        }

      tolerance = GSL_MAX_DBL (epsabs, epsrel * fabs (area));

      if (errsum > tolerance)
        {
          if (roundoff_type1 >= 6 || roundoff_type2 >= 20)
            {
              error_type = 2;   /* round off error */
            }

          /* set error flag in the case of bad integrand behaviour at
             a point of the integration range */

          if (subinterval_too_small (a1, a2, b2))
            {
              error_type = 3;
            }
        }

      update (workspace, a1, b1, area1, error1, a2, b2, area2, error2);

      retrieve (workspace, &a_i, &b_i, &r_i, &e_i);

      iteration++;

    }
  while (iteration < limit && !error_type && errsum > tolerance);

  *result = sum_results (workspace);
  *abserr = errsum;

  if (errsum <= tolerance)
    {
      return GSL_SUCCESS;
    }
  else if (error_type == 2)
    {
      GSL_ERROR ("roundoff error prevents tolerance from being achieved",
                 GSL_EROUND);
    }
  else if (error_type == 3)
    {
      GSL_ERROR ("bad integrand behavior found in the integration interval",
                 GSL_ESING);
    }
  else if (iteration == limit)
    {
      GSL_ERROR ("maximum number of subdivisions reached", GSL_EMAXITER);
    }
  else
    {
      GSL_ERROR ("could not integrate function", GSL_EFAILED);
    }
}

/* END_QAG_H */


/* integration/qk61.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//#include <config.h>
//#include <gsl/gsl_integration.h>

/* Gauss quadrature weights and kronrod quadrature abscissae and
   weights as evaluated with 80 decimal digit arithmetic by
   L. W. Fullerton, Bell Labs, Nov. 1981. */

static const double xgk[31] =   /* abscissae of the 61-point kronrod rule */
{
  0.999484410050490637571325895705811,
  0.996893484074649540271630050918695,
  0.991630996870404594858628366109486,
  0.983668123279747209970032581605663,
  0.973116322501126268374693868423707,
  0.960021864968307512216871025581798,
  0.944374444748559979415831324037439,
  0.926200047429274325879324277080474,
  0.905573307699907798546522558925958,
  0.882560535792052681543116462530226,
  0.857205233546061098958658510658944,
  0.829565762382768397442898119732502,
  0.799727835821839083013668942322683,
  0.767777432104826194917977340974503,
  0.733790062453226804726171131369528,
  0.697850494793315796932292388026640,
  0.660061064126626961370053668149271,
  0.620526182989242861140477556431189,
  0.579345235826361691756024932172540,
  0.536624148142019899264169793311073,
  0.492480467861778574993693061207709,
  0.447033769538089176780609900322854,
  0.400401254830394392535476211542661,
  0.352704725530878113471037207089374,
  0.304073202273625077372677107199257,
  0.254636926167889846439805129817805,
  0.204525116682309891438957671002025,
  0.153869913608583546963794672743256,
  0.102806937966737030147096751318001,
  0.051471842555317695833025213166723,
  0.000000000000000000000000000000000
};

/* xgk[1], xgk[3], ... abscissae of the 30-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 30-point gauss rule */

static const double wg[15] =    /* weights of the 30-point gauss rule */
{
  0.007968192496166605615465883474674,
  0.018466468311090959142302131912047,
  0.028784707883323369349719179611292,
  0.038799192569627049596801936446348,
  0.048402672830594052902938140422808,
  0.057493156217619066481721689402056,
  0.065974229882180495128128515115962,
  0.073755974737705206268243850022191,
  0.080755895229420215354694938460530,
  0.086899787201082979802387530715126,
  0.092122522237786128717632707087619,
  0.096368737174644259639468626351810,
  0.099593420586795267062780282103569,
  0.101762389748405504596428952168554,
  0.102852652893558840341285636705415
};

static const double wgk[31] =   /* weights of the 61-point kronrod rule */
{
  0.001389013698677007624551591226760,
  0.003890461127099884051267201844516,
  0.006630703915931292173319826369750,
  0.009273279659517763428441146892024,
  0.011823015253496341742232898853251,
  0.014369729507045804812451432443580,
  0.016920889189053272627572289420322,
  0.019414141193942381173408951050128,
  0.021828035821609192297167485738339,
  0.024191162078080601365686370725232,
  0.026509954882333101610601709335075,
  0.028754048765041292843978785354334,
  0.030907257562387762472884252943092,
  0.032981447057483726031814191016854,
  0.034979338028060024137499670731468,
  0.036882364651821229223911065617136,
  0.038678945624727592950348651532281,
  0.040374538951535959111995279752468,
  0.041969810215164246147147541285970,
  0.043452539701356069316831728117073,
  0.044814800133162663192355551616723,
  0.046059238271006988116271735559374,
  0.047185546569299153945261478181099,
  0.048185861757087129140779492298305,
  0.049055434555029778887528165367238,
  0.049795683427074206357811569379942,
  0.050405921402782346840893085653585,
  0.050881795898749606492297473049805,
  0.051221547849258772170656282604944,
  0.051426128537459025933862879215781,
  0.051494729429451567558340433647099
};

void
gsl_integration_qk61 (const gsl_function * f, double a, double b,
                      double *result, double *abserr,
                      double *resabs, double *resasc)
{
  double fv1[31], fv2[31];
  gsl_integration_qk (31, xgk, wg, wgk, fv1, fv2, f, a, b, result, abserr, resabs, resasc);
}

/* END_QK61_H */


/* integration/qk.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//#include <config.h>
#include <float.h>
//#include <math.h>
//#include <gsl/gsl_integration.h>
//#include "err.c"

void
gsl_integration_qk (const int n, 
                    const double xgk[], const double wg[], const double wgk[],
                    double fv1[], double fv2[],
                    const gsl_function * f, double a, double b,
                    double *result, double *abserr,
                    double *resabs, double *resasc)
{

  const double center = 0.5 * (a + b);
  const double half_length = 0.5 * (b - a);
  const double abs_half_length = fabs (half_length);
  const double f_center = GSL_FN_EVAL (f, center);

  double result_gauss = 0;
  double result_kronrod = f_center * wgk[n - 1];

  double result_abs = fabs (result_kronrod);
  double result_asc = 0;
  double mean = 0, err = 0;

  int j;

  if (n % 2 == 0)
    {
      result_gauss = f_center * wg[n / 2 - 1];
    }

  for (j = 0; j < (n - 1) / 2; j++)
    {
      const int jtw = j * 2 + 1;        /* j=1,2,3 jtw=2,4,6 */
      const double abscissa = half_length * xgk[jtw];
      const double fval1 = GSL_FN_EVAL (f, center - abscissa);
      const double fval2 = GSL_FN_EVAL (f, center + abscissa);
      const double fsum = fval1 + fval2;
      fv1[jtw] = fval1;
      fv2[jtw] = fval2;
      result_gauss += wg[j] * fsum;
      result_kronrod += wgk[jtw] * fsum;
      result_abs += wgk[jtw] * (fabs (fval1) + fabs (fval2));
    }

  for (j = 0; j < n / 2; j++)
    {
      int jtwm1 = j * 2;
      const double abscissa = half_length * xgk[jtwm1];
      const double fval1 = GSL_FN_EVAL (f, center - abscissa);
      const double fval2 = GSL_FN_EVAL (f, center + abscissa);
      fv1[jtwm1] = fval1;
      fv2[jtwm1] = fval2;
      result_kronrod += wgk[jtwm1] * (fval1 + fval2);
      result_abs += wgk[jtwm1] * (fabs (fval1) + fabs (fval2));
    };

  mean = result_kronrod * 0.5;

  result_asc = wgk[n - 1] * fabs (f_center - mean);

  for (j = 0; j < n - 1; j++)
    {
      result_asc += wgk[j] * (fabs (fv1[j] - mean) + fabs (fv2[j] - mean));
    }

  /* scale by the width of the integration region */

  err = (result_kronrod - result_gauss) * half_length;

  result_kronrod *= half_length;
  result_abs *= abs_half_length;
  result_asc *= abs_half_length;

  *result = result_kronrod;
  *resabs = result_abs;
  *resasc = result_asc;
  *abserr = rescale_error (err, result_abs, result_asc);

}

/* END_QK_H */

