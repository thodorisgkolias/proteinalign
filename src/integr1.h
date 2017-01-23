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
#include "func.h"
#include "error.h"
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

#include "initialise.h"
#include "set_initial.h"
#include "qpsrt.h"
#include "util.h"

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

