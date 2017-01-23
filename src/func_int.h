 /* func_int.h
 * 
 * Functions and structure for integration 
 * with Deriv.lap.invC.c
 * Last update : 27/8/2016
 */

 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <float.h>
//#include "gsl_functions.h"

struct f_params {
     double L1;
     double L2;
     double L3;
     double L4;
};

// Int lims (0,1) //

double f11(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l2 - l1;
     double b1 = l3 - l2;
     double b2 = l4 - l2;
     double bb1 = b1 + (1 - pow(x,2)) * u;
     double bb2 = b2 + (1 - pow(x,2)) * u;
     double temp = -l2 + (1 - pow(x,2)) * u;
     double temp2 = sqrt(bb1 * bb2);
     return 2*( (exp(temp) / (bb1 * temp2) +
                exp(temp) / (bb2 * temp2)) / 
                2 - (exp(temp) / temp2)) *
                sqrt(1 - pow(x,2));
}


// Int lims (0,1/sqrt(2)) //

double f21(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l4 - l3;
     double b1 = l1 - l4;
     double b2 = l2 - l4;
     double bb1 = b1 + pow(x,2) * u;
     double bb2 = b2 +  pow(x,2) * u;
     double bb3 = b1 + (1 - pow(x,2))*u;
     double bb4 = b2 + (1 - pow(x,2))*u;
     double temp = -l4 + pow(x,2) * u;
     double temp2 = sqrt(bb1 * bb2);
     double temp3 = -l4 + (1-pow(x,2)) * u; 
     double temp4 = sqrt(bb3 * bb4);
     return 2*( (exp(temp) / (bb1 * temp2)) +
                (exp(temp3) / (bb3 * temp4))) / sqrt(1 - pow(x,2));
}     

// Int lims (0,1) //

double f12(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l2 - l1;
     double b1 = l3 - l2;
     double b2 = l4 - l2;
     double bb1 = b1 +  pow(x,2) * u;
     double bb2 = b2 +  pow(x,2) * u;
     double temp = -l2 +  pow(x,2) * u;
     double temp2 = sqrt(bb1 * bb2);
     return 2*( (exp(temp) / (bb1 * temp2) +
                exp(temp) / (bb2 * temp2)) / 
                2 - (exp(temp) / temp2)) *
                sqrt(1 - pow(x,2));
}

// Int lims (0,1/sqrt(2)) //
double f22(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l4 - l3;
     double b1 = l1 - l4;
     double b2 = l2 - l4;
     double bb1 = b1 + pow(x,2) * u;
     double bb2 = b2 +  pow(x,2) * u;
     double bb3 = b1 + (1 - pow(x,2))*u;
     double bb4 = b2 + (1 - pow(x,2))*u;
     double temp = -l4 + pow(x,2) * u;
     double temp2 = sqrt(bb1 * bb2);
     double temp3 = -l4 + (1-pow(x,2)) * u; 
     double temp4 = sqrt(bb3 * bb4);
     return 2*( (exp(temp) / (bb2 * temp2)) +
                (exp(temp3) / (bb4 * temp4))) / sqrt(1 - pow(x,2));
}     


// Int lims (0,1/sqrt(2)) //

double f13(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l2 - l1;
     double b1 = l3 - l2;
     double b2 = l4 - l2;
     double bb1 = b1 +  pow(x,2) * u;
     double bb2 = b2 +  pow(x,2) * u;
     double bb3 = b1 + (1 - pow(x,2)) * u;
     double bb4 = b2 + (1 - pow(x,2)) * u;
     double temp = -l2 +  pow(x,2) * u;
     double temp2 = sqrt(bb1 * bb2);
     double temp3 = -l2 + (1 - pow(x,2)) * u;
     double temp4 = sqrt(bb3 * bb4);
     
     return 2* ( (exp(temp) / (bb1 * temp2)) +
           (exp(temp3) / (bb3 * temp4))
     ) / sqrt(1 - pow(x,2));
     
}


// Int lims (0,1) //

double f23(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l4 - l3;
     double b1 = l1 - l4;
     double b2 = l2 - l4;
     double bb1 = b1 + (1 - pow(x,2) * u);
     double bb2 = b2 + (1 - pow(x,2) * u);
     double temp = -l4 + (1 - pow(x,2) * u);
     double temp2 = sqrt(bb1 * bb2);
     
     return sqrt(1 - pow(x,2)) * (
               (exp(temp) / (bb1 * temp2)) +
                    (exp(temp) / (bb2 * temp2))  - 
                    2 * (exp(temp)/ temp2)
     );
}

// Int lims (0,1/sqrt(2)) //

double f14(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l2 - l1;
     double b1 = l3 -l2;
     double b2 = l4 -l2;
     double bb1 = b1 + pow(x,2) * u;
     double bb2 = b2 + pow(x,2) * u;
     double bb3 = b1 + (1 - pow(x,2)) * u;
     double bb4 = b2 + (1 - pow(x,2)) * u;
     double temp = -l2 + pow(x,2) * u;
     double temp2 = sqrt(bb1 * bb2);
     double temp3 = -l2 + (1 - pow(x,2)) * u;
     double temp4 = sqrt(bb3 * bb4);
     
     return 2 * (
          (exp(temp) / (bb2 * temp2)) +
               (exp(temp3) / (bb4 * temp4))
     ) / sqrt(1 - pow(x,2));
}


// Int lims (0,1) //

double f24(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l4 - l3;
     double b1 = l1 - l4;
     double b2 = l2 - l4;
     double bb1 = b1 + pow(x,2) * u;
     double bb2 = b2 + pow(x,2) * u;
     double temp = -l4 + pow(x,2) * u;
     double temp2 = sqrt(bb1 * bb2);
     
     return 2 * (
          ((exp(temp) / (bb1*temp2)) + 
               (exp(temp) / (bb2 * temp2))
          ) / 2 - 
               (exp(temp)/ temp2)
     ) * sqrt(1 - pow(x,2));
     
}     

// Int lims (0,1/sqrt(2))

double part1(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l2-l1;
     double b1 = l3 - l2;
     double b2 = l4 - l2;
     double bb1 = b1 + pow(x,2) * u;
     double bb2 = b2 + pow(x,2) * u;
     double bb3 = b1 + (1 - pow(x,2)) * u;
     double bb4 = b2 + (1 - pow(x,2)) * u;
     double temp = -l2 + pow(x,2) * u;
     double temp2 = sqrt(bb1 * bb2);
     double temp3 = -l2 + (1 - pow(x,2)) * u;
     double temp4 = sqrt(bb3 * bb4);
     
     return 2 * (
          (exp(temp) / temp2) + 
               (exp(temp3) / temp4)
     ) / sqrt(1 - pow(x,2));
}



// Int lims (0,1/sqrt(2))

double part2(double x, void *p) {
     struct f_params * params
     = (struct f_params *)p;
     double l1 = (params->L1);
     double l2 = (params->L2);
     double l3 = (params->L3);
     double l4 = (params->L4);
     double u = l4-l3;
     double b1 = l1 - l4;
     double b2 = l2 - l4;
     double bb1 = b1 + pow(x,2) * u;
     double bb2 = b2 + pow(x,2) * u;
     double bb3 = b1 + (1 - pow(x,2)) * u;
     double bb4 = b2 + (1 - pow(x,2)) * u;
     double temp = -l4 + pow(x,2) * u;
     double temp2 = sqrt(bb1 * bb2);
     double temp3 = -l4 + (1 - pow(x,2)) * u;
     double temp4 = sqrt(bb3 * bb4);
     
     return 2 * (
               (exp(temp) / temp2) + 
                    (exp(temp3) / temp4)
     ) / sqrt(1 - pow(x,2));
}




