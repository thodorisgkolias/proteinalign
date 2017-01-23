/* EM_C.c
 * EM algortihm in C
 * Returns rotated data, mean, sigma, likelihood
 * Uses helmertized landmarks to remove location
 * Fixed issues with memory leakage
 * Possible addition of simple.rotation(M1,Mproc)
 * Last Update : 29/9/2016
 */







#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "svdfun.h"
#include "mycfun.h"
#include "matdet.h"
#include "shapecfun.h"
#include <limits.h>
#include <float.h>
#include "integr1.h"
#include "func_int.h"



// Create a 3-dimensional diagonal matrix
struct matrix *diag(double *w) {
     struct matrix *temp;
     temp = matrix_blank(3,3);
     temp -> e[0][0] = w[0];
     temp -> e[1][1] = w[1];
     temp -> e[2][2] = w[2];
     return temp;
     matrix_free(temp);
}

// myrotsimlaplace in C
struct vector *dointeg(double *l) {
     
     double f[5], sum, ss[3], L[4], v, G[4];
     sum = l[0]+l[1]+l[2];
     ss[0] = 2*l[0] - sum;
     ss[1] = 2*l[1] - sum;
     ss[2] = 2*l[2] - sum;
     L[0] = 0;
     L[1] = -ss[0] + sum;
     L[2] = -ss[1] + sum;
     L[3] = -ss[2] + sum;
     v = 2*pow(M_PI,2);
     // Integration part//
     struct f_params params;
     params.L1 = L[0];
     params.L2 = L[1];
     params.L3 = L[2];
     params.L4 = L[3];
     gsl_integration_workspace * w
          = gsl_integration_workspace_alloc (1000);
     gsl_function F11, F21, F12, F22, F13, F23, F14, F24;
     gsl_function P1, P2;
     P1.function = &part1; P2.function = &part2;
     F11.function = &f11; F21.function = &f21; F12.function = &f12;
     F22.function = &f22; F13.function = &f13; F23.function = &f23;
     F14.function = &f14; F24.function = &f24;
     F11.params = &params; F21.params = &params; F12.params = &params;
     F22.params = &params; F13.params = &params; F23.params = &params;
     F14.params = &params; F24.params = &params;
     P1.params = &params; P2.params = &params;
     double resultp1, resultp2;
     double result11, result21, result12, result22, result13;
     double result23, result14, result24, error;
     
     gsl_integration_qag (&F11, 0, 1, 0, 1e-7, 1000,6,
                          w, &result11, &error);
     gsl_integration_qag (&F21, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &result21, &error);
     gsl_integration_qag (&F12, 0, 1, 0, 1e-7, 1000,6,
                          w, &result12, &error);
     gsl_integration_qag (&F22, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &result22, &error);
     gsl_integration_qag (&F13, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &result13, &error);
     gsl_integration_qag (&F23, 0, 1, 0, 1e-7, 1000,6,
                          w, &result23, &error);
     gsl_integration_qag (&F14, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &result14, &error);
     gsl_integration_qag (&F24, 0, 1, 0, 1e-7, 1000,6,
                          w, &result24, &error);
     gsl_integration_qag (&P1, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &resultp1, &error);
     gsl_integration_qag (&P2, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &resultp2, &error);
     
     
     f[0] = 2*M_PI * (resultp1 - resultp2);
     
     
     f[1] = 2*M_PI * (result11 + 0.5*result21);
     f[2] = 2*M_PI * (result12 + 0.5*result22);
     f[3] = 2*M_PI * (-0.5*result13 - result23);
     f[4] = 2*M_PI * (-0.5*result14 - (result24));
     /*printf ("result11= % .18f\n", result11);
     printf ("result21= % .18f\n", result21);
     printf ("result12= % .18f\n", result12);
     printf ("result22= % .18f\n", result22);
     printf ("result13= % .18f\n", result13);
     printf ("result23= % .18f\n", result23);
     printf ("result14= % .18f\n", result14);
     
     */
     
     
     gsl_integration_workspace_free (w);
     // End integration// 
     
     
     G[0] = f[0]/v*2/M_PI; G[1] = -f[2]/v*2/M_PI; 
     G[2] = -f[3]/v*2/M_PI; G[3] = -f[4]/v*2/M_PI;
     
     struct vector *out;
     out = vector_blank(4);
     out -> e[0] = log(G[0]) + sum;
     out -> e[1] = 1-2*(G[2]+G[3])/G[0];
     out -> e[2] = 1-2*(G[1]+G[3])/G[0];
     out -> e[3] = 1-2*(G[2]+G[1])/G[0];
     return out;
     vector_free(out);
}



// Sort function for 3-dimensional vector
void sort3(double *x) {
     double tmp;
     
     if (x[0] < x[1]) {
          tmp = x[0];
          x[0] = x[1];
          x[1] = tmp;
     }
     if (x[1] < x[2]) {
          if (x[0] < x[2]) {
               tmp = x[2];
               x[2] = x[1];
               x[1] = x[0];
               x[0] = tmp;
          }
          else {
               tmp = x[1];
               x[1] = x[2];
               x[2] = tmp;
          }
     }
}

// rotationfit in C
double rotationfit2 (struct matrix *X, struct matrix *M, double sigma,double l0) {
     
     int i,j,k,m;
     struct matrix *rot, *svdmat, *Mt, *V, *temp;
     struct vector *D;
     D = vector_blank(3);
     V = matrix_blank(3,3);
     m = M->ncol;
     k = M->nrow;     
     Mt = matrix_transpose(M);
     svdmat = matrix_mult(X,Mt);
     matrix_free(Mt);
     for (i=0;i<3;i++) {
         for (j=0;j<3;j++){
               svdmat ->e[i][j] = svdmat ->e[i][j]/(pow(sigma,2));
          }
     }
          // Do SVD
     dsvdcmp(svdmat->e,3,3,D->e,V->e);
    
     // Pass D values to double
     float ch1,ch2;
     ch1=Det(svdmat);
     ch2=Det(V);  
     double l[3];
     for (i=0;i<3;i++) {
          l[i] = D->e[i];
     }
     vector_free(D);
     // Rank values of D
     int r[3];
     if (l[0] > l[1] & l[0] > l[2]) {
          //printf("IN1\n");
          r[0] = 0;
          if (l[1] > l[2]) {
               r[1] = 1;
               r[2] = 2;
          } else {
               r[1] = 2;
               r[2] = 1;
          }
     }
     if (l[1] > l[0] & l[1] > l[2]) {
          //printf("IN2\n");
          r[1] = 0;
          if (l[0] > l[2]) {
               r[0] = 1;
               r[2] = 2;
          } else {
               r[0] = 2;
               r[2] = 1;
          }
     }
     if (l[2] > l[0] & l[2] > l[1]) {
          //printf("IN3\n");
          r[2] = 0;
          if (l[1] > l[0]) {
               r[1] = 1;
               r[0] = 2;
          } else {
               r[1] = 2;
               r[0] = 1;
          }
     }
     
     // Sort values of D
     sort3(l);
     // Sign
     if( ch1*ch2<0) {
          l[0] = l[0];
          l[1] = l[1];
          l[2] = l[2]*(-1);
     }
     
     // Reorder U,V
     struct matrix *tempV, *tempU;
     tempV = matrix_blank(3,3);
     tempU = matrix_blank(3,3);
     for (i=0;i<3;i++) {
          for (j=0;j<3;j++) {
               tempV ->e[i][j] = V ->e[i][r[j]];
               tempU ->e[i][j] = svdmat ->e[i][r[j]];
          }
     }
     matrix_free(svdmat);
     matrix_free(V);
     
     // Change sign
     if (ch1<0) {
          struct matrix *Dtemp, *U1;
          double vtemp[3];
          vtemp[0] = 1;
          vtemp[1] = 1;
          vtemp[2] = -1;
          Dtemp = diag(vtemp);
          U1 = matrix_copy(tempU);
          matrix_free(tempU);
          tempU = matrix_mult(U1,Dtemp);
          matrix_free(Dtemp);
          matrix_free(U1);
     }
     if (ch2<0) {
          struct matrix *Dtemp, *V1;
          double vtemp[3];
          vtemp[0] = 1;
          vtemp[1] = 1;
          vtemp[2] = -1;
          V1 = matrix_copy(tempV);
          matrix_free(tempV);
          Dtemp = diag(vtemp);
          tempV = matrix_mult(V1,Dtemp);
          matrix_free(Dtemp);
          matrix_free(V1);
     }
     
     
     // Do integration
     struct vector *l1;
     l1 = dointeg(l);
     // Reorder values from integration
     double temp1[3];
     for (i=0;i<3;i++) {
          temp1[i] = l1->e[i+1]; 
     }
     l0 = l1 ->e[0];
     l0 = l0 - log(2/M_PI);
     vector_free(l1);
     // Create diagonal matrix //
     struct matrix *D1;
     D1 = diag(temp1);
     // Rotate X
     rot = matrix_mult(tempV,D1);
     matrix_free(D1);
     matrix_free(tempV);
     D1 = matrix_transpose(tempU);
     matrix_free(tempU);
     V = matrix_mult(rot,D1);
     matrix_free(D1);
     temp = matrix_mult(V,X);
     matrix_free(V);
     matrix_free(rot);
     // Replace data X with rotated
     matrix_copy_inplace(temp,X);
     matrix_free(temp);
     return l0;
     
}


struct mult_ret {
     struct matrix3d *new1;
     double sig1;
     double lik1;
};



// EM function in C
struct mult_ret EM(struct matrix3d *data) {
     
     // define structures //
     struct matrix *M0, *M1, *temp, *diff, *templ;
     struct matrix3d *rotdata, *data_orig;
     double tol1, tol2, d, SS,  sigma, sigma2, lik, lik_rot;
     int k, m, n, i, j;
     tol1 = 1e-05;
     tol2 = 1e-05;
     d = 100;
     k = data -> nrow;
     m = data -> ncol;
     n = data -> nblo;
     data_orig = matrix3d_blank(k,m,n);
     temp = matrix_blank(k,m);
     for (i=0;i<n;i++) {
          get_matrix_inplace(data,i,temp);
          update_entry_matrix(data_orig,temp,i);
          
     }
    matrix_free(temp);
    
     // Remove location //
     struct matrix3d *res_gpa;
     struct matrix *meangpa;
     // GPA mean
     res_gpa = GPA1(data,0,0,1e-05,1e-05);
     meangpa = mean3dmatrix(res_gpa);
     temp = preshape_transform(meangpa,0);
     M0 = matrix_transpose(temp);
     matrix_free(meangpa);
     matrix3d_free(res_gpa);
     matrix_free(temp);
    
     struct matrix3d *data1;
     data1 = matrix3d_blank(k-1,m,n);
     // Helmetized landmarks
     temp = matrix_blank(k,m);
     for (i=0;i<n;i++) {
          get_matrix_inplace(data_orig,i,temp);
          templ = preshape(temp,0);
          update_entry_matrix(data1,templ,i);
          matrix_free(templ);
        
     }
     matrix_free(temp);
    
     // Find SS //
     double temp3[n];
     temp = matrix_blank(k-1,m);
     for (i=0;i<n;i++) {
          get_matrix_inplace(data1,i,temp);
          temp3[i] = pow(norm(temp),2);
     }
     matrix_free(temp);     
     SS = 0;
     for (i=0;i<n;i++) {
          SS +=temp3[i];
     }     
     SS = SS/(data1->nblo);
     sigma = sigma2 = 1.0; 
     lik_rot = lik = 0;

     /////////// Begin EM //////////////
    
     int count;
     count = 0;
     rotdata = matrix3d_blank(k-1,m,n);

     struct matrix *temp2;
     while (d>tol1 & count<200) {
          count +=1;
          lik =0;
          for (i=0;i<n;i++) {
              temp2 = matrix_blank(k-1,m);
               // Choose data
               get_matrix_inplace(data1,i,temp2);
               temp = matrix_transpose(temp2);
               matrix_free(temp2);
               // Do Rotation
               lik_rot = rotationfit2(temp,M0,sigma,0);
               // Update Rotated data
               temp2 = matrix_transpose(temp);
               matrix_free(temp);
               update_entry_matrix(rotdata,temp2,i);
               // Find likelihood
               get_matrix_inplace(data1,i,temp2);
               lik = lik+ lik_rot - pow(norm(temp2),2)/(2*pow(sigma,2)) -
                    pow(norm(M0),2)/(2*pow(sigma,2)) - 
                   log(sigma)*k*m;
               matrix_free(temp2);
          }
   
    // New Mean
    M1 = mean3dmatrix(rotdata);
    temp = matrix_transpose(M0);
    // Find Difference
    diff = matrix_blank(M1->nrow,M1->ncol);
    for (i=0;i<M1->nrow;i++) {
        for (j=0;j<M1->ncol;j++) {
            diff ->e[i][j] = M1->e[i][j] - temp->e[i][j];
        }
    }
    matrix_free(temp);
    d = norm(diff);
    //printf ("Norm= % .18f\n", d);
    
    matrix_free(diff);
    matrix_free(M0);
   
    //Update Mean
    M0 = matrix_transpose(M1);
    matrix_free(M1);
    // Update Sigma
    sigma2 = (SS - pow(norm(M0),2))/(k*m);
    sigma = sqrt(sigma2);
    
    }
    
     struct mult_ret res;
     res.lik1 = lik;
     res.sig1 = sigma;
     res.new1 = rotdata;
     matrix3d_free(data_orig);
     matrix3d_free(data1);
     matrix_free(M0);
     
     return res;
     matrix3d_free(rotdata);
     
  
}
     

void EM_C(double *vectortemp, int *nrow, int *ncol, int *nblo,double *out,double *sig, double *lik) {
          struct matrix3d *temp;
          temp=matrix3d_blank((int) *nrow, (int) *ncol,(int) *nblo);
          int i,j,k,n,nr,nc,nb;
          nr=(int) *nrow;
          nc=(int) *ncol;
          nb=(int) *nblo;
          for (k=0;k<temp->nblo;k++) {
               for (j=0;j<temp->ncol;j++) {
                    for (i=0;i<temp->nrow;i++) {
                         temp->e[i][j][k]=vectortemp[(k)*nc*nr+(nr)*j+i];      
                    }
               }
          }   
          struct mult_ret res;
          res=EM(temp);
          temp = res.new1;
          
          *sig = res.sig1;
          *lik = res.lik1;
          //lik = lik11;
          
          for (k=0;k<temp->nblo;k++) {
               for (j=0;j<temp->ncol;j++) {
                    for (i=0;i<temp->nrow;i++) {
                         out[(k)*nc*nr+(nr)*j+i]=temp->e[i][j][k];
                    }
               }
          }  
          matrix3d_free(temp);
     }



// myrotsimlaplac2 in C
void dointeg2(double *l, double *out) {
     
     double f[5], sum, ss[3], L[4], v, G[4];
     sum = l[0]+l[1]+l[2];
     ss[0] = 2*l[0] - sum;
     ss[1] = 2*l[1] - sum;
     ss[2] = 2*l[2] - sum;
     L[0] = 0;
     L[1] = -ss[0] + sum;
     L[2] = -ss[1] + sum;
     L[3] = -ss[2] + sum;
     v = 2*pow(M_PI,2);
     // Integration part//
     struct f_params params;
     params.L1 = L[0];
     params.L2 = L[1];
     params.L3 = L[2];
     params.L4 = L[3];
     gsl_integration_workspace * w
          = gsl_integration_workspace_alloc (1000);
     gsl_function F11, F21, F12, F22, F13, F23, F14, F24;
     gsl_function P1, P2;
     P1.function = &part1; P2.function = &part2;
     F11.function = &f11; F21.function = &f21; F12.function = &f12;
     F22.function = &f22; F13.function = &f13; F23.function = &f23;
     F14.function = &f14; F24.function = &f24;
     F11.params = &params; F21.params = &params; F12.params = &params;
     F22.params = &params; F13.params = &params; F23.params = &params;
     F14.params = &params; F24.params = &params;
     P1.params = &params; P2.params = &params;
     double resultp1, resultp2;
     double result11, result21, result12, result22, result13;
     double result23, result14, result24, error;
     
     gsl_integration_qag (&F11, 0, 1, 0, 1e-7, 1000,6,
                          w, &result11, &error);
     gsl_integration_qag (&F21, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &result21, &error);
     gsl_integration_qag (&F12, 0, 1, 0, 1e-7, 1000,6,
                          w, &result12, &error);
     gsl_integration_qag (&F22, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &result22, &error);
     gsl_integration_qag (&F13, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &result13, &error);
     gsl_integration_qag (&F23, 0, 1, 0, 1e-7, 1000,6,
                          w, &result23, &error);
     gsl_integration_qag (&F14, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &result14, &error);
     gsl_integration_qag (&F24, 0, 1, 0, 1e-7, 1000,6,
                          w, &result24, &error);
     gsl_integration_qag (&P1, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &resultp1, &error);
     gsl_integration_qag (&P2, 0, 1/sqrt(2), 0, 1e-7, 1000,6,
                          w, &resultp2, &error);
     
     
     f[0] = 2*M_PI * (resultp1 - resultp2);
     
     
     f[1] = 2*M_PI * (result11 + 0.5*result21);
     f[2] = 2*M_PI * (result12 + 0.5*result22);
     f[3] = 2*M_PI * (-0.5*result13 - result23);
     f[4] = 2*M_PI * (-0.5*result14 - (result24));
     //printf ("result11= % .18f\n", result11);
     // printf ("result21= % .18f\n", result21);
     // printf ("result12= % .18f\n", result12);
     // printf ("result22= % .18f\n", result22);
     // printf ("result13= % .18f\n", result13);
     // printf ("result23= % .18f\n", result23);
     // printf ("result14= % .18f\n", result14);
      
      
     
     
     gsl_integration_workspace_free (w);
     // End integration//
     
     
     G[0] = f[0]/v*2/M_PI; G[1] = -f[2]/v*2/M_PI;
     G[2] = -f[3]/v*2/M_PI; G[3] = -f[4]/v*2/M_PI;
     
     
     out[0] = log(G[0]) + sum;
     out[1] = 1-2*(G[2]+G[3])/G[0];
     out[2] = 1-2*(G[1]+G[3])/G[0];
     out[3] = 1-2*(G[2]+G[1])/G[0];
     //struct vector *out;
     //out = vector_blank(4);
     //out -> e[0] = log(G[0]) + sum;
     //out -> e[1] = 1-2*(G[2]+G[3])/G[0];
     //out -> e[2] = 1-2*(G[1]+G[3])/G[0];
     //out -> e[3] = 1-2*(G[2]+G[1])/G[0];
     //return out;
     //vector_free(out);
}

     
     

/*
int main(void) {
 struct matrix3d *data;
 data = matrix3d_blank(10,3,2);
 data -> e[0][0][0] =   7.534;
 data -> e[1][0][0] =  -3.922;
 data -> e[2][0][0] = -0.651;
 data -> e[3][0][0] = -12.201;
 data -> e[4][0][0] = -7.173;
 data -> e[5][0][0] =   -8.673;
 data -> e[6][0][0] =  -4.889;
 data -> e[7][0][0] = -0.108;
 data -> e[8][0][0] = -9.669;
 data -> e[9][0][0] = -5.350;

 data -> e[0][1][0] = -3.776 ;
 data -> e[1][1][0] =  -3.881;
 data -> e[2][1][0] = -2.752;
 data -> e[3][1][0] =  3.975;
 data -> e[4][1][0] = -2.314;
 data -> e[5][1][0] = -6.154;
 data -> e[6][1][0] = -6.162 ;
 data -> e[7][1][0] = -1.143;
 data -> e[8][1][0] = -0.447;
 data -> e[9][1][0] = 2.653;

 data -> e[0][2][0] =  -4.109;
 data -> e[1][2][0] = 4.044;
 data -> e[2][2][0] =  2.466;
 data -> e[3][2][0] = -3.121;
 data -> e[4][2][0] =2.811;
 data -> e[5][2][0] =  0.203;
 data -> e[6][2][0] = 0.598;
 data -> e[7][2][0] =  7.430;
 data -> e[8][2][0] = 4.998;
 data -> e[9][2][0] = 5.692;

 data -> e[0][0][1] = 30.479;
 data -> e[1][0][1] = 26.772;
 data -> e[2][0][1] =28.605;
 data -> e[3][0][1] =  25.829 ;
 data -> e[4][0][1] =  26.235;
 data -> e[5][0][1] =21.419;
 data -> e[6][0][1] = 23.212;
 data -> e[7][0][1] = 31.720;
 data -> e[8][0][1] = 26.850;
 data -> e[9][0][1] = 31.677;

 data -> e[0][1][1] = 35.369;
 data -> e[1][1][1] =  33.436;
 data -> e[2][1][1] = 33.965;
 data -> e[3][1][1] = 19.825;
 data -> e[4][1][1] = 30.058;
 data -> e[5][1][1] = 30.253;
 data -> e[6][1][1] = 32.762;
 data -> e[7][1][1] = 36.289;
 data -> e[8][1][1] = 29.021;
 data -> e[9][1][1] = 30.275;

 data -> e[0][2][1] = 22.374;
 data -> e[1][2][1] =  9.197;
 data -> e[2][2][1] = 12.503;
 data -> e[3][2][1] =  8.494;
 data -> e[4][2][1] =  7.497;
 data -> e[5][2][1] =  9.620;
 data -> e[6][2][1] = 11.891;
 data -> e[7][2][1] =  9.176;
 data -> e[8][2][1] = 3.898;
 data -> e[9][2][1] = 6.639;

 EM(data);
 
}
*/

/*
int main(void) {
     double data[60] = {7.534,-3.922,-0.651,-12.201,-7.173,-8.673,-4.889,-0.108,-9.669,-5.350,-3.776,-3.881,-2.752,3.975,-2.314,-6.154,-6.162,-1.143,-0.447,2.653,-4.109,4.044,2.466,-3.121,2.811,0.203,0.598,7.430,4.998,5.692,30.479,26.772,28.605,25.829,26.235,21.419,23.212,31.720,26.850,31.677,35.369,33.436,33.965,19.825,30.058,30.253,32.762,36.289,29.021,30.275,22.374,9.197,12.503,8.494,7.497,9.620,11.891,9.176,3.898,6.639};
     double *vec;
     vec = data;
     int nr[1],nc[1],nb[1];
     nr[0] = 10;
     nc[0] = 3;
     nb[0] = 2;
     int *nrow, *ncol, *nblo;
     nrow = nr;
     ncol = nc;
     nblo = nb;
     double lik1[1], sig1[1];
     double *out, *lik, *sig;
     lik1[0] = 0;
     sig1[0] = 1;
     lik = lik1;
     sig = sig1;
     out = data;
     
     EM_C(vec,nrow,ncol,nblo,out,sig,lik);     

     }
*/

