/*
  These is the library of general functions.
 Commends updated 10/04/2004 Alfred Kume
 */

struct vector
{
   int nrow;
   double *e;
};

struct matrix
{
   int nrow;
   int ncol;
   double **e;
};

struct matrix3d
{
   int nrow;
   int ncol;
   int nblo;
   double ***e;
};





struct vector *vector_alloc(FILE *fp,int len)
/* Function reads teh vector of length len from teh data file with pointer fp */
{
   int i;
   struct vector *temp;
   temp = NULL;
   temp = (struct vector *)calloc(1,sizeof(struct vector));
   if (temp == NULL)
   {
           //printf("Not enough memory to continue");
          return(NULL);
   }
/*
 * Read length of vector
 */
 temp->nrow=len;
   /*if (fscanf(fp,"%d",&(temp->nrow))!=1)
   {
      printf("File not of correct format");
          return(NULL);
   }
   */
   temp->e = (double *)calloc(temp->nrow,sizeof(double));
   if (temp->e == NULL )
   {
      //printf("Not enough memory to store vector");
          return(NULL);
   }
   for (i=0;i<temp->nrow;i++)
   {
/*
 * Read elements of vector
 */
      if (fscanf(fp,"%lf",&(temp->e[i])) != 1)
          {
        // printf("File not of correct format");
                 return(NULL);
          }
   }
   return temp;
}

struct matrix *matrix_alloc(FILE *fp, int nrow, int ncol)
/* Function reads teh matrix nrow X nclo from teh data file with pointer fp */

{
   int i,j;
   struct matrix *temp;
   temp = NULL;
   temp = (struct matrix *)calloc(1,sizeof(struct matrix));
   if (temp == NULL)
   {
    //  printf("Not enough memory to store matrix");
          return(NULL);
   }
/* 
 * Read dimensions of matrix
 */
   temp->nrow=nrow;
   temp->ncol=ncol;
   temp->e = (double **)calloc(temp->nrow,sizeof(double *));
   if (temp->e == NULL )
   {
     // printf("Not enough memory to store matrix");
          return(NULL);
   }
   for (i=0;i<temp->nrow;i++)
   {
      temp->e[i] = (double *)calloc(temp->ncol,sizeof(double));
      if (temp->e[i]  == NULL )
          {
        // printf("Not enough memory to store matrix");
                 return(NULL);
          }
   }
   for (i=0;i<temp->nrow;i++)
   {
      for (j=0;j<temp->ncol;j++)
      {
/*
 * Read elements of matrix along the columns
 */
         if (fscanf(fp,"%lf",&(temp->e[i][j])) != 1)
                 {
           // printf("File not of correct format");
                        return(NULL);
                 }
      }
   }
   return temp;
}


struct matrix3d *matrix3d_alloc(FILE *fp, int nrow, int ncol,int nblo)
/* Function reads a 3d matrix nrow X nclo X nblo from teh data file with pointer fp */

{
   int i,j,k;
   struct matrix3d *temp;
   temp = NULL;
   temp = (struct matrix3d *)calloc(1,sizeof(struct matrix3d));
   if (temp == NULL)
   {
    //  printf("Not enough memory to store matrix");
          return(NULL);
   }
/* 
 * Read dimensions of matrix
 */
   temp->nrow=nrow;
   temp->ncol=ncol;
   temp->nblo=nblo;
   temp->e = (double ***)calloc(temp->nrow,sizeof(double **));
   if (temp->e == NULL )
   {
    //  printf("Not enough memory to store matrix");
          return(NULL);
   }
   for (i=0;i<temp->nrow;i++)
   {
      
      temp->e[i] = (double **)calloc(temp->ncol,sizeof(double *));
      if (temp->e[i]  == NULL )
          {
      //   printf("Not enough memory to store matrix");
                 return(NULL);
          }
          
          for (j=0;j<temp->ncol;j++)
          {
      temp->e[i][j] = (double *)calloc(temp->nblo,sizeof(double ));
      if (temp->e[i][j]  == NULL )
          {
     //    printf("Not enough memory to store matrix");
                 return(NULL);
          }
          }
   }
   
   for (k=0;k<temp->nblo;k++)
{
   
    for (j=0;j<temp->ncol;j++)
    {
     for (i=0;i<temp->nrow;i++)
        {
/*
 * Read elements of matrix along the columns
 */
         if (fscanf(fp,"%lf",&(temp->e[i][j][k])) != 1)
                 {
       //     printf("File not of correct format");
                        return(NULL);
                 }
      }
   }
}   
   
   return temp;
}


/*
 * Function sets up a 2dmatrix struct from a block of a 3dmatrix
 */

struct matrix *get_matrix(struct matrix3d *mat3d,int nbl)
{
   int i,j;
   struct matrix *temp;
   temp =  NULL;
   temp = (struct matrix *)calloc(1,sizeof(struct matrix));
   if (temp == NULL)
   {
       //    printf("MCMC Error 0026 : Not enough memory to continue");
          return(NULL);
   }
   temp->nrow = mat3d->nrow;
   temp->ncol = mat3d->ncol;
   temp->e = (double **)calloc(temp->nrow,sizeof(double *));
   if (temp->e == NULL )
   {
        //   printf("MCMC Error 0027 : Not enough memory to continue");
          return(NULL);
   }
   for (i=0;i<temp->nrow;i++)
   {
      temp->e[i] = (double *)calloc(temp->ncol,sizeof(double));
      if (temp->e[i]  == NULL )
          {
            //      printf("MCMC Error 0028 : Not enough memory to continue");
                 return(NULL);
          }
   }
   for (i=0;i<temp->nrow;i++)
   {
      for (j=0;j<temp->ncol;j++)
      {
/*
 * Fill in entries 
 */
         temp->e[i][j] = mat3d->e[i][j][nbl];
      }
   }
   return temp;
}

/*

*/

void get_matrix_inplace(struct matrix3d *mat3d,int nbl, struct matrix *temp)
{
 int i,j; 
   for (i=0;i<temp->nrow;i++)
   {
      for (j=0;j<temp->ncol;j++)
      {
/*
 * Fill in entries 
 */
         temp->e[i][j] = mat3d->e[i][j][nbl];
      }
   }
}


void update_entry_matrix(struct matrix3d *mat3d,struct matrix *mat2d,int nbl)
{
   int i,j;
      for (i=0;i<mat3d->nrow;i++)
   {
      for (j=0;j<mat3d->ncol;j++)
      {
/*
 * Fill in zeroes in matrix 
 */
         mat3d->e[i][j][nbl] = mat2d->e[i][j];
      }
   }

}

/*
 * Function frees memory in a 3dmatrix
 */
void matrix3d_free(struct matrix3d *mat)
{
  
   int i,j;
   if (mat != NULL)
   {
          if(mat->nrow > 0)
          {
          if (mat->e != NULL)
          {
                  
                  for (i=0;i<mat->nrow;i++)
                  {
                  for (j=0;j<mat->ncol;j++)
                  {
                  
                        if(mat->e[i][j] != NULL)
                                free(mat->e[i][j]);
                  }
                  
                        if(mat->e[i] != NULL)
                                free(mat->e[i]);
                  }
                  free(mat->e);
          }
      free(mat);
          mat = NULL;
          }
   }

}








/*
 * Function frees memory in a matrix
 */
void matrix_free(struct matrix *mat)
{
  
   int i;
   if (mat != NULL)
   {
          if(mat->nrow > 0)
          {
          if (mat->e != NULL)
          {
                  
                  for (i=0;i<mat->nrow;i++)
                  {
                        if(mat->e[i] != NULL)
                                free(mat->e[i]);
                  }
                  free(mat->e);
          }
      free(mat);
          mat = NULL;
          }
   }

}

/*
 * Function frees a list of matrices
 */
void matrix_list_free(struct matrix **matlist, int len)
{
   int i;
   if(matlist != NULL)
   {
           for (i = 0; i < len; i++)
           {
             if(matlist[i] != NULL)
            matrix_free(matlist[i]);
                 else 
                         break;
           }
        free(matlist);
                matlist = NULL;
   }
}


/*
 * Function frees a vector
 */
void vector_free(struct vector *vec)
{
   if(vec != NULL)
   {
                if(vec->e != NULL)
                        free(vec->e);
                free(vec);
                vec = NULL;
   }
}

/*
 * Function frees a list of vectors
 */ 
void vector_list_free(struct vector **veclist, int len)
{
   int i;
   if (veclist != NULL)
   {
      for (i = 0; i < len; i++)
      {
             if (veclist[i] != NULL)
             {
                        vector_free(veclist[i]);
             }
      }
      free(veclist);
          veclist = NULL;
  }
}




/*
 * Function creates a vector of zeros of length nrow 
 */
struct vector *vector_blank(int nrow)
{
   int i;
   struct vector *temp;
   temp =  NULL;
   temp = (struct vector *)calloc(1,sizeof(struct vector));
   if (temp == NULL)
   {
        //   printf("MCMC Error 0024 : Not enough memory to continue");
          return(NULL);
   }
   temp->nrow = nrow;
   temp->e = (double *)calloc(temp->nrow,sizeof(double));
   if (temp->e == NULL )
   {
        //   printf("MCMC Error 0025 : Not enough memory to continue");
          return(NULL);
   }
   for (i=0;i<temp->nrow;i++)
   {
/*
 * Fill in zeroes in matrix 
 */
         temp->e[i] = 0;
   }
   return temp;
}

/*
 * Function copys vector original and returns the copy
 */      
struct vector *vector_copy(struct vector *original)
{
   int i;
   struct vector *temp;
   temp = vector_blank(original->nrow);
   if(temp == NULL)
          return(NULL);
   for (i=0;i<temp->nrow;i++)
   {
/*
 * Copy elements of vector
 */
         temp->e[i] = original->e[i];
   }
   return temp;
}
/*
 * Function creates a matrix of zeros of dimensions nrow * ncol
 */
struct matrix *matrix_blank(int nrow,int ncol)
{
   int i,j;
   struct matrix *temp;
   temp =  NULL;
   temp = (struct matrix *)calloc(1,sizeof(struct matrix));
   if (temp == NULL)
   {
        //   printf("MCMC Error 0026 : Not enough memory to continue");
          return(NULL);
   }
   temp->nrow = nrow;
   temp->ncol = ncol;
   temp->e = (double **)calloc(temp->nrow,sizeof(double *));
   if (temp->e == NULL )
   {
        //   printf("MCMC Error 0027 : Not enough memory to continue");
          return(NULL);
   }
   for (i=0;i<temp->nrow;i++)
   {
      temp->e[i] = (double *)calloc(temp->ncol,sizeof(double));
      if (temp->e[i]  == NULL )
          {
            //      printf("MCMC Error 0028 : Not enough memory to continue");
                 return(NULL);
          }
   }
   for (i=0;i<temp->nrow;i++)
   {
      for (j=0;j<temp->ncol;j++)
      {
/*
 * Fill in zeroes in matrix 
 */
         temp->e[i][j] = 0.0;
      }
   }
   return temp;
}

/*
Function creates a 3d matrix nrow X nclo X nblo with zero entries 
 */
struct matrix3d *matrix3d_blank(int nrow, int ncol,int nblo)


{
   int i,j,k;
   struct matrix3d *temp;
   temp = NULL;
   temp = (struct matrix3d *)calloc(1,sizeof(struct matrix3d));
   if (temp == NULL)
   {
   //   printf("Not enough memory to store matrix");
          return(NULL);
   }
/* 
 * Read dimensions of matrix
 */
   temp->nrow=nrow;
   temp->ncol=ncol;
   temp->nblo=nblo;
   temp->e = (double ***)calloc(temp->nrow,sizeof(double **));
   if (temp->e == NULL )
   {
     // printf("Not enough memory to store matrix");
          return(NULL);
   }
   for (i=0;i<temp->nrow;i++)
   {
      
      temp->e[i] = (double **)calloc(temp->ncol,sizeof(double *));
      if (temp->e[i]  == NULL )
          {
      //   printf("Not enough memory to store matrix");
                 return(NULL);
          }
          
          for (j=0;j<temp->ncol;j++)
          {
      temp->e[i][j] = (double *)calloc(temp->nblo,sizeof(double ));
      if (temp->e[i][j]  == NULL )
          {
      //   printf("Not enough memory to store matrix");
                 return(NULL);
          }
          }
   }
   
   for (k=0;k<temp->nblo;k++)
{
   for (i=0;i<temp->nrow;i++)
   {
      for (j=0;j<temp->ncol;j++)
      {
/*
 * Input zeros to the  elements of matrix along the columns
 */
       temp->e[i][j][k]=0;
      }
   }
}   
   
   return temp;
}





/*
 * Function copys matrix original and returns the copy
 */      
struct matrix *matrix_copy(struct matrix *original)
{
   int i,j;
   struct matrix *temp;
   temp = matrix_blank(original->nrow,original->ncol);
   if(temp == NULL)
          return(NULL);
   for (i=0;i<temp->nrow;i++)
   {
      for (j=0;j<temp->ncol;j++)
      {
/*
 * Copy elements of matrix
 */
         temp->e[i][j] = original->e[i][j];
      }
   }
   return temp;
}

/*
 * Function copys matrix original and returns the copy
 */      
void matrix_copy_inplace(struct matrix *original, struct matrix *copy)
{
   int i,j;
   for (i=0;i<copy->nrow;i++)
   {
      for (j=0;j<copy->ncol;j++)
      {
/*
 * Copy elements of matrix
 */
         copy->e[i][j] = original->e[i][j];
      }
   }
}

/*
 * Function creates an n by n identity matrix
 */
struct matrix *matrix_identity(int n)
{
   int i;
   struct matrix *temp;
   temp = matrix_blank(n,n);
   if(temp == NULL)
          return(NULL);
   for (i=0;i<temp->nrow;i++)
   {
      temp->e[i][i] = 1;
   }
   return temp;
}


/*
 * Function adds matrices mat1 and mat2 and returns the result
 */
struct matrix *matrix_add(struct matrix *mat1, struct matrix *mat2)
{
   int i,j;
   struct matrix *temp;
   if (mat1->ncol != mat2->ncol)
   {
        //   printf("MCMC Error 0029 : Wrong parameters for matrix_add routine");
          return(NULL);
   }
   if (mat1->nrow != mat2->nrow)
   {
        //   printf("MCMC Error 0030 : Wrong parameters for matrix_add routine");
          return(NULL);
   }
   temp = matrix_copy(mat1);
   if(temp == NULL)
      return(NULL);
   for (i=0;i<temp->nrow;i++)
   {
      for (j=0;j<temp->ncol;j++)
      {
/*
 * Add element from matrix mat2 to element from matrix mat1 
 */
         temp->e[i][j] += mat2->e[i][j];
      }
   }
   return temp;
}

/*
 * Function subtracts matrices mat2 from mat1 and returns the result
 * NOT USED IN OUR CODE!!!
 */
struct matrix *matrix_sub(struct matrix *mat1, struct matrix *mat2)
{
   int i,j;
   struct matrix *temp;
   if (mat1->ncol != mat2->ncol)
   {
        //   printf("MCMC Error 0031 : Wrong parameters for matrix_sub routine");
          return(NULL);
   }
   if (mat1->nrow != mat2->nrow)
   {
        //   printf("MCMC Error 0032 : Wrong parameters for matrix_sub routine");
          return(NULL);
   }
   temp = matrix_copy(mat1);
   if(temp == NULL)
      return(NULL);
   for (i=0;i<temp->nrow;i++)
   {
      for (j=0;j<temp->ncol;j++)
      {
/*
 * Subtract element from matrix mat2 to element from matrix mat1 
 */
         temp->e[i][j] -= mat2->e[i][j];
      }
   }
   return temp;
}


/*
 * Function subtracts matrix mat2 from mat1 and updates mat1 with the result!!!
 */
void matrix_sub_inplace(struct matrix *mat1, struct matrix *mat2)
{
   int i,j;
   if (mat1->ncol != mat2->ncol)
   {
       //    printf("MCMC Error 0031 : Wrong parameters for matrix_sub routine");
   }
   if (mat1->nrow != mat2->nrow)
   {
        //   printf("MCMC Error 0032 : Wrong parameters for matrix_sub routine");
   }

   for (i=0;i<mat1->nrow;i++)
   {
      for (j=0;j<mat1->ncol;j++)
      {
/*
 * Subtract element from matrix mat2 to element from matrix mat1 
 */
         mat1->e[i][j] -= mat2->e[i][j];
      }
   }

}




/*
 * Function multiplies matrices mat1 and mat2 and returns the result
 */
struct matrix *matrix_mult(struct matrix *mat1, struct matrix *mat2)
{
   int i,j,terms;
   struct matrix *temp;
   if (mat1->ncol != mat2->nrow)
   {
       //    printf("MCMC Error 0033 : Wrong parameters for matrix_mult routine");
          return(NULL);
   }
   temp = matrix_blank(mat1->nrow,mat2->ncol);
   if(temp == NULL)
          return(NULL);
   for(i=0;i<temp->nrow;i++)
   {
      for(j=0;j<temp->ncol;j++)
      {
         for(terms=0;terms<mat1->ncol;terms++)
         {
            temp->e[i][j] += mat1->e[i][terms] * mat2->e[terms][j];
         }
      }
   }
   return temp;
}

/*
 * Function multiplies matrix mat by scalar s and returns the result
 * NOT USED IN MLN!!!
 */
struct matrix *ms_mult(struct matrix *mat, double s)
{
   int i,j;
   struct matrix *temp;
   temp = matrix_copy(mat);
   if(temp == NULL)
      return(NULL);
   for(i=0;i<temp->nrow;i++)
      for(j=0;j<temp->ncol;j++)
         temp->e[i][j] = mat->e[i][j] * s;
   return temp;
}

/*
 * Function finds the transpose of a matrix
 */
struct matrix *matrix_transpose(struct matrix *original)
{
   int i,j;
   struct matrix *temp;
   temp = matrix_blank(original->ncol,original->nrow);
   if(temp == NULL)
      return(NULL);
   for(i=0;i<temp->nrow;i++)
   {
      for(j=0;j<temp->ncol;j++)
      {
         temp->e[i][j] += original->e[j][i];
      }
   }
   return temp;
}
/*
 * Function finds the (first 2 dimensional) mean matrix of a 3d matrix
 */

struct matrix *mean3dmatrix(struct matrix3d *mat3d)
{
int i,j,k;
struct matrix *temp;

/*      
temp=get_matrix(mat3d,0);     
for (i=1;i<mat3d->nblo;i++)
temp=matrix_add(temp,get_matrix(mat3d,i));
temp=ms_mult(temp, 1/((double) mat3d->nblo));
*/
temp=matrix_blank(mat3d->nrow,mat3d->ncol);
for (i=0;i<mat3d->ncol;i++)
    for (j=0;j<mat3d->nrow;j++)
       for (k=0;k<mat3d->nblo;k++)
           temp->e[j][i]+=mat3d->e[j][i][k];

           
for (i=0;i<mat3d->ncol;i++)
    for (j=0;j<mat3d->nrow;j++)      
           temp->e[j][i]/=(double)mat3d->nblo;
           

return temp;
}

/* Function finds the (first 2 dimensional) sum matrix of a 3d matrix and its entries are 
 stored at mat
 */

void sum3dmatrix_inplace(struct matrix3d *mat3d, struct matrix *mat)
{
int i,j,k;
/*get_matrix_inplace(mat3d,0,mat);*/
get_matrix_inplace(mat3d,0,mat);
for (i=1;i<mat3d->nblo;i++)
   for (j=0;j<mat3d->ncol;j++)
       for (k=0;k<mat3d->nrow;k++)
            mat->e[k][j]+=mat3d->e[k][j][i];
/*
for (j=0;j<mat3d->ncol;j++)
       for (k=0;k<mat3d->nrow;k++)
          mat->e[k][j]/=((double)mat3d->nblo);
*/              
}




/*
 Calculates the covariance matrix of the column vectors of mat 
 if ind is 1 then it calculates the correlation matrix 
 */
struct matrix *covar(struct matrix *A, int k)                        
{
struct matrix *B,*temp;
int i,j;
double s;
temp=matrix_blank(A->ncol,A->ncol);
B=matrix_blank(A->nrow,A->ncol);
matrix_copy_inplace(A,B);
/*B=A;*/
for (i=0;i<B->ncol;i++)
  {
s=0;
for (j=0;j<B->nrow;j++)
     {
     s+=B->e[j][i];
     }
s=s/((float)B->nrow);
   
   for (j=0;j<B->nrow;j++)
    B->e[j][i]-=s;
   
 }
if (k==1)
         {
         for (i=0;i<B->ncol;i++)
         {
s=0;
for (j=0;j<B->nrow;j++)
                {
     s+=B->e[j][i]*B->e[j][i];
                }
s=sqrt(s);  
      for (j=0;j<B->nrow;j++)
                        {
                       /* printf("%lf \n", B->e[j][i]);*/
    
                                    B->e[j][i]=((s) == (0) ?  (0) :(B->e[j][i]/s));   
                        }          
  
         }
        }                  


matrix_free(temp);
/*temp=matrix_blank(A->ncol,A->nrow); */  

temp=matrix_blank(B->ncol,B->ncol);
for (i=0;i<temp->ncol;i++)
  {for (k=0;k<i+1;k++)
   {
   temp->e[i][k]=0;
        for (j=0;j<B->nrow;j++)         
            {
            temp->e[i][k]+=B->e[j][i]*B->e[j][k];            
            }           
  temp->e[k][i]=temp->e[i][k];
   }
  }
matrix_free(B);
/*     
tmp1=matrix_transpose(B);        
temp=matrix_mult(tmp1,B);
matrix_free(B);
matrix_free(tmp1);
*/
return(temp);
}     

/*
 Calculates the covariance matrix of the column vectors of mat 
 if ind is 1 then it calculates the correlation matrix 
 the result is input into Mat whose elements are 0.
  */
void covar_inplace(struct matrix *A, int k,struct matrix *Mat)                        
{
struct matrix *B;
int i,j,w;
double s;
B=matrix_blank(A->nrow,A->ncol);
matrix_copy_inplace(A,B);
/*B=A;*/
for (i=0;i<B->ncol;i++)
  {
s=0;
for (j=0;j<B->nrow;j++)
     {
     s+=B->e[j][i];
     }
s=s/((float)B->nrow);
   
   for (j=0;j<B->nrow;j++)
    B->e[j][i]-=s;
   
 }
if (k==1)
         {
         for (i=0;i<B->ncol;i++)
         {
s=0;
for (j=0;j<B->nrow;j++)
                {
     s+=B->e[j][i]*B->e[j][i];
                }
s=sqrt(s);  
      for (j=0;j<B->nrow;j++)
                        {
                       /* printf("%lf \n", B->e[j][i]);*/
    
                                    B->e[j][i]=((s) == (0) ?  (0) :(B->e[j][i]/s));   
                        }          
  
         }
        }                  

/*temp=matrix_blank(A->ncol,A->nrow); */       
for (i=0;i<Mat->ncol;i++)
  {
  for (w=0;w<Mat->ncol;w++)
   {
   Mat->e[i][w]=0;
        for (j=0;j<B->nrow;j++)         
            Mat->e[i][w]+= B->e[j][i]*B->e[j][w];
                        
  
   }
  }
matrix_free(B);
}     

/*
finds the index of the vector's minimum 
*/
double vector_min(struct vector *w)
{
double s;
int j;
s=w->e[0];
for (j=1;j<w->nrow;j++)
     ((s>w->e[j])?(s=w->e[j]):(s=s));
/*
    {  
 if (s<w->e[j])     
 {
  s=w->e[j];
  }
 } 
 */
return(s);
}
/*
finds the index of the vector's minimum 
*/
int vector_min_index(struct vector *w)
{
int s,j;
s=0;
for (j=1;j<w->nrow;j++)
     ((w->e[s]>w->e[j])?(s=j):(s=s));
return(s);
}

/*
finds the index of the vector's meaximum 
*/
int vector_max_index(struct vector *w)
{
int s,j;
s=0;
for (j=1;j<w->nrow;j++)
     ((w->e[s]<w->e[j])?(s=j):(s=s));
return(s);
}


/*
 * Function multiplies matrices transpose(mat2) and mat1 and returns the result
 */
struct matrix *matrix_tranp_mult(struct matrix *mat1, struct matrix *mat2)
{
   int i,j,terms;
   struct matrix *temp;
   if (mat1->ncol != mat2->ncol)
   {
        //   printf("MCMC Error 0033 : Wrong parameters for matrix_mult routine");
          return(NULL);
   }
   temp = matrix_blank(mat1->ncol,mat1->ncol);
   if(temp == NULL)
          return(NULL);
   for(i=0;i<temp->nrow;i++)
   {
      for(j=0;j<temp->ncol;j++)
      {
         for(terms=0;terms<mat1->nrow;terms++)
         {
            temp->e[i][j] += mat2->e[terms][i] * mat1->e[terms][j];
         }
      }
   }
   return temp;
}



/*
 * Function performs the same as scale function of R,Splus and returns the result
 if mean=1 removes the means of each columns otherwise not
 if scale=1  rescales each columns so that their variance is 1 otherwise not 
 */
struct matrix *scale(struct matrix *mat1, int mean, int scale)
{
   double s,v;
   int i,j;
   struct matrix *temp;
   
   temp = matrix_blank(mat1->nrow,mat1->ncol);
   for(j=0;j<temp->ncol;j++)
   { 
   s=0;
   v=0.0;
     for(i=0;i<temp->nrow;i++)
        {
        s+=mat1->e[i][j];
        temp->e[i][j]=mat1->e[i][j];
        }
        
        s=s/(double)temp->nrow;
        
   if(mean==1)
        {
        for(i=0;i<temp->nrow;i++)
        temp->e[i][j]-=s;
        }
   if(scale==1)
        {
        for(i=0;i<temp->nrow;i++)
            v+=(mat1->e[i][j]-s)*(mat1->e[i][j]-s);
                    
        v=sqrt(v)/sqrt((double)temp->nrow);
        for(i=0;i<temp->nrow;i++)
        temp->e[i][j]/=v;
        }
             
 
   }
   return temp;
}





