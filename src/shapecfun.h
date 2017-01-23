/*
These is the library of shape functions named as the corresponding functions at R package for shapes
Commends updated 10/04/2004
 
 */

/*
  helmert matrix with dimension nx(n+1)
 */

struct matrix *helmert(int n)
{             
int i,j;
struct matrix *temp;
temp = matrix_blank(n,n+1);
for (j=0;j<n;j++)
{
for (i=0;i<j+1;i++)
{
temp->e[j][i]=-1/sqrt((double)(j+1)*(j+2));
}
temp->e[j][j+1]=(double)(j+1)/sqrt((double)(j+1)*(j+2));
}
return(temp);                              
}
 
 

/*
 * Function creates the centralising matrix with dimension nxn
 */

struct matrix *centralising(int n)
{             
int i,j;
struct matrix *temp;
temp = matrix_blank(n,n);
for (j=0;j<n;j++)
{
for (i=0;i<n;i++)
{

temp->e[j][i]=-1/(double)n;
if (i==j) 
  temp->e[j][i]=1-1/(double)n;
  
}
}

return(temp);                              
}


/*
 * Function calculates the norm of the  matrix and returns the result
 */

double norm(struct matrix *mat)
{        
int i,j;     
double s;
s=0;
   for (i=0;i<mat->nrow;i++)
   {
      for (j=0;j<mat->ncol;j++)
      {
/*
 * Copy elements of matrix
 */
         s+=mat->e[i][j]* mat->e[i][j];
      }
   }
return(sqrt(s));                              
}
 
/*
 * Function creates the preshape matrix ((k-1)xk) after the helmert transformation  
 is carried out. This does not involve matrix multplication which is useful if  
 k is too large
 */
struct matrix *preshape_transform(struct matrix *mat,int rescale)
{ 
int i,j;
struct matrix *temp, *temp1;
double s; 
temp = matrix_blank(mat->nrow-1,mat->ncol); 
 
for (j=0;j<mat->ncol;j++)
{   
s=mat->e[0][j];
   for (i=0;i<temp->nrow ;i++)
   {
   temp->e[i][j]=1/sqrt((double)(i+1)*(i+2))*((double)(i+1)*mat->e[i+1][j]-s);
   s+=mat->e[i+1][j]; 
   }

}
 
if (rescale==1)
{  
temp1=temp;
temp=ms_mult(temp1,1/norm(temp));
matrix_free(temp1);
}
 
return(temp);                              
 
} 
 
 /*
 * Function creates the the icon of the corresponding preshape mat((k-1)xm).
   This does not involve matrix inverse helmert matrix multplication which is useful if k is too large
 */
struct matrix *preshapetoicon_transform(struct matrix *mat)
{ 
int i,j,s;
struct matrix *temp; 
temp = matrix_blank(mat->nrow+1,mat->ncol);  
for (j=0;j<mat->ncol;j++)
{   
   for (i=0;i<temp->nrow;i++)
   {
   if(i==0)
   temp->e[i][j]=0;
   else
   temp->e[i][j]=sqrt((double)(i)/(double)(i+1))*mat->e[i-1][j];
     if(i<mat->nrow)
     {  
      for (s=i;s<temp->nrow-1;s++)
        temp->e[i][j]-=1/sqrt((double)(s+1)*(s+2))*mat->e[s][j];     
     }     
   }
  
}
return(temp);                               
} 

/*
 * Function produces the pre-shape of configuration matrix of dimension kxm and 
 returns the result 
 if rescale is 1 than the rescaling is done
 */

struct matrix *preshape(struct matrix *mat,int rescale)
{
struct matrix *temp,*temp1;
temp1=helmert(mat->nrow-1); 
temp=matrix_mult(temp1,mat);
matrix_free(temp1);
if (rescale==1)
{  
temp1=temp;
temp=ms_mult(temp1,1/norm(temp));
matrix_free(temp1);
}
return(temp);
}


/*
 * Function produces the riemannian shape distance of 2 pre-shapes k-1xm and 
 returns the result
 */

double riemdist_preshape(struct matrix *mat1,struct matrix *mat2 )
{     
int i;
struct matrix *q,*v,*temp;
struct vector *w;
double s;                           
temp=matrix_transpose(mat1);
q=matrix_mult(temp,mat2);
matrix_free(temp);
temp=matrix_transpose(q);
v=q;
q=matrix_mult(temp,v);   
matrix_free(temp);
matrix_free(v);
w=vector_blank(q->ncol);   
v=matrix_blank(q->nrow,q->ncol);   
dsvdcmp(q->e,q->nrow,q->ncol,w->e,v->e);
s=0;
for (i=0;i<w->nrow;i++)
{
/*printf("%2.5lf \n",w->e[i]);*/
s+=sqrt(w->e[i]);
}

temp=matrix_transpose(mat1);
if (Det(matrix_mult(temp,mat2))<0)
   {
   /* s-=2*sqrt(w->e[w->nrow-1]);*/
    s-=2*sqrt(vector_min(w));
  /*  printf("%2.5lf \t\n",vector_min(w));*/
        }
    
matrix_free(temp);
matrix_free(q);
matrix_free(v);
vector_free(w);  
return(acos((s) < (1) ? (s) : (1)));
}



/*
 * Function produces the riemannian shapes&reflection distance of 2 pre-shapes k-1xm and 
 returns the result
 */

double riemdist_preshape_reflect(struct matrix *mat1,struct matrix *mat2 )
{     
int i;
struct matrix *q,*v,*temp;
struct vector *w;
double s;                           
temp=matrix_transpose(mat1);
q=matrix_mult(temp,mat2);
matrix_free(temp);
temp=matrix_transpose(q);
v=q;
q=matrix_mult(temp,v);
matrix_free(temp); 
matrix_free(v);
w=vector_blank(q->ncol);   
v=matrix_blank(q->nrow,q->ncol);   
dsvdcmp(q->e,q->nrow,q->ncol,w->e,v->e);
/*vector_print(w);*/
s=0;
for (i=0;i<w->nrow;i++)
{
/*printf("%2.5lf \n",w->e[i]);*/
s+=sqrt(w->e[i]);
}
matrix_free(q);
matrix_free(v);
vector_free(w);  
return(acos((s) < (1) ? (s) : (1)));
}

/*
 Function produces the riemannian shape distance of 2 configurations kxm and 
 returns the result
 */

double riemdist(struct matrix *mat1,struct matrix *mat2)
{
float s;
struct matrix *tmp1, *tmp2;        
if (1000<mat1->nrow) 
{ 
tmp1=preshape_transform(mat1,1);
tmp2=preshape_transform(mat2,1);          
}  
else 
{
tmp1=preshape(mat1,1);
tmp2=preshape(mat2,1);         
} 
s=riemdist_preshape(tmp1,tmp2);
matrix_free(tmp1);
matrix_free(tmp2);
return(s);
}                             


/*
 Function produces the riemannian shape distance of 2 configurations kxm and 
 returns the result
 */

double riemdist_reflect(struct matrix *mat1,struct matrix *mat2)
{
float s;
struct matrix *tmp1, *tmp2;
if (1000<mat1->nrow) 
{ 
tmp1=preshape_transform(mat1,1);
tmp2=preshape_transform(mat2,1);          
}  
else 
{
tmp1=preshape(mat1,1);
tmp2=preshape(mat2,1);         
} 
s=riemdist_preshape_reflect(tmp1,tmp2);
matrix_free(tmp1);
matrix_free(tmp2);
return(s);
}                             


/*
This function updates the preshape B such that A and B are pocrustes matched. 
*/

void procrustes(struct matrix *A, struct matrix *B)
{
struct vector *w;
w=vector_blank(A->ncol);
struct matrix *temp,*v,*temp1,*temp2;

v=matrix_blank(A->ncol,A->ncol);
temp=matrix_blank(A->ncol,A->ncol);                
temp1=matrix_transpose(B);
temp2=matrix_mult(temp1,A);
matrix_copy_inplace(temp2,temp);
matrix_free(temp1);
matrix_free(temp2);
dsvdcmp(temp->e,A->ncol,A->ncol,w->e,v->e);
temp1=matrix_transpose(v); 
temp2=temp;
temp=matrix_mult(temp2,temp1);
/*matrix_print(temp);*/
matrix_free(temp2);
matrix_free(temp1);
temp1=matrix_mult(B,temp);

matrix_copy_inplace(temp1,B);
/*matrix_print(temp); */
matrix_free(temp1);
matrix_free(temp);
matrix_free(v);
vector_free(w);
}

/*
This function updates the preshape B such that A and B are pocrustes ROTATION matched. 
*/

void procrustes_rot(struct matrix *A, struct matrix *B)
{
int s,i;
float ch1,ch2;
struct vector *w;
w=vector_blank(A->ncol);
struct matrix *temp,*v,*temp1,*temp2;
/*struct matrix *temp1,*temp2;*/

v=matrix_blank(A->ncol,A->ncol);

temp1=matrix_transpose(B);
temp=matrix_mult(temp1,A);
/*matrix_copy_inplace(matrix_mult(matrix_transpose(B),A),temp);*/
matrix_free(temp1);
    //printf("SVD__MAT\n");
    //matrix_print(temp);
dsvdcmp(temp->e,A->ncol,A->ncol,w->e,v->e);
   // printf("____Values of D____\n");
    
   // printf ("d1= % .18f\n", w ->e[0]);
   // printf ("d2= % .18f\n", w ->e[1]);
  //  printf ("d3= % .18f\n", w ->e[2]);
ch1=Det(temp);
ch2=Det(v);

if (ch1*ch2<0)
{
s=vector_min_index(w);
/*printf("s is %d \n",s);*/
for (i=0;i<v->nrow;i++)
     v->e[i][s]*=-1;

}



/*
if (ch1<0 && ch2>0)
{
temp1=temp;  
temp=matrix_mult(temp1,I);
matrix_free(temp1);
}
if (ch1>0 && ch2<0)
{
temp1=v;
v=matrix_mult(temp1,I);
matrix_free(temp1);
}
*/
temp1=temp;
temp2=matrix_transpose(v);
temp=matrix_mult(temp1,temp2);
matrix_free(temp1);
matrix_free(temp2);
/*matrix_print(temp);*/
/*printf("%2.2lf\n",Det(temp));*/
temp1=temp;
temp=matrix_mult(B,temp1);
matrix_free(temp1);
/*matrix_print(matrix_mult(matrix_transpose(A),temp));*/

matrix_copy_inplace(temp,B);

/*matrix_print(temp); */
matrix_free(temp);
matrix_free(v);
vector_free(w);
}



/* function opens a 3d matrix to a 2d one wrt the 3-d dimension */
struct matrix *matrix3d_to_2d(struct matrix3d *mat3d)
{
int i,j,k;
struct matrix *temp;
temp=matrix_blank(mat3d->nrow*mat3d->ncol,mat3d->nblo);
for (k=0;k<mat3d->nblo;k++) 
    {
    for(i=0;i<mat3d->nrow;i++)
       {
       for(j=0;j<mat3d->ncol;j++)
          {
          temp->e[j*(mat3d->nrow)+i][k]=mat3d->e[i][j][k];
          }
       
       }
    }    
return(temp); 
}    
/*
This function calculates the rescaling constants in the procrustes GPA algorithm
the out put is the matrix of dimension nxn whose 1 vector is that of the 
required constants
*/
struct vector *scales_1 (struct matrix3d *mat)
{
double nor;
int i,j;
struct matrix *temp,*v,*temp1;
struct vector *w ,*scals, *z;
temp=matrix3d_to_2d(mat);
nor=norm(temp);

z=vector_blank(mat->nblo);     
for(i=0;i<mat->nblo;i++)
 { 
for(j=0;j<temp->nrow;j++)
  {
  z->e[i]+=temp->e[j][i]*temp->e[j][i];
  }
 }
temp1=covar(temp,1);
/*matrix_print(temp1);*/
matrix_free(temp);
w=vector_blank(mat->nblo);
v=matrix_blank(mat->nblo,mat->nblo);
dsvdcmp(temp1->e,mat->nblo,mat->nblo,w->e,v->e);
/*matrix_print(v);*/
scals=vector_blank(mat->nblo);
/*
for(i=0;i<mat->nblo;i++)
 scals->e[i]=nor/sqrt(z->e[i])*sqrt(temp1->e[i][0]*v->e[i][0]);
 */
j=vector_max_index(w);
for(i=0;i<mat->nblo;i++)
{
 scals->e[i]=nor/sqrt(z->e[i])*sqrt(temp1->e[i][j]*v->e[i][j]);
/* scals->e[i]=sqrt(temp1->e[i][j]*v->e[i][j]);*/
/*printf("\n \t %lf", sqrt(temp1->e[i][j]*v->e[i][j]));*/
}

/*vector_print(z);*/
matrix_free(temp1);
matrix_free(v);
vector_free(w);  
vector_free(z);
return(scals);
}

/*
This function calculates the rescaling constants in the procrustes GPA algorithm
the out put is the matrix of dimension nxn whose 1 vector is that of the 
required constants
*/
struct vector *scales_2 (struct matrix3d *mat)
{
double nor;
int i,j;
struct matrix *temp,*v,*temp1;
struct vector *w ,*scals,*u,*z;
temp=matrix3d_to_2d(mat);
nor=norm(temp);

z=vector_blank(mat->nblo);     
for(i=0;i<mat->nblo;i++)
 { 
for(j=0;j<temp->nrow;j++)
  {
  z->e[i]+=temp->e[j][i]*temp->e[j][i];
  }
 }
temp1=scale(temp,0,1);  
matrix_free(temp);
temp=temp1;

v=matrix_transpose(temp);
/*matrix_free(temp1);*/
temp1=matrix_tranp_mult(v,v);
/*matrix_print(temp1);*/
matrix_free(v);
w=vector_blank(temp->nrow);
v=matrix_blank(temp->nrow,temp->nrow);

dsvdcmp(temp1->e,temp->nrow,temp->nrow,w->e,v->e);
/*matrix_print(temp1);*/
/*matrix_print(v);*/
/*vector_print(w);*/
matrix_free(v);
v=matrix_transpose(temp1);
matrix_free(temp1);
temp1=v;
/*matrix_free(v);*/
v=matrix_mult(temp1,temp);
matrix_free(temp); 
matrix_free(temp1);
u=vector_blank(v->nrow);
for(j=0;j<v->nrow;j++)
      { 
      for(i=0;i<mat->nblo;i++)
            u->e[j]+=v->e[j][i]*v->e[j][i];
      u->e[j]*=w->e[j]; 
      }
/*matrix_print(v);*/
j=vector_max_index(u);

scals=vector_blank(mat->nblo);
float t;
t=0.0;
for(i=0;i<mat->nblo;i++)
{
 scals->e[i]=nor/sqrt(z->e[i])*fabs(v->e[j][i])*sqrt(w->e[j]/u->e[j]);
 /*scals->e[i]=fabs(v->e[j][i])*sqrt(w->e[j]/u->e[j]);*/
/*printf("\n \t %f \n",(v->e[j][i])*(v->e[j][i])*(w->e[j]/u->e[j]));*/
}

/*vector_print(z);*/
/*matrix_print(scale(v,1,1));*/
/*matrix_print(scale(matrix_transpose(v),1,1));*/
matrix_free(v);
vector_free(w);  
vector_free(z);
vector_free(u);
return(scals);
}

struct vector *scales(struct matrix3d *mat)
{
int k,l;
struct vector *z;
k=mat->nblo;
l=mat->ncol*mat->nrow;
/*l=20000;*/
if (k>l)
z=scales_2(mat);
else
z=scales_1(mat);
return(z);
}

/* the same as cmt in R*/

void cnt3old(struct matrix3d *mat)
{
int i;
struct matrix *cent, *temp, *temp1;
cent=centralising(mat->nrow);
for (i=0;i<mat->nblo;i++)
{ 
temp1=get_matrix(mat,i);
temp=matrix_mult(cent,temp1);
update_entry_matrix(mat,temp,i);
matrix_free(temp);
matrix_free(temp1);
}

matrix_free(cent);
}

/* the same as dif in R*/



/* the same as cmt in R but more efficient in high number of landmarks*/

void cnt3(struct matrix3d *mat)
{
int i,j,k;
double s;

for (i=0;i<mat->nblo;i++)
    for (j=0;j<mat->ncol;j++)    
     {
      s=0;
      for (k=0;k<mat->nrow;k++)
        s+=mat->e[k][j][i];

s=s/(double)mat->nrow;
        
for (k=0;k<mat->nrow;k++)
        mat->e[k][j][i]-=s;     

     }
}


/* the same as dif in R*/

double dif (struct matrix3d *mat) 
{
double s,s0;
int i,j,k;
struct matrix *temp,*temp1; 
s=0;
temp=mean3dmatrix(mat);
for (i=0;i<mat->nblo;i++)
    for (j=0;j<mat->ncol;j++)    
      for (k=0;k<mat->nrow;k++)
        s+=(mat->e[k][j][i]-temp->e[k][j])*(mat->e[k][j][i]-temp->e[k][j]);

/*temp1=preshape(temp,2);*/
/*s0=norm(temp1);*/
temp1=matrix_blank(2,2);
s0=norm(temp);
s=s/(s0*s0); 

   
matrix_free(temp);
matrix_free(temp1);
return(s/(double)mat->nblo);
}
/* the same as rgpa in R*/

struct matrix3d *rgpa(struct matrix3d *mat, double toler)
{
double d1,d2;
int i;
struct matrix *temp1,*temp2,*temp3;

d1=10000;
d2=dif(mat);
/*cnt3(mat);*/
temp1=matrix_blank(mat->nrow,mat->ncol);
temp2=matrix_blank(mat->nrow,mat->ncol);
temp3=matrix_blank(mat->nrow,mat->ncol);

        while (fabs(d1-d2)>toler)   
       {   
       d1=d2;
for (i=0;i<mat->nblo;i++)
    {                  
    sum3dmatrix_inplace(mat,temp2);
    get_matrix_inplace(mat,i,temp3);
    matrix_copy_inplace(temp3,temp1);
    matrix_sub_inplace(temp2,temp3);  
  /* procrustes_rot(temp2,temp1); */    
    procrustes(temp2,temp1);      
    update_entry_matrix(mat,temp1,i);
    }
    d2=dif(mat);              
       }
matrix_free(temp1);
matrix_free(temp2);     
matrix_free(temp3);
/*
matrix_free(temp1);
matrix_free(temp2);
*/
return(mat);
}


/* the same as rgpa in R*/
/* but the rotation matching*/
struct matrix3d *rgpa_rot(struct matrix3d *mat, double toler)
{
double d1,d2;
int i;
struct matrix *temp1,*temp2,*temp3;

d1=10000;
d2=dif(mat);
/*cnt3(mat);*/
temp1=matrix_blank(mat->nrow,mat->ncol);
temp2=matrix_blank(mat->nrow,mat->ncol);
temp3=matrix_blank(mat->nrow,mat->ncol);

        while (fabs(d1-d2)>toler)   
       {   
       d1=d2;
for (i=0;i<mat->nblo;i++)
    {                  
    sum3dmatrix_inplace(mat,temp2);
    get_matrix_inplace(mat,i,temp3);
    matrix_copy_inplace(temp3,temp1);
    matrix_sub_inplace(temp2,temp3);
        //printf("Matrix1\n");
        //matrix_print(temp2);
        //printf("Matrix2\n");
        //matrix_print(temp1);
  procrustes_rot(temp2,temp1);      
  /*  procrustes(temp2,temp1); */     
    update_entry_matrix(mat,temp1,i);
    }
    d2=dif(mat);              
       }
matrix_free(temp1);
matrix_free(temp2);     
matrix_free(temp3);
/*
matrix_free(temp1);
matrix_free(temp2);
*/
return(mat);
}

/* the same as sgpa in R*/

void sgpa(struct matrix3d *mat)
{
int i;
struct vector *z;
struct matrix *temp3, *temp;
z=scales(mat);
for (i=0;i<mat->nblo;i++)
{
temp=get_matrix(mat,i);
temp3=ms_mult(temp,z->e[i]);  
update_entry_matrix(mat,temp3,i);
matrix_free(temp3);
matrix_free(temp);
}
vector_free(z);
/*matrix_free(temp3);*/
}




/*
Calculates the GPA in general
*/

struct matrix3d *matrix3d_proc_algor(struct matrix3d *orig,int rescale,double tol1,double tol2)
{
double d1,d2,rho;
int i;
struct matrix3d *mat;
struct matrix *temp, *temp1,*temp2,*temp3,*temp4; 
struct vector *z;
mat=orig;
cnt3(mat);
temp=matrix_blank(mat->nrow,mat->ncol);
/*temp=mean3dmatrix(mat); */
/*temp1=matrix_blank(mat->nrow,mat->ncol);*/
rho=1.52;
while (rho>tol2)
{   
matrix_free(temp);
temp=mean3dmatrix(mat);
d1=10000;    
d2=norm(temp)*(double)mat->nblo;          
       while (d1-d2>tol1)   
       { 
       d1=d2;
for (i=0;i<mat->nblo;i++)
    {      
      /* get_matrix_inplace(mat,i,temp1); */             
   /*
    matrix_copy_inplace(ms_mult(mean3dmatrix(mat),(double)mat->nblo),temp2);
    matrix_free(temp);
    */
  /*  temp2=ms_mult(mean3dmatrix(mat),(double)mat->nblo);*/
       temp2=matrix_blank(mat->nrow,mat->ncol);  
       sum3dmatrix_inplace(mat,temp2);
       temp3=get_matrix(mat,i);
       temp4=matrix_sub(temp2,temp3);
       matrix_free(temp2);             
       procrustes(temp4,temp3);       
       update_entry_matrix(mat,temp3,i);    
       matrix_free(temp3);
   /* matrix_free(temp1);*/ 
       matrix_free(temp4);    
    }                
    temp2=mean3dmatrix(mat);
    d2=norm(temp2)*(double)mat->nblo;       
    matrix_free(temp2); 
/*temp2=mean3dmatrix(mat);*/    
/*temp=mean3dmatrix(mat);*/        
/*    matrix_copy_inplace(temp2,temp); */
    /*matrix_free(temp2); */ 
        }
        
    if (rescale==1)
   {
  z=scales(mat);
for (i=0;i<mat->nblo;i++)
                {
                temp4=get_matrix(mat,i);
                temp3=ms_mult(temp4,z->e[i]);  
                update_entry_matrix(mat,temp3,i);
                matrix_free(temp3);
                matrix_free(temp4);
                }
                
                vector_free(z);

    }                              
temp1=mean3dmatrix(mat);    
rho=riemdist(temp,temp1);
/*temp=mean3dmatrix(mat);*/        
matrix_free(temp1);
/*printf("%lf",rho);*/
}      
/*printf("%lf",rho);*/
/*printf("f");*/
matrix_free(temp);
/*matrix_free(temp1);*/
/*matrix_free(temp2);*/

/*
if (rescale==1)
{
matrix_free(temp3);
vector_free(z);
}
*/
return(mat);
}

/* the algorithm constructed as in R*/
struct matrix3d *GPA(struct matrix3d *orig, int rescale, double tol1, double tol2)
{
double rho;
struct matrix3d *mat;
struct matrix *temp, *temp1;
/*
mat=matrix3d_blank(orig->nrow,orig->ncol,orig->nblo); 
for (i=0;i<mat->nblo;i++)
update_entry_matrix(mat,get_matrix(orig,i),i);
*/
mat=orig;
cnt3(mat);

mat=rgpa(mat,tol1);
temp=mean3dmatrix(mat);  

if (rescale==1)
sgpa(mat);
/*
mat1=mat;
mat=rgpa(mat1,tol1);
matrix3d_free(mat1);
*/
mat=rgpa(mat,tol1);
temp1=mean3dmatrix(mat); 
rho=riemdist(temp,temp1);

while (rho>tol2)
{
matrix_copy_inplace(temp1,temp);
matrix_free(temp1);
if (rescale==1)
sgpa(mat);     
/*
mat1=mat;
mat=rgpa(mat1,tol1);
matrix3d_free(mat1);
*/
mat=rgpa(mat,tol1);
temp1=mean3dmatrix(mat); 
rho=riemdist(temp,temp1);
/*printf("%lf",rho);*/
}

/*printf("%lf",rho);*/
matrix_free(temp);
/*matrix_free(temp1);*/
return(mat);
}




/* the algorithm constructed as in R*/
struct matrix3d *GPA1(struct matrix3d *orig, int rescale, int reflect, double tol1, double tol2)
{
double rho,d1,d2;
struct matrix3d *mat;

/*
mat=matrix3d_blank(orig->nrow,orig->ncol,orig->nblo); 
for (i=0;i<mat->nblo;i++)
update_entry_matrix(mat,get_matrix(orig,i),i);
*/
mat=orig;
cnt3(mat);
 
if (reflect==1)
{
mat=rgpa(mat,tol1);
/*printf("got here");*/
d1=dif(mat); 
/*printf("dif=");
*printf("%lf",d1);
*printf("\n");*/ 
if (rescale==1)
sgpa(mat);
/*
mat1=mat;
mat=rgpa(mat1,tol1);
matrix3d_free(mat1);
*/
mat=rgpa(mat,tol1);
d2=dif(mat);
/*printf("dif="); 
*printf("%lf",d2);
*printf("\n");*/ 
rho=d1-d2;
/*printf("%f",rho);*/
while (rho>tol2)
{
d1=d2;
if (rescale==1)
sgpa(mat);     
/*printf("sgpa done");*/
/*printf("\n");*/
/*
mat1=mat;
mat=rgpa(mat1,tol1);
matrix3d_free(mat1);
*/
mat=rgpa(mat,tol1);
/*printf("rgpa done");
printf("\n");*/

d2=dif(mat);
/*printf("dif="); 
*printf("%lf",d2);
*printf("\n");*/
rho=d1-d2;
}
}

else 
{
mat=rgpa_rot(mat,tol1);
d1=dif(mat); 
/*printf("dif="); 
*printf("%lf",d1);
*printf("\n"); */
if (rescale==1)
sgpa(mat);
/*
mat1=mat;
mat=rgpa(mat1,tol1);
matrix3d_free(mat1);
*/
mat=rgpa_rot(mat,tol1);
d2=dif(mat); 
rho=d1-d2;
/*printf("%f",rho);*/
while (rho>tol2)
{
d1=d2;
if (rescale==1)
sgpa(mat);     
/*
mat1=mat;
mat=rgpa(mat1,tol1);
matrix3d_free(mat1);
*/
mat=rgpa_rot(mat,tol1);
d2=dif(mat); 
/*printf("dif=");
*printf("%lf",d2);
*printf("\n");*/
rho=d1-d2;
}
}
/*printf("%lf",rho);*/
/*matrix_free(temp);*/
/*matrix_free(temp1);*/
return(mat);
}
