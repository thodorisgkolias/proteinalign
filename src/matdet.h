struct matrix *cofactor(struct matrix *A, int k)
 {
  struct matrix *B;
  int i, j;    
    B=matrix_blank(A->nrow - 1,A->ncol - 1);
    for (j = 0 ; j < B->nrow ; j++ )
    for (i = 0 ; i < B->nrow ; i++ ){
    if ( i < k )
      B->e[j][i] = A->e[j+1][i];
    else
      B->e[j][i] = A->e[j+1][i+1];
  }
  return(B);
}

double  Det(struct matrix *A ){
  struct matrix *B;
  double det;
  int i;  

  det = 0;

  if (A->nrow == 1)
    det = A->e[0][0];
  else 
    for (i = 0 ; i < A->nrow ; i++ )
    {
      B = cofactor(A,i);
      if (i%2 == 0)
        det = det + A->e[0][i] * Det(B);
      else
        det = det - A->e[0][i] * Det(B);
    matrix_free(B);
    }
return(det);
}

