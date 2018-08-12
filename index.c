#include <stdio.h>
#include <math.h>

#include "library/matrix.c"

int main(){
  Matrix A,B;
  int i,j;
  int q,w,e,r;
  scanf("%d %d",&A.row,&A.col);
  allocateSpace(&A);
  for(i=0;i<A.row;i++){
    for(j=0;j<A.col;j++){
      Complex x;
      scanf("%d %d %d %d",&q,&w,&e,&r);
      x = newComplex(q,w,e,r);
      setElement(&A,i,j,x);
    }
  }

  scanf("%d %d",&B.row,&B.col);
  allocateSpace(&B);
  for(i=0;i<B.row;i++){
    for(j=0;j<B.col;j++){
      Complex y;
      scanf("%d %d %d %d",&q,&w,&e,&r);
      y = newComplex(q,w,e,r);
      setElement(&B,i,j,y);
    }
  }

  matrixView(A);
  printf("\n");
  matrixView(B);
  printf("\n");
  matrixPolarView(matrixMultiply(matrixInverse(A),B));

  return 0;

}
