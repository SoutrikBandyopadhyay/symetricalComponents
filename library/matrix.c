#include "newComplex.c"
#include <malloc.h>


typedef struct _mat{
  int row;
  int col;
  Complex *elements;
} Matrix;

//GLOBAL VARIABLES

// Complex ZERO = {{0,1},{0,1}};
// Complex ONE = {{1,1},{0,1}};
Complex ZERO = {0,0};
Complex ONE = {1,0};


void allocateSpace(Matrix *A);
void initializeZero(Matrix *A);
void initializeToN(Matrix *A,Complex n);
void initializeUnityMatrix(Matrix *A);

Complex getElement(Matrix A, int i, int j);
void setElement(Matrix *A, int i, int j, Complex n);

Matrix matrixAdd(Matrix A, Matrix B);
Matrix multiplyScalar(Matrix A, Complex k);
Matrix matrixSubtract(Matrix A, Matrix B);
Matrix matrixMultiply(Matrix A, Matrix B);
Matrix transpose(Matrix A);
Complex cofactor(Matrix A, int i, int j);
Complex determinant(Matrix A);

Matrix adjoint(Matrix A);
Matrix matrixInverse(Matrix A);

Matrix matrixRead();

void allocateSpace(Matrix *A){
  A->elements = (Complex *) malloc(A->row * A->col * sizeof(Complex));
}

void initializeZero(Matrix *A){
  int i,j;
  for(i=0; i<A->row; i++){
    for(j=0; j<A->col; j++){
      complexSetValue((A->elements + i* A->col + j),ZERO);
    }
  }
}

void initializeToN(Matrix *A, Complex n){
  int i,j;
  for(i=0; i<A->row; i++){
    for(j=0; j<A->col; j++){
      complexSetValue((A->elements + i* A->col + j),n);
    }
  }
}

void initializeUnityMatrix(Matrix *A){
  int i,j;
  for(i=0; i<A->row; i++){
    for(j=0; j<A->col; j++){
      if(i==j){
        complexSetValue((A->elements + i* A->col + j),ONE);
      }
      else{
        complexSetValue((A->elements + i* A->col + j),ZERO);
      }
    }
  }
}

Complex getElement(Matrix A, int i, int j){
  return *(A.elements + i * A.col + j);
}

void setElement(Matrix *A, int i, int j, Complex n){
  complexSetValue((A->elements + i * A->col + j),n);
}

Matrix matrixAdd(Matrix A, Matrix B){
  if(A.row == B.row && A.col == B.col){
    Matrix C;
    C.row = A.row;
    C.col = A.col;

    allocateSpace(&C);
    initializeZero(&C);

    int i,j;
    Complex x;
    for(i=0;i<A.row;i++){
      for(j=0;j<A.col;j++){
        x = complexAdd(*(A.elements + i * A.col + j),*(B.elements + i * B.col + j));
        setElement(&C,i,j,x);
      }
    }
    return C;
  }
}

Matrix multiplyScalar(Matrix A, Complex k){
  Matrix B;
  B.row = A.row;
  B.col = A.col;
  int i,j;
  Complex x;
  allocateSpace(&B);
  for(i=0;i<A.row;i++){
    for(j=0;j<A.col;j++){
      x = complexMultiply(*(A.elements + i * A.col + j),k);
      setElement(&B,i,j,x);
    }
  }
  return B;
}

Matrix matrixSubtract(Matrix A, Matrix B){
  // Complex NEGATIVEONE = {{-1,1},{0,1}};
  Complex NEGATIVEONE = {-1,0};

  return matrixAdd(A,multiplyScalar(B,NEGATIVEONE));
}

Matrix matrixMultiply(Matrix A, Matrix B){

  if(A.col == B.row){
    Matrix C;
    C.row = A.row;
    C.col = B.col;

    int i,j,k;
    Complex x;
    allocateSpace(&C);
    initializeZero(&C);
    for(i=0;i<A.row;i++){
      for(j=0;j<B.col;j++){
        for(k=0;k<A.col;k++){
          x = complexMultiply(*(A.elements + i * A.col + k),*(B.elements + k * B.col + j));
          setElement(&C,i,j,complexAdd(getElement(C,i,j),x));
        }
      }
    }
    return C;
  }
}

Matrix transpose(Matrix A){
  Matrix B;
  B.row = A.col;
  B.col = A.row;
  allocateSpace(&B);
  int i,j;
  Complex x;
  for(i=0;i<A.row;i++){
    for(j=0;j<A.col;j++){
      x = getElement(A,i,j);
      setElement(&B,j,i,x);
    }
  }
  return B;
}

Complex cofactor(Matrix A, int i, int j){
  Matrix B;
  B.row = A.row - 1 ;
  B.col = A.col - 1 ;
  int x,y,li,lj;
  allocateSpace(&B);
  li = 0;
  lj = 0;
  Complex Z;
  for(x=0;x<A.row;x++){
    if(x!=i){
      for(y=0;y<A.col;y++){
        if(y!=j){
          Z = getElement(A,x,y);
          setElement(&B,li,lj,Z);
          lj++;
        }
      }
      li++;
    }
    lj = 0;
  }

  if((i+j)%2 == 1){
    return complexNegate(determinant(B));
  }
  else{
    return determinant(B);
  }
}

Complex determinant(Matrix A){
  if(A.row == 1){
    return getElement(A,0,0);
  }
  else{
    Complex det;
    det = ZERO;
    Complex x;
    int j;

    for(j=0;j<A.col;j++){
      x = complexMultiply(getElement(A,0,j),cofactor(A,0,j));
      det = complexAdd(det,x);
    }
    return det;
  }
}

Matrix adjoint(Matrix A){
  Matrix B;
  B.row = A.row;
  B.col = A.col;
  allocateSpace(&B);
  int i,j;
  for(i=0;i<A.row;i++){
    for(j=0;j<A.col;j++){
      setElement(&B,i,j,cofactor(A,i,j));
    }
  }
  return B;
}

Matrix matrixInverse(Matrix A){
  return multiplyScalar(transpose(adjoint(A)),complexInverse(determinant(A)));
}

void matrixView(Matrix A){
  int i,j;
  for(i=0;i<A.row;i++){
    for(j=0;j<A.col;j++){
      complexView(*(A.elements + i * A.col + j));
      printf("\t");
    }
    printf("\n");
  }
}

void matrixPolarView(Matrix A){
  int i,j;
  Complex Z;
  Polar P;
  for(i=0;i<A.row;i++){
    for(j=0;j<A.col;j++){
      Z = *(A.elements + i * A.col + j);
      P = complexToPolar(Z);
      // P.theta *= 57.29577951;
      //
      // if(Z.real.n<0 && Z.imag.n<0){
      //   P.theta -= 180;
      // }
      // if(Z.real.n<0 && Z.imag.n>=0){
      //   P.theta += 180;
      // }


      polarView(P);
      printf("\t");
    }
    printf("\n");
  }
}

//
// Matrix matrixRead(){
//   Matrix A,B;
//   int i,j;
//   int q,w,e,r;
//   printf("\nEnter rows and cols\n");
//   scanf("%d %d",&A.row,&A.col);
//   allocateSpace(&A);
//   for(i=0;i<A.row;i++){
//     for(j=0;j<A.col;j++){
//       Complex x;
//       printf("\nEnter element of %d row and %d col\n",i,j);
//
//       scanf("%d %d %d %d",&q,&w,&e,&r);
//       x = newComplex(q,w,e,r);
//       setElement(&A,i,j,x);
//     }
//   }
//   return A;
// }
