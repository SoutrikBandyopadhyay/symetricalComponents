#include <stdio.h>
#include <math.h>
#include <string.h>

#include <malloc.h>

float pi = 3.141592654;

typedef struct _cmp{
  float real;
  float imag;
} Complex;

typedef struct _pol{
  float mod;
  float theta;
} Polar;


void complexSetValue(Complex *a, Complex old){
  a->real = old.real;
  a->imag = old.imag;

}

float degToRad(const float angle){
  return 0.01745329252 * angle;
}

float radToDeg(const float angle){
  return angle/0.01745329252;
}

Complex newComplex(const float a, const float b){
  Complex new;
  new.real = a;
  new.imag = b;
  return new;
}

void complexView(const Complex a){
  if(a.imag>=0){
    printf("%0.2f + i%0.2f",a.real,a.imag);
  }else{
    printf("%0.2f - i%0.2f",a.real,-a.imag);
  }
}

Polar newPolar(const float r, const float phi){
  Polar new;
  new.mod = r;
  new.theta = degToRad(phi);
  return new;
}

void polarView(const Polar a){
  float angle;
  angle = radToDeg(a.theta);
  if(angle>180){
    angle = angle-360;
  }
  if(angle<-180){
    angle = angle+360;

  }

  if(angle==0){
    printf("%0.2f",a.mod);
  }else{
    printf("%0.2f/__%0.2f",a.mod,angle);
  }
}

Complex complexNegate(const Complex Z){
  Complex ans;
  ans.real = -Z.real;
  ans.imag = -Z.imag;
  return ans;
}

Complex complexConjugate(const Complex Z){
  Complex ans;
  ans.real = Z.real;
  ans.imag = -Z.imag;
  return ans;
}

Complex complexAdd(const Complex Left, const Complex Right){
  Complex ans;
  ans.real = Left.real+Right.real;
  ans.imag = Left.imag+Right.imag;
  return ans;
}

Complex complexSubtract(const Complex Left, const Complex Right){
  return complexAdd(Left,complexNegate(Right));
}

Complex complexMultiply(const Complex Left, const Complex Right){
  Complex ans;
  ans.real = (Left.real*Right.real)-(Left.imag*Right.imag);
  ans.imag = (Left.real*Right.imag)+(Left.imag*Right.real);
  return ans;
}

Complex complexInverse(const Complex Z){
  Complex Zcon;
  Complex ans;
  float den;
  float denInv;

  Zcon = complexConjugate(Z);
  den = complexMultiply(Z,Zcon).real;

  denInv = 1/den;

  ans.real = Zcon.real*denInv;
  ans.imag = Zcon.imag*denInv;

  return ans;
}

Complex complexDivide(const Complex Left, const Complex Right){
  return complexMultiply(Left,complexInverse(Right));
}

Polar complexToPolar(const Complex Z){
  Polar ans;
  // float pi = 3.141592654;

  ans.mod = sqrt(pow(Z.real,2)+pow(Z.imag,2));
  ans.theta = atan(Z.imag/Z.real);
  if(Z.real<0){
    if(Z.imag>0){
      ans.theta-=pi;
    }else{
      ans.theta+=pi;
    }
  }
  return ans;
}

Complex polarToComplex(const Polar P){
  Complex ans;
  ans.real = P.mod*cos(P.theta);
  ans.imag = P.mod*sin(P.theta);
  return ans;
}


typedef struct _mat{
  int row;
  int col;
  Complex *elements;
} Matrix;

//GLOBAL VARIABLES

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
      polarView(P);
      printf("\t");
    }
    printf("\n");
  }
}

int main(){
  printf("_______________________________________________________________");
  printf("\n\nProgram to compute Symetrical Components of Unsymetrical System");
  printf("\n_______________________________________________________________");

  Matrix A;
  Matrix alphaMatrix;
  Matrix unsymetric;
  float r,theta;
  int i;
  Polar Va,Vb,Vc;
  A.row = 3;
  A.col = 1;

  alphaMatrix.row = 3;
  alphaMatrix.col = 3;

  allocateSpace(&A);
  allocateSpace(&alphaMatrix);

  printf("\n\nEnter the Current/Voltage values for the Unsymetrical System");
  for (i = 0; i < 3; i++) {
    char a[3];
    // a = {"R","Y","B"};
    strcpy(a,"RYB");
    printf("\nEnter the Magnitude for %c phase: ",a[i]);
    scanf("%f",&r);
    printf("\nEnter the Phase Angle(degrees) for %c phase: ",a[i]);
    scanf("%f",&theta);

    Polar V;
    V = newPolar(r,theta);
    setElement(&A,i,0,polarToComplex(V));
  }
  // matrixView(A);
  // matrixPolarView(A);

  Polar alpha_p = {1,(pi*2)/3};
  Complex alpha,alpha2;


  alpha = polarToComplex(alpha_p);
  alpha2 = complexMultiply(alpha,alpha);
  // polarView(complexToPolar(ONE));

  setElement(&alphaMatrix,0,0,ONE);
  setElement(&alphaMatrix,0,1,ONE);
  setElement(&alphaMatrix,0,2,ONE);

  setElement(&alphaMatrix,1,0,ONE);
  setElement(&alphaMatrix,1,1,alpha);
  setElement(&alphaMatrix,1,2,alpha2);

  setElement(&alphaMatrix,2,0,ONE);
  setElement(&alphaMatrix,2,1,alpha2);
  setElement(&alphaMatrix,2,2,alpha);

  unsymetric = multiplyScalar(matrixMultiply(alphaMatrix,A),newComplex(1.0/3.0,0));
  printf("\nThe Values of Components (0 , +ve , -ve) are respectively - \n");
  matrixPolarView(unsymetric);

  return 0;
}
