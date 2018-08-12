#include <stdio.h>
#include <math.h>

#include "library/matrix.c"

float pi = 3.141592654;

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
    printf("\nEnter the Magnitude for phase {%d}",i+1);
    scanf("%f",&r);
    printf("\nEnter the Phase Angle(degrees) for phase {%d}",i+1);
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

  // matrixPolarView(alphaMatrix);
  // printf("\n");
  // matrixPolarView(A);
  // matrixView(A);

  unsymetric = multiplyScalar(matrixMultiply(alphaMatrix,A),newComplex(1.0/3.0,0));

  matrixPolarView(unsymetric);

  return 0;
}
