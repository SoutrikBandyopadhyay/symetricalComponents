// #include <stdio.h>
// #include <math.h>

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
  float pi = 3.141592654;

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





// int main(int argc, char const *argv[]) {
//   Complex a,b;
//   // Polar b;
//   a = newComplex(-2,-3);
//   b = newComplex(15,10);
//   // complexView(complexDivide(a,b));
//   complexView(polarToComplex(complexToPolar(a)));
//   return 0;
// }
