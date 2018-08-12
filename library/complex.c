#include "rational.c"
#include "float.c"


typedef struct _cmp{
  Rational real;
  Rational imag;
} Complex;

typedef struct _pol{
  float mod;
  float theta;
} Polar;


Complex newComplex(const int a, const int b, const int c, const int d);
Polar newPolar(const float r, const float phi);

void complexSetValue(Complex * a, Complex old);

Complex complexNegate(const Complex Z);

Complex complexAdd(const Complex Left, const Complex Right);
Complex complexSubtract(const Complex Left, const Complex Right);
Complex complexMultiply(const Complex Left, const Complex Right);
Complex complexInverse(const Complex Z);
Complex complexDivide(const Complex Left, const Complex Right);

Complex complexConjugate(const Complex Z);



Polar complexToPolar(const Complex Z);
Complex polarToComplex(const Polar P);



void complexView(const Complex Z);



//____________Functions______________


void complexSetValue(Complex *a, Complex old){
  a->real = old.real;
  a->imag = old.imag;

}

float degToRad(const float angle){
  return 0.01745329252 * angle;
}

Complex newComplex(const int a, const int b, const int c, const int d){
  Complex new;
  new.real = newRational(a,b);
  new.imag = newRational(c,d);
  return new;
}

Polar newPolar(const float r, const float phi){
  Polar new;
  new.mod = r;
  new.theta = degToRad(phi);
  return new;
}


Complex complexNegate(const Complex Z){
  Complex ans;
  Rational m;
  m.n = -1;
  m.d = 1;
  ans.real = rationalMultiply(Z.real,m);
  ans.imag = rationalMultiply(Z.imag,m);
  return ans;
}

Complex complexConjugate(const Complex Z){
  Complex ans;
  Rational m;
  m.n = -1;
  m.d = 1;
  ans.real = Z.real;
  ans.imag = rationalMultiply(Z.imag,m);
  return ans;

}


Complex complexAdd(const Complex Left, const Complex Right){
  Complex ans;
  ans.real = rationalAdd(Left.real,Right.real);
  ans.imag = rationalAdd(Left.imag,Right.imag);
  return ans;
}

Complex complexSubtract(const Complex Left, const Complex Right){
  return complexAdd(Left,complexNegate(Right));
}

Complex complexMultiply(const Complex Left, const Complex Right){
  Complex ans;
  ans.real = rationalSubtract(rationalMultiply(Left.real,Right.real),rationalMultiply(Left.imag,Right.imag));
  ans.imag = rationalAdd(rationalMultiply(Left.real,Right.imag),rationalMultiply(Left.imag,Right.real));
  return ans;
}

Complex complexInverse(const Complex Z){
  Complex Zcon;
  Complex ans;

  Rational den;
  Rational denInv;

  Zcon = complexConjugate(Z);
  den = complexMultiply(Z,Zcon).real;

  denInv = rationalInverse(den);

  ans.real = rationalMultiply(Zcon.real,denInv);
  ans.imag = rationalMultiply(Zcon.imag,denInv);

  return ans;
}

Complex complexDivide(const Complex Left, const Complex Right){
  return complexMultiply(Left,complexInverse(Right));
}



Polar complexToPolar(const Complex Z){
  Polar ans;
  ans.mod = sqrt(rationalPow(Z.real,2)+rationalPow(Z.imag,2));
  ans.theta = atan(rationalToFloat(rationalDivide(Z.imag,Z.real)));
  return ans;
}



Complex polarToComplex(const Polar P){
  Complex ans;
  ans.real = floatToRational(P.mod*cos(P.theta));
  ans.imag = floatToRational(P.mod*sin(P.theta));
  return ans;
}


void complexView(const Complex Z){
  if(Z.imag.n == 0){
    printf("{%d/%d}",Z.real.n,Z.real.d);
  }
  else{
    if(Z.imag.n>0){
      printf("{%d/%d + %d/%d i}",Z.real.n,Z.real.d,Z.imag.n,Z.imag.d);
    }
    else{
      printf("{%d/%d - %d/%d i}",Z.real.n,Z.real.d,abs(Z.imag.n),Z.imag.d);
    }
  }
}
