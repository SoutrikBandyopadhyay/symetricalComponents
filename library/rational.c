// #include <stdio.h>
// #include <math.h>
#include <stdbool.h>

typedef struct _rat {
  int n;
  int d;
} Rational;




//______________Helper Function___________
static int hcf(const int a , const int b);

Rational newRational(const int n, const int d);


Rational rationalNegate(const Rational R);

Rational rationalProper(const Rational R);

Rational rationalAdd(const Rational Left, const Rational Right);
Rational rationalSubtract(const Rational Left, const Rational Right);
Rational rationalMultiply(const Rational Left, const Rational Right);
Rational rationalInverse(const Rational R);
Rational rationalDivide(const Rational Left, const Rational Right);



void rationalView(Rational R);

// bool rationalEquals(const Rational Left, const Rational Right);
// bool rationalLessThan(const Rational Left, const Rational Right);
// bool rationalLessThanOrEqual(const Rational Left, const Rational Right);
// bool rationalGreaterThan(const Rational Left, const Rational Right);
// bool rationalGreaterThanOrEqual(const Rational Left, const Rational Right);

Rational newRational(const int n, const int d){
  Rational new;
  new.n = n;
  new.d = d;
  return rationalProper(new);
}


Rational rationalNegate(const Rational R){
  Rational ans;
  ans.n = -1 * R.n;
  ans.d = R.d;
  ans = rationalProper(ans);
  return ans;
}

Rational rationalProper(const Rational R){
  Rational ans;
  ans.n = R.n/hcf(R.n,R.d);
  ans.d = R.d/hcf(R.n,R.d);
  if(ans.d < 0){
    ans.n *= -1;
    ans.d *= -1;
  }
  if(ans.n == 0){
    ans.d = 1;
  }
  return ans;
}

static int hcf(const int a , const int b){
    int i, gcd;
    gcd = 1;

    for(i=1; i <= abs(a) && i <= abs(b); ++i)
    {
      if(a%i==0 && b%i==0){
        gcd = i;
      }
    }
    return gcd;
}


Rational rationalAdd(const Rational Left, const Rational Right){
  Rational ans;
  ans.d = Left.d * Right.d;
  ans.n = Left.n * Right.d + Right.n * Left.d;

  return rationalProper(ans);
}

Rational rationalSubtract(const Rational Left, const Rational Right){
  return rationalAdd(Left,rationalNegate(Right));
}

Rational rationalMultiply(const Rational Left, const Rational Right){
  Rational ans;
  ans.n = Left.n * Right.n;
  ans.d = Left.d * Right.d;
  return rationalProper(ans);
}

Rational rationalInverse(const Rational R){
  Rational ans;
  ans.n = R.d;
  ans.d = R.n;
  return rationalProper(ans);
}

Rational rationalDivide(const Rational Left, const Rational Right){
  return rationalMultiply(Left,rationalInverse(Right));
}


void rationalView(Rational R){
  if(R.n != 0 && R.d != 1){
    printf("\n%d/%d",R.n,R.d);
  }
  else{
    printf("\n%d",R.n);
  }
}
