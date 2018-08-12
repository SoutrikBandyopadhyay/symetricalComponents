

Rational floatToRational(float f);
float rationalToFloat(Rational R);
float rationalSqrt(const Rational R);

float rationalPow(const Rational R, const float n);

Rational floatToRational(float f){
  Rational ans;
  if(f!=0){
    ans.n = (int)(f*1000000);
    ans.d = 1000000;
  }
  else{
    ans.n = 0;
    ans.d = 1;
  }
  return rationalProper(ans);
}

float rationalToFloat(Rational R){
  return (((float)R.n)/((float)R.d));
}


float rationalSqrt(const Rational R){
  float f;
  f = sqrt(((float)R.n)/((float)R.d));
  return f;
}


float rationalPow(const Rational R, const float n){
  float f;
  f = pow((((float)R.n)/((float)R.d)),n);
  return f;
}
