#ifndef ACL_MODEL_MATH_HPP
#define ACL_MODEL_MATH_HPP
// TMB's normal CDF has analytic AD derivatives. Evaluate right tails via
// symmetry to avoid subtracting two probabilities rounded to one.
template<class Type>
Type normal_interval(Type upper, Type lower) {
  Type left = pnorm(upper) - pnorm(lower);
  Type right = pnorm(-lower) - pnorm(-upper);
  Type probability = CppAD::CondExpGt(lower, Type(0), right, left);
  return CppAD::CondExpLt(probability, Type(0), Type(0), probability);
}
template<class Type>
Type positive_probability(Type p) {
  return CppAD::CondExpLt(p, Type(1e-300), Type(1e-300), p);
}
template<class Type>
matrix<Type> length_age_key(vector<Type> borders, vector<Type> ages,
                           Type Linf, Type vbk, Type t0, Type cv) {
  int L = borders.size() + 1, A = ages.size();
  matrix<Type> pla(L, A);
  for (int a = 0; a < A; ++a) {
    Type mean = Linf * (Type(1) - exp(-vbk * (ages(a) - t0)));
    Type sd = mean * cv;
    vector<Type> z = (borders - mean) / sd;
    pla(0,a) = positive_probability(pnorm(z(0)));
    for (int l = 1; l < L-1; ++l)
      pla(l,a) = positive_probability(normal_interval(z(l), z(l-1)));
    pla(L-1,a) = positive_probability(pnorm(-z(L-2)));
    Type total = pla.col(a).sum();
    pla.col(a) /= total;
  }
  return pla;
}
#endif
