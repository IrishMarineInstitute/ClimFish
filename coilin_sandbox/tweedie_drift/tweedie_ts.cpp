#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_IVECTOR(idx); // indexes year
  DATA_INTEGER(Y); // maximum year index
  PARAMETER_VECTOR(logmu);
  PARAMETER(logsdmu);
  PARAMETER(d); // drift term
  PARAMETER(logphi);
  PARAMETER(logitp);
  // process on mu
  Type sdmu = exp(logsdmu);
  Type nll = 0;  
  for(int i=1; i<Y; i++){
    nll -= dnorm(logmu(i), logmu(i-1) + d, sdmu, true);
  }
  vector<Type> mu = exp(logmu);
  Type p = 1.0001 + 0.9999 * invlogit(logitp); // bound 1 < p < 2 for point mass at zero
  Type phi = exp(logphi);  
  for(int i=0; i<y.size(); i++){
    nll -= dtweedie(y(i), mu(idx(i)), phi, p, true);
  }
  return nll;
}
