#define TMB_LIB_INIT R_init_compResidual
#include <TMB.hpp>
#include "../inst/include/distr.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  Type nll=Type(0);
  DATA_INTEGER(code)
  
  if(code==0){ // multivariate-normal
    DATA_VECTOR(obs);
    DATA_VECTOR(mu);
    DATA_VECTOR_INDICATOR(keep,obs);
    DATA_MATRIX(Sigma);
    using namespace density;
    MVNORM_t<Type> mvnorm(Sigma);
    nll += mvnorm(obs-mu, keep);
  }
  if(code==1){ // multinomial
    DATA_INTEGER(dim);    
    DATA_VECTOR(obs);
    DATA_VECTOR(pred);
    DATA_IVECTOR(idx);  
    DATA_VECTOR_INDICATOR(keep, obs);
    vector<Type> p(dim);
    for(int i=0; i<idx.size(); ++i){
      p = pred.segment(idx(i),dim);
      p /= sum(p);
      nll += -dmultinom_osa(vector<Type>(obs.segment(idx(i),dim)),p,keep.segment(idx(i),dim),true);
    }
  }
  if(code==2){ // Dirichlet
    DATA_INTEGER(dim);    
    DATA_VECTOR(obs);
    DATA_VECTOR(alpha);
    DATA_IVECTOR(idx);  
    DATA_VECTOR_INDICATOR(keep, obs);
    for(int i=0; i<idx.size(); ++i){
      nll += -ddirichlet_osa(vector<Type>(obs.segment(idx(i),dim)), vector<Type>(alpha.segment(idx(i),dim)), keep.segment(idx(i),dim), true);
    }
  }
  if(code==3){ // Dirichlet-multinomial  
    DATA_INTEGER(dim);
    DATA_VECTOR(obs);
    DATA_VECTOR(alpha);
    DATA_IVECTOR(idx);
    DATA_VECTOR_INDICATOR(keep,obs);
    for(int i=0; i<idx.size(); ++i){
      nll += -ddirmultinom_osa(vector<Type>(obs.segment(idx(i),dim)), vector<Type>(alpha.segment(idx(i),dim)), keep.segment(idx(i),dim), true);
    }
    
  }
  if(code==4){ // Logistic-normal
    DATA_VECTOR(obs);
    DATA_VECTOR(mu);
    DATA_INTEGER(do_mult); // 0 = Additive logistic-normal, 1 = Multiplicative logistic-normal
    DATA_VECTOR_INDICATOR(keep,obs);
    DATA_MATRIX(Sigma);
    nll += -dlogisticnormal_osa(obs, mu, Sigma, keep, do_mult, true);
  }
  PARAMETER(dummy);
  nll += dummy*dummy;  
  return nll;
}
