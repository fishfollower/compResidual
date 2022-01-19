#define TMB_LIB_INIT R_init_compResidual
#include <TMB.hpp>

namespace my_atomic {
  /*
   *  Modified from R source pbinom.c
   *  Mathlib : A C Library of Special Functions
   *  Copyright (C) 1998 Ross Ihaka
   *  Copyright (C) 2000-2015  The R Core Team
   *  Copyright (C) 2004-2015  The R Foundation
   */
  template<class T> int R_finite(T x) { return std::isfinite(asDouble(x)); }
  template<class T> int isnan(T x) { return std::isnan(asDouble(x)); }

#undef ML_ERROR
#undef MATHLIB_ERROR
#undef MATHLIB_WARNING
#undef MATHLIB_WARNING2
#undef MATHLIB_WARNING3
#undef MATHLIB_WARNING4
#undef MATHLIB_WARNING5
#undef ML_POSINF
#undef ML_NEGINF
#undef ML_NAN
#undef M_SQRT_2dPI
#undef ISNAN
# define ML_ERROR(x, s) /* nothing */
# define MATHLIB_ERROR(fmt,x) /* nothing */
# define MATHLIB_WARNING(fmt,x) /* nothing */
# define MATHLIB_WARNING2(fmt,x,x2) /* nothing */
# define MATHLIB_WARNING3(fmt,x,x2,x3) /* nothing */
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) /* nothing */
# define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) /* nothing */
#define ML_POSINF	R_PosInf
#define ML_NEGINF	R_NegInf
#define ML_NAN		R_NaN
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#define ISNAN(x) (isnan(x)!=0)

#define ML_ERR_return_NAN return R_NaN

#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0) /* 1 */
# define attribute_hidden __attribute__ ((visibility ("hidden")))


  template<class Float>
  attribute_hidden
  Float pbinom0_raw(Float x, Float n, Float p, Float lower_tail_, Float log_p_)
  {
    int lower_tail = (int)trunc(lower_tail_);
    int log_p = (int)trunc(log_p_);

#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n) || ISNAN(p))
      return x + n + p;
    if (!R_FINITE(n) || !R_FINITE(p)) ML_ERR_return_NAN;

#endif
    // if(R_nonint(n)) {
    //   MATHLIB_WARNING(("non-integer n = %f"), n);
    //   //ML_ERR_return_NAN;
    //   return 0;
    // }
    n = (Float)(int)trunc(n);
    /* PR#8560: n=0 is a valid value */
    if(n < 0 || p < 0 || p > 1) ML_ERR_return_NAN;

    if (x < 0) return R_DT_0;
    //x = floor(x + 1e-7);
    if (n <= x) return R_DT_1;
    return atomic::toms708::pbeta((Float)p, (Float)(x + 1), (Float)(n - x), (int)!lower_tail, (int)log_p);
  }
  
  template<class Float>
  Float pbinom0(Float x, Float n, Float p, Float lower_tail, Float log_p)
  {
    return pbinom0_raw(x,n,p,lower_tail,log_p);
  }

  TMB_BIND_ATOMIC(pbinom1,11100,pbinom0(x[0], x[1], x[2], x[3], x[4]))

}

template<class Type>
Type pbinom(Type x, Type n, Type p, int lower_tail, int log_p){
  CppAD::vector<Type> tx(6);
  tx[0] = x;
  tx[1] = n;
  tx[2] = p;
  tx[3] = lower_tail;
  tx[4] = log_p;
  tx[5] = 0; // extra argument for derivative order
  Type res = my_atomic::pbinom1(tx)[0];
  return res;
}

template<class Type>
vector<int> order(vector<Type> k){
  int n=k.size();
  vector<int> o(n);
  o.setZero();
  int at=-1;
  for(int i=0; i<n;++i){
    if(k(i)>0.5){o(++at) = i;}  
  }
  at=n;  
  for(int i=n-1; i>=0;--i){
    if(k(i)<0.5){o(--at) = i;}  
  }
  return o;
}

template<class Type>
vector<Type> rmultinom(Type N, vector<Type> p)
{
  //multinomial
  int dim = p.size();
  vector<Type> x(dim);
  int Nint = CppAD::Integer(N);
  x.setZero();
  for(int i = 0; i < Nint; i++)
  {
    Type y = runif(0.0,1.0);
    for(int a = 0; a < dim; a++) if(y < p.head(a+1).sum())
    {
      x(a) += 1.0;
      break;
    }
  }
  return x;
}

template <class Type>
Type dmultinom_osa(vector<Type> x, vector<Type> p, data_indicator<vector<Type>, Type> keep, int give_log=0)
{
  vector<Type> k=keep;
  vector<Type> l=keep.cdf_lower;
  vector<Type> h=keep.cdf_upper;
  vector<int> o=order(k);
  x=x(o); p=p(o); k=k(o); l=l(o); h=h(o);
  Type logres=0;
  Type nUnused=sum(x);
  Type pUsed=0;
  Type cdf;
  for(int i=0; i<x.size(); ++i){
    if(i!=(x.size()-1)){
      logres += k(i)*dbinom(x(i),nUnused,p(i)/(Type(1)-pUsed),true);
      cdf = pbinom(x(i),nUnused,p(i)/(Type(1)-pUsed),true,false);
      nUnused -= x(i);
      pUsed += p(i);
    }else{ // last index 
      logres += k(i)*Type(0);
      cdf = Type(1);
    }
    cdf = squeeze(cdf);
    logres += l[i] * log( cdf );       // NaN protected
    logres += h[i] * log( 1.0 - cdf ); // NaN protected
  }
  if(give_log){
    return logres;
  }else{ 
    return exp(logres);
  }
}

template <class Type>
vector<Type> rdirichlet(vector<Type> alpha){
  vector<Type> x=rgamma(alpha,Type(1));
  return x/sum(x);
}

template <class Type>
Type ddirichlet(vector<Type> x, vector<Type> alpha, int give_log=0)
{
  Type logB = lgamma(alpha).sum() - lgamma(alpha.sum());
  Type logres=((alpha-Type(1))*log(x)).sum() - logB;
  if(give_log){
    return logres;
  }else{
    return exp(logres);
  }
}

template <class Type>
Type ddirichlet_osa(vector<Type> x, vector<Type> alpha, data_indicator<vector<Type>, Type> keep, int give_log=0)
{
  vector<Type> k=keep;
  vector<Type> l=keep.cdf_lower;
  vector<Type> h=keep.cdf_upper;
  vector<int> o=order(k);
  x=x(o); alpha=alpha(o); k=k(o); l=l(o); h=h(o);
  
  int n = alpha.size();
  Type cdf;
  Type sx = x.sum();
  Type sa = alpha.sum();
  sa -= alpha(0);
  Type logres=k(0)*dbeta(x(0),alpha(0),sa,true);
  cdf = pbeta(x(0),alpha(0),sa);
  cdf = squeeze(cdf);
  logres += l(0) * log( cdf );       
  logres += h(0) * log( 1.0 - cdf ); 
  
  for(int i=1; i<(n-1); ++i){
    sx -= x(i-1);
    sa -= alpha(i);
    logres += k(i)*(dbeta(x(i)/sx,alpha(i),sa,true)-log(sx));
    cdf = pbeta(x(i)/sx,alpha(i),sa);
    cdf = squeeze(cdf);
    logres += l(i) * log( cdf );       
    logres += h(i) * log( 1.0 - cdf ); 
  }
  logres += k(n-1)*Type(0);
  cdf=Type(1);
  cdf = squeeze(cdf);
  logres += l(n-1) * log( cdf );       
  logres += h(n-1) * log( 1.0 - cdf ); 
  
  if(give_log){
    return logres;
  }else{
    return exp(logres);
  }
}

// simulation function 
template<class Type>
vector<Type> rdirmultinom(Type N, vector<Type> alpha) //dirichlet generated from iid gammas
{
  vector<Type> dp = rdirichlet(alpha);
  vector<Type> obs = rmultinom(N,dp);
  return(obs);
}

//the usual D-M
template<class Type>
Type ddirmultinom(vector<Type> obs, vector<Type> alpha, int do_log)
{
  int dim = obs.size();
  Type N = obs.sum();
  Type phi=sum(alpha);
  Type ll = lgamma(N + 1.0) + lgamma(phi) - lgamma(N + phi);
  for(int a = 0; a < dim; a++) ll += -lgamma(obs(a) + 1.0) + lgamma(obs(a) + alpha(a)) - lgamma(alpha(a));
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type>
Type pbetabinom(Type x, Type N, Type alpha, Type beta, int do_log)
{
  //just sum up the probabilities at i = 0,...,x. Nothing fancy.
  int x_int = CppAD::Integer(x);
  vector<Type> p(2), aa(2);
  aa[0] = alpha;
  aa[1] = beta;
  Type Fx = 0.0;
  vector<Type> obs(2);
  for(int i = 0; i <= x_int; i++){
    obs[0] = i;
    obs[1] = N-i;
    Fx += ddirmultinom(obs, aa, 0); //dbetabinomial, just two categories.
  }
  if(do_log == 1) return log(Fx);
  else return Fx;
}

//the D-M as a series of conditional beta-binomials and added args for osa residuals
template<class Type> 
Type ddirmultinom_osa(vector<Type> obs, vector<Type> alpha, data_indicator<vector<Type>, Type> keep, int do_log = 0)
{

  vector<Type> k=keep;
  vector<Type> l=keep.cdf_lower;
  vector<Type> h=keep.cdf_upper;

  vector<int> o=order(k);
  obs=obs(o); alpha=alpha(o); k=k(o); l=l(o); h=h(o);
 
  int dim = obs.size();
  // Type N = obs.sum(); 
  Type ll = 0.0;
  vector<Type> alphas_a(2), obs_a(2);
  for(int a = 1; a < dim; a++){
    obs_a(0) = obs(a-1);
    obs_a(1) = obs.tail(dim-a).sum();
    alphas_a(0) = alpha(a-1);
    alphas_a(1) = alpha.tail(dim-a).sum();
    ll += k(a-1) * ddirmultinom(obs_a, alphas_a, 1); //beta-binomial, just two categories
    Type cdf = pbetabinom(obs_a(0), obs_a.sum(), alphas_a(0), alphas_a(1),0);
    cdf = squeeze(cdf);
    ll += l(a-1) * log( cdf );       
    ll += h(a-1) * log( 1.0 - cdf ); 
  }
  if(do_log == 1) return ll;
  else return exp(ll);
}

template<class Type>
vector<Type> rlogisticnormal(vector<Type> mu, matrix<Type> S, int do_mult)
{
  using namespace density;
  MVNORM_t<Type> mvnorm(S);
  vector<Type> x = mvnorm.simulate() + mu;
  vector<Type> obs(x.size()+1);
  if(do_mult == 1){
    Type accumulate = 1;
    for(int i = 0; i < x.size(); i++) {
      accumulate *= 1 + exp(x(i));
      obs(i) = exp(x(i))/accumulate;
    }
  }
  else {
    for(int i = 0; i < x.size(); i++) obs(i) = exp(x(i))/(1 + sum(exp(x)));
  }
  obs(obs.size()-1) = obs.head(obs.size()-1).sum();
  return(obs);
}

//the usual logistic normal
template<class Type>
Type dlogisticnormal(vector<Type> obs, vector<Type> mu,  matrix<Type> S, int do_mult, int do_log)
{
  using namespace density;
  MVNORM_t<Type> mvnorm(S);

  vector<Type> x(obs.size()-1);
  if(do_mult == 1){
    x = log(obs.head(obs.size()-1));
    for(int i = 0; i < x.size(); i++) x(i) -= log(1-x.head(i+1).sum());
  }
  else x = log(obs.head(obs.size()-1)) - log(obs(obs.size()-1));
  
  Type nll = mvnorm(x-mu);
  nll += log(obs).sum(); //jacobian

  if(do_log == 1) return -nll;
  else return exp(-nll);
}

//the logistic normal with added args for osa residuals
template<class Type>
Type dlogisticnormal_osa(vector<Type> obs, vector<Type> mu,  matrix<Type> S, vector<Type> keep, int do_mult, int do_log)
{
  using namespace density;
  MVNORM_t<Type> mvnorm(S);

  vector<Type> x(obs.size()-1);
  if(do_mult == 1){
    x = log(obs.head(obs.size()-1));
    for(int i = 0; i < x.size(); i++) x(i) -= log(1-x.head(i+1).sum());
  }
  else x = log(obs.head(obs.size()-1)) - log(obs(obs.size()-1));
  Type nll = mvnorm(x-mu, keep.head(x.size()));
  if(sum(keep)>x.size()) nll += log(obs).sum(); //jacobian, do it only when osa residuals not being calculated?

  if(do_log == 1) return -nll;
  else return exp(-nll);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  Type nll=Type(0);
  DATA_INTEGER(code)
  
  if(code==0){ // multivariate normal

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
  if(code==4){ // Additive logistic normal 

  } 
  if(code==5){ // Multiplicative logistic normal   

  }
  PARAMETER(dummy);
  nll += dummy*dummy;  
  return nll;
}
