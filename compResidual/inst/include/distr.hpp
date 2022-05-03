namespace dists {
// TODO: Move to TMB
template<class Type>
Type pbinom(Type x, Type n, Type p) {
  return 1. - pbeta(p, x+1, n-x);
}
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
  Type nUnused=asDouble(sum(x));
  Type pUsed=0;
  Type cdf;
  for(int i=0; i<x.size(); ++i){
    if(i!=(x.size()-1)){
      logres += k(i)*dbinom(x(i),nUnused,p(i)/(Type(1)-pUsed),true);
      cdf = dists::pbinom(x(i),nUnused,p(i)/(Type(1)-pUsed));
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
  Type sx = 1; // was: x.sum();
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
  Type ll = 0.0;
  vector<Type> alphas_a(2), obs_a(2);
  Type alp_sum = alpha.sum();
  Type obs_sum = asDouble(obs.sum());
  for(int a = 0; a < dim-1; a++){
    obs_sum -= obs[a];
    alp_sum -= alpha[a];
    obs_a(0) = obs(a);
    obs_a(1) = obs_sum;
    alphas_a(0) = alpha(a);
    alphas_a(1) = alp_sum;
    ll += k(a) * ddirmultinom(obs_a, alphas_a, 1); //beta-binomial, just two categories
    Type cdf = pbetabinom(obs_a(0), obs_a.sum(), alphas_a(0), alphas_a(1), 0);
    cdf = squeeze(cdf);
    ll += l(a) * log( cdf );
    ll += h(a) * log( 1.0 - cdf );
  }
  if(do_log == 1) return ll;
  else return exp(ll);
}
