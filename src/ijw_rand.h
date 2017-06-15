/** @file */
// @file      Time-stamp: <2013-12-10 14:48:58 nijw>
#ifndef GSL_RAND_H__
#define GSL_RAND_H__

#include <cassert>
#include <stdexcept>  // for underflow_error
#include "utilityfunctionals.h"
//#include <Rmath.h>
#include <Rcpp.h>


/// A class to act as a cumulative sum functional
/*template <typename T>
class cumsum {
public:
  cumsum():sum(0){}
  T& operator()(const T &item) {
    sum += item;
    return sum;
  }
private:
  T sum;
};*/




/** The global random number generator             */
class rng;

/** wrapper for the gsl rng class             */
class rng {
public:
  /** this is set for the R random number generator.  Note that
   * there should only be a single member of this class, and that
   * this should be enforced (using a singleton template)        */
  rng() {
    GetRNGstate();
  }
  ~rng() {
    PutRNGstate();
  }
  double next() {
    return unif_rand();
  }
  unsigned long int rint(unsigned long int mx) {
    return (unsigned long)(next()*mx);
  }
  double normal() {
    return norm_rand();
  }
  double rGamma(double a, double b) {
    return R::rgamma(a, b);
  }
  int rBinom(double pp, int n) {
    return static_cast<int>(R::rbinom(static_cast<double>(n),pp));
  }
  int rpoisson(double mu) {
    return int(R::rpois(mu));
  }
  double sexp() {
    return exp_rand();
  }
  double sexp(double mu) {
   return mu*exp_rand();
  }
  int rgeometric(double p) {
    return int(R::rgeom(p));
  }
  double rlnorm(double zeta, double sigma) {
    return R::rlnorm(zeta, sigma);
  }
  std::vector<double> rdirichlet(const std::vector<double> &a) {
    size_t d=a.size();
    std::vector<double> bb(d);
    for (size_t ii=0;ii<d;ii++) bb[ii] = R::rgamma(a[ii],1.0);
    double sm=std::accumulate(bb.begin(),bb.end(),0.0);
    for (size_t ii=0;ii<d;ii++) bb[ii] /= sm;
    return bb;
  }
  unsigned int rhypergeometric(unsigned int n1, unsigned int n2, unsigned int t) {
    return static_cast<unsigned int>(R::rhyper(static_cast<double>(n1),
                                            static_cast<double>(n2),static_cast<double>(t)));
  }
  //
  double operator()(double mn, double mx) {
    return mn+(mx-mn)*next();
  }
  double operator()() {
    return next();
  }
  unsigned long operator()(unsigned long N) {
    return rint(N);
  }

  std::vector<double> normal(unsigned long n) {
    std::vector<double> a(n);
    for (unsigned int i=0;i<n;i++) a[i]= normal();
    return a;
  }

  /** Sample a sorted pair of integers in [from,to] with first>second */
  std::pair<int,int> sample2intsorted(int from, int to) {
    std::pair<int,int> p;
    double len=(double)(to-from+1);
    p.second = from + (int)(next()*len);
    p.first = from + (int)(next()*(len-1.));
    if (p.first>=p.second) {
      p.first++;
    } else {
      std::swap(p.first,p.second);
    }
    return p;
  }
  /** Sample a sorted pair of integers in [0,to) with first>second */
  std::pair<int,int> sample2intsorted(int to) {
    std::pair<int,int> p;
    double len=(double)(to);
    p.second = static_cast<int>(next()*len);
    p.first =  static_cast<int>(next()*(len-1.));
    if (p.first>=p.second) {
      p.first++;
    } else {
      std::swap(p.first,p.second);
    }
    return p;
  }
  /** Sample a pair of (different) integers from  [0,to) */
  std::pair<int,int> sample2int(int to) {
    std::pair<int,int> p;
    p.first = rint(to-1);
    p.second=rint(to-1);
    if (p.first>=p.second) {
      p.first++;
    }
    return p;
  }

 /** produce a sorted random selection k integers from [0:n) */
  std::vector<int> integer_permutations(int k, int n) {
    std::vector<int> a(n);
    for (int i=0;i<n;i++) a[i]=i;
    for (int i=0;i<k;i++) {
      int swp=rint(n);
      int tmp=a[swp];
      a[swp]=a[i];
      a[i]=tmp;
    }
    a.resize(k);
    return a;
  }
 /** produce a random permutation of k integers  from [0:n) (k=1..n) */
  std::vector<int> integer_choose(int k, int n) {
    std::vector<int> a = integer_permutations(k,n);
    sort(a.begin(),a.end());
    return a;
  }

private:
  // functions that are not defined for safety
  rng(const rng &a);
  rng& operator=(const rng &a);
};



template <typename T> size_t gen_from_cump(const T &cp, rng &r);
template <typename T> int gen_from_pr(const T &p, int lst, rng &r);
template <typename T>
int gen_from_pb(const T *const p, int lst, rng &r);

/** Generate from cumulative probabilities in p  */
template <typename T>
size_t gen_from_cump(const T &cp, rng &r)
// a routine to pick a number from 0 to n-1 based on cumulative
// "probabilities" in p
// this only takes deques and vectors...
{
  if (cp.back()<=0.0)
    throw std::runtime_error("sum of probabilities is "
			       "equal to zero in  gen_from_cump");
  double relprob = r();
  size_t where = static_cast<int>(relprob*(cp.size())); // approximate start
  relprob *= cp.back();
  for (;;) {
    if  (relprob <= cp[where]) {
      if (where==0) return 0;
      if (relprob > cp[where-1]) return where;
      else where--;
    } else {
      where++;
      assert(where < cp.size());
      if (where==cp.size()-1) return where;
    }
  }
}
/** Generate from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
int gen_from_p(const T &p, rng &r)
{
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  //std::partial_sum(p.begin(),p.end(),cprob.begin());
  return gen_from_cump(cprob,r);
}
/** Generate n samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
std::vector<int> gen_from_p(int n,const T &p, rng &r)
{
  std::vector<int> ret;
  ret.reserve(n);
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  //std::partial_sum(p.begin(),p.end(),cprob.begin());
  for (int i=0;i<n;i++) ret.push_back(gen_from_cump(cprob,r));
  return ret;
}
/** Generate n samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
void gen_from_pn(int *res,int n,T *first, int len, rng &r)
{
  std::vector<double> cprob(len);
  std::transform(first,first+len,cprob.begin(),cumsum<double>());
  for (int i=0;i<n;i++) res[i] = gen_from_cump(cprob,r);
}

/** Generate 2 samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
std::pair<int,int> gen_2_from_p(const T &p, rng &r)
{
  std::pair<int,int> ret;
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  //std::partial_sum(p.begin(),p.end(),cprob.begin());
  ret.first=gen_from_cump(cprob,r);
  ret.second=gen_from_cump(cprob,r);
  return ret;
}
/** Generate 2 samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
std::pair<int,int> gen_2_from_p_diff(const T &p, rng &r)
{
  std::pair<int,int> ret;
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  //std::partial_sum(p.begin(),p.end(),cprob.begin());
  do {
	ret.first=gen_from_cump(cprob,r);
	ret.second=gen_from_cump(cprob,r);
  } while (ret.first == ret.second);
  return ret;
}
/** Generate n samples from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename INTTYPE, typename T>
void gen_from_p(INTTYPE *a,int n,const T &p, rng &r)
{
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  //std::partial_sum(p.begin(),p.end(),cprob.begin());
  for (int i=0;i<n;i++) a[i]=gen_from_cump(cprob,r);
}
/** Generate from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
int gen_from_p(const T &p, rng &r,double &psample)
{
  std::vector<double> cprob(p.size());
  std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  int i= gen_from_cump(cprob,r);
  if (i>0) psample=(cprob[i]-cprob[i-1])/cprob.back();
  else psample=cprob[i]/cprob.back();
  assert(psample>=0.0&&psample<=1.0);
  return i;
}
/** Generate from a vector of probabilities.
 * psample gives the probability of the sample selected      */
template <typename T>
int gen_from_p(int n,const T &p, rng &r,double &psample)
{
  std::vector<double> cprob(p.size());
    std::transform(p.begin(),p.end(),cprob.begin(),cumsum<double>());
  int i= gen_from_cump(cprob,r);
  if (i>0) psample=(cprob[i]-cprob[i-1])/cprob.back();
  else psample=cprob[i]/cprob.back();
  assert(psample>=0.0&&psample<=1.0);
  return i;
}
/** this assumes that T has a random access iterator */
/** gives values in [frst,last]                      */
template <typename T>
int gen_from_p(const T &p,int frst, int lst, rng &r)
{
  std::vector<double> cprob(lst-frst+1);
    std::transform(p.begin()+frst,p.begin()+lst+1,cprob.begin(),cumsum<double>());
    //std::partial_sum(p.begin()+frst,p.begin()+lst+1,cprob.begin());
  return frst+gen_from_cump(cprob,r);
}
/** this assumes that T has a random access iterator */
/** gives values in [0,last)                      */
template <typename T>
int gen_from_pr(const T &p, int lst, rng &r)
{
  std::vector<double> cprob(lst);
    std::transform(p.begin(),p.begin()+lst,cprob.begin(),cumsum<double>());
    // std::partial_sum(p.begin(),p.begin()+lst,cprob.begin());
  return gen_from_cump(cprob,r);
}
/** this assumes that T has a random access iterator */
/** gives values in [0,last)                      */
template <typename T>
int gen_from_pr(const T &p, int lst, rng &r, double &psample)
{
  std::vector<double> cprob(lst);
  std::transform(p.begin(),p.begin()+lst,cprob.begin(),cumsum<double>());
   //std::partial_sum(p.begin(),p.begin()+lst,cprob.begin());
  int i = gen_from_cump(cprob,r);
  if (i>0) psample=(cprob[i]-cprob[i-1])/cprob.back();
  else psample=cprob[i]/cprob.back();
  assert(psample>=0.0&&psample<=1.0);
  return i;
}

/** generate proportional to numbers in [first,end)
 * don't need to be random access                           */
template <typename ITOR>
int gen_from_p(ITOR first, ITOR it_end, rng &r, double &psample)
{
  std::vector<double> cprob(std::distance(first,it_end));
  //std::partial_sum(first,it_end,cprob.begin());
  std::transform(first,it_end,cprob.begin(),cumsum<double>());
  int i= gen_from_cump(cprob,r);
  if (i>0) psample=(cprob[i]-cprob[i-1])/cprob.back();
  else psample=cprob[i]/cprob.back();
  assert(psample>=0.0&&psample<=1.0);
  return i;
}
/** this assumes that T has a random access iterator */
/** gives values in [0,last]                         */
/** I would like to get rid .....                    */
template <typename T>
int gen_from_p(const T *p, int lst, rng &r)
{
  std::vector<double> cprob(lst+1);
  std::transform(p,p+lst+1,cprob.begin(),cumsum<double>());
   //std::partial_sum(p,p+lst+1,cprob.begin());
  return gen_from_cump(cprob,r);
}
/** this assumes that T has a random access iterator */
/** gives values in [0,last)                         */
/** I would like to get rid .....                    */
template <typename T>
int gen_from_pb(const T *const p, int lst, rng &r)
{
  std::vector<double> cprob(lst);
   std::transform(p,p+lst,cprob.begin(),cumsum<double>());
   //std::partial_sum(p,p+lst,cprob.begin());
  return gen_from_cump(cprob,r);
}
/** generate proportional to numbers in [first,end)
 * don't need to be random access                           */
template <typename ITOR>
int gen_from_p(ITOR first, ITOR it_end, rng &r)
{
  std::vector<double> cprob(std::distance(first,it_end));
   std::transform(first,it_end,cprob.begin(),cumsum<double>());
   //std::partial_sum(first,it_end,cprob.begin());
  return gen_from_cump(cprob,r);
}
/** generate proportional to numbers in [first,end)
 * don't need to be random access                           */
// template <typename ITOR>
// int gen_from_p(ITOR first, ITOR it_end, rng &r, double &prob)
// {
//   std::vector<double> cprob(std::distance(first,it_end));
//   std::transform(first,it_end,cprob.begin(),cumsum<double>());
//   int i=gen_from_cump(cprob,r);
//   prob=cprob[i]/cprob.back();
//   return i;
// } 
/** template class for permuting a set of T's              */
template <class T>
void permute(std::vector<T> &x,rng &r)
{
  for (size_t i=0;i<x.size();i++) {
    unsigned long which = r.rint(x.size());
    std::swap(x[which],x[i]);
  }
}
/** template class for permuting a set of T's              */
template <class T>
void permute(T &x,rng &r)
{
  for (size_t i=0;i<x.size();i++) {
    unsigned long which = r.rint(x.size());
    std::swap(x[which],x[i]);
  }
}

/** template class for permuting a set of T's              */
template <class T>
void permute(T *x,int size,rng &r)
{
  for (int i=0;i<size;i++) {
    int which = r.rint(size);
    std::swap(x[which],x[i]);
  }
}

/** template class for permuting a set of T's              */
template <class T,class Y>
void permute2(T *x1,Y *x2,int size,rng &r)
{
  for (int i=0;i<size;i++) {
    int which = r.rint(size);
    std::swap(x1[which],x1[i]);
    std::swap(x2[which],x2[i]);
  }
}



#endif
