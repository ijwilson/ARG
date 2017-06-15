/** @file
 * Some very simple utility functions (mainly for mathematics)
 */
#ifndef CPPUTIL_H_
#define CPPUTIL_H_
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <sstream>
#include <numeric>
#include <cassert>

#ifndef USE_R
#include "gsl/gsl_cdf.h"
#endif


#ifdef CHECK
#define LOG mylog
#else
#define LOG log
#endif


using std::cerr;
using std::ostringstream;
using std::vector;

enum sex{male,female};


#ifdef WIN32
// const double  M_PI = 3.14159265358979323846;
//const double  M_LN2 =0.69314718055994530942;
double lgamma(double xx);
#endif

// some mathematical constants
// #ifndef USE_R
// const double  M_1_SQRT_2PI=	0.398942280401432677939946059934;
// const double  M_SQRT_32=	5.656854249492380195206754896838;
// #endif
const double LMY_SQRT_2PI =-0.918938533204672669540968854562;

#ifdef CHECK
double mylog(double x);
#define mymalloc malloc;
#endif

/// A very simple error class
class error {
public:
  /** Only consists of a constructor that prints a message and then exits 
   * The extra constructors are templates when we want to pass something 
   * (that can be output to an ostream using <<) to the error              */
  error(const char *text) {
    cerr << text;
    exit(EXIT_FAILURE);
  }
  template<typename T>
  error(const char *text, const T &p1,const char *text2=0) {
    ostringstream oss;
    oss << text  << p1 << text2 << std::endl;
    cerr << oss.str().c_str();
    exit(EXIT_FAILURE);
  }
  template<typename T,typename R>
  error(const char *text, const T p1,const char *text2, const R p2, const char *text3=0) {
    ostringstream oss;
    oss << text << p1 << text2 << p2 << text3 << std::endl;
    cerr << oss.str().c_str();
    exit(EXIT_FAILURE);
  }
};


/// class that simulates a factorial function- but includes a lookup table.
class lfactorial {
public:
  lfactorial() {
    for (int i=0;i<32;i++) a[i]=lgamma(static_cast<double>(i+1));
  }
  double operator()(int n) const {
    if (n<=1) return 0.0;
    if (n<=31) return a[n];
    else return lgamma(static_cast<double>(n+1));
  };
 private:
  double a[32];
};

int pow2(int x);
/** returns the sum of the powers of 2 up to n
    sumpow2(6) = 2+4+8+16+32+64 */
int sumpow2(int x);


// functions
double max(double x, double y);
double max(double *x, int n);
// sorting
void hpsort(unsigned long n, double  ra[]);
void hpsort(unsigned long n, int ra[]);
// indexing
void indexi(unsigned long n, int arr[], int indx[]);
// quantile calculation
void quantiles(double *res,double *x, int nx, double *p,int np);
// mathematical functions
/***************************************************************************/
template<typename T> 
T log_D(const vector<T> &b)
{
  double sum=0.0,temp=0.0;
  for (size_t i=0;i<b.size();i++) {
    sum += b[i];
    temp += lgamma(b[i]);
  }
  return temp-lgamma(sum);
}
double log_D(const double *b, int n);

// functionoids
extern lfactorial lfactrl;

#ifndef USE_R
template <typename T>
double cumnorm(T ox, double mean, double sd) {
  return gsl_cdf_gaussian_P(ox-mean,sd) ;
}
#endif


template <typename T>
T sumd(T* a,int len) {
  T tmp(0.0);
  for (int i=0;i<len;i++) tmp+=a[i];
  return tmp;
}

/** is an integer the square of a natural number */
bool IsSquare(int n, int &siz);


//two pass mean and variance calculation
template <typename T>
void avevar(const std::vector<T> &xx, double &mean, double &v) {
  mean=std::accumulate(xx.begin(),xx.end(),0.0)/static_cast<double>(xx.size());
  double ep=0.0,var=0.0;
  for (size_t ii=0;ii<xx.size();ii++) {
    double s=static_cast<double>(xx[ii])-mean;
    ep+=s;
    var+=s*s;
  }
  v = (var-ep*ep/static_cast<double>(xx.size()))/static_cast<double>(xx.size()-1);

}

//two pass mean and variance calculation
template <typename T>
double sample_mean(const std::vector<T> &xx) {
  return std::accumulate(xx.begin(),xx.end(),0.0)/static_cast<double>(xx.size());
}

//two pass mean and variance calculation
template <typename T>
void wavevar(const std::vector<T> weights, const std::vector<T> &xx, double &mean, double &v,bool sub1=false) {
  assert(weights.size()==xx.size());
  double sumw = std::accumulate(weights.begin(),weights.end(),0.0);
  double sum = 0.0;
  for (size_t ii=0;ii<xx.size();ii++) sum += weights[ii]*xx[ii];
  
  mean=sum/sumw;

  double ep=0.0,var=0.0;
  for (size_t ii=0;ii<xx.size();ii++) {
    double s=static_cast<double>(xx[ii])-mean;
    ep+=weights[ii]*s;
    var+=weights[ii]*s*s;
  }
  if (!sub1) v = (var-ep*ep/sumw)/sumw;  // note the lack of a -1 below!
  else v = (var-ep*ep/sumw)/(sumw-1.0); 
}
      
// euclidean distance
template <typename T>
double EuclideanDistance(const std::vector<T> &xx, const std::vector<T> &yy) {
  assert(xx.size()==yy.size());
  double dis = 0.0;
  for (size_t ii=0;ii<xx.size();ii++) {
      dis += (xx[ii]-yy[ii])*(xx[ii]-yy[ii]);
  }
  return sqrt(dis);
}
      
// rescale
template <typename T>
  void rescale( std::vector<T> &xx, const std::vector<std::pair<T,T> > &s) {
  for (size_t ii=0;ii<xx.size();ii++) {
      xx[ii] -= s[ii].first;
      xx[ii] /= s[ii].second;
    }
 }


#endif
