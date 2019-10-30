/** @file */
#ifndef GSL_DISTRIBUTIONS_H__
#define GSL_DISTRIBUTIONS_H__

#include "ijw_rand.h"
#include <Rmath.h>
//#include "util.H"     // for warning
using std::vector;
using std::ostream;

/** a wrapper class around the gsl distributions                */
/** First a generic distribution
 * which has a distribution on T with parameters PARAMETERS     */
template<typename T,typename PARAMETERS>
class distribution {
public:
  virtual T sample(rng &r) const =0;
  virtual double log_pdf(const T &x) const =0;
  virtual distribution *clone() const =0;
  virtual ~distribution(){};
  virtual ostream& print(ostream &o) const =0;
  const PARAMETERS &parameters() const {
    return _p;
  }
protected:
  PARAMETERS _p;
};

/** A class for discrete distribtions of a single parameter that
* can be represented by an integer that have parameters represented
* by pointers to doubles -- i.e. an array of doubles                 */
class ddistribution:
  public distribution<int,double* > {
public:
  ddistribution(int npars) {
    _p = new double(npars);
  }
  ddistribution(double par1, double par2) {
    _p = new double[2];
    _p[0]=par1;
    _p[1]=par2;
  }
  ddistribution(double par) {
    _p = new double[1];
    _p[0]=par;
  }
  ddistribution() { // for improper distribution
  }
  virtual int sample(rng &r) const =0;
  virtual vector<int> sample(rng &r,int n) const =0;
  virtual bool proper() const {
    return true;
  }
  virtual bool isconstant() const {
    return false;
  }
  ~ddistribution() {
    delete[] _p;
  }
  virtual double log_posterior(int x1, int x2) const {
    return 0.0;
  }
  virtual double log_pdf(const int &x) const =0;
};
/** Class for a univariate continuous distribution on the real line     */
class ctsdistribution:
  public distribution<double,double* > {
public:
  ctsdistribution(int npars) {
    _p = new double(npars);
  }
  ctsdistribution(double par1, double par2) {
    _p = new double[2];
    _p[0]=par1;
    _p[1]=par2;
  }
  ctsdistribution(double par) {
    _p = new double[1];
    _p[0]=par;
  }
  ctsdistribution() { // for improper distribution
  }
  virtual std::vector<double> parameters() const=0;
  virtual double sample(rng &r) const =0;
  virtual ctsdistribution *clone() const =0;
  /** try to sample                          */
  void sample(rng &r,vector<double> &x) const  {
    for (size_t i=0;i<x.size();i++) {
      x[i]= sample(r);
     }
  }
  /** Sample to a vector                     */
  vector<double> sample(rng &r, int n) const {
    vector<double> tmp(n);
    for (int i=0;i<n;i++) tmp[i]=sample(r);
    return tmp;
  }
  virtual bool proper() const {
    return true;
  }
  virtual bool isconstant() const {
    return false;
  }
  virtual double max() const {
    return -1E20;
  }
  virtual double min() const {
    return 1E20;
  }
  ~ctsdistribution() {
    delete[] _p;
  }
  virtual double log_posterior(int x1, int x2) const {
    return 0.0;
  }
  virtual double log_pdf(const double &x) const =0;
};

// template<typename T>
// class Constant:
//   public distribution<T,T> {
// public:
//   Constant(const T &val) {
//     _p = new T[1];
//     _p[0]=val;
//   }
//   ~Constant(){
//     delete[] val;
//   };
//   T sample(rng &r) {
//     return _p[0];
//   }
//   virtual bool isconstant() const {
//     return true;
//   }
//   Constant *clone() {
//     Constant *tmp=new Constant(_p[0]);
//     return tmp;
//   }
//   double log_pdf(const double &x) const {
//     assert(fabs(x-_p[0])<1E-8);
//     return 0.0;
//   }
// }
/** an improper distribution                                             */
class improperuniform:
  public ctsdistribution {
public:
  improperuniform()
    :ctsdistribution()  {
  }
  ~improperuniform(){};
  double sample(rng &r) const  {
    warning("cannot generate from an improper distribution\n");
    return 1.0;
  }
  improperuniform *clone() const {
    improperuniform *tmp=new improperuniform();
    return tmp;
  }
  virtual bool proper() const  {
    return false;
  }
  double log_pdf(const double &x) const {
    // note that a specifically written version of this for
    // logarithms may be preferable
    return 0.0;
  }
  ostream& print(ostream &o)const {
    return o << "uniform";
  }
  std::vector<double> parameters() const {
    return std::vector<double>(0);
  }

};
/** a positive improper distribution                                             */
class positiveuniform:
  public ctsdistribution {
public:
  positiveuniform()
    :ctsdistribution()  {
  }
  ~positiveuniform(){};
  double sample(rng &r) const {
    throw std::domain_error("cannot generate from an improper distribtion");
    return 1.0;
  }
  vector<double> sample(rng &r,int n) const {
    throw std::domain_error("cannot generate from an improper distribtion");
    return vector<double>(0);
  }
  positiveuniform *clone() const {
    positiveuniform *tmp=new positiveuniform();
    return tmp;
  }
  virtual bool proper() const  {
    return false;
  }
  double log_pdf(const double &x) const {
    // note that a specifically written version of this for
    // logarithms may be preferable
    assert(x>=0.0);
    return 0.0;
  }
 ostream& print(ostream &o) const {
   return o << "positiveuniform";
  }
  std::vector<double> parameters() const {
    return std::vector<double>(0);
  }
};
/** Constant values                                                */
class Constant:
  public ctsdistribution {
public:
  Constant(double val)
    :ctsdistribution(val)  {
  }
  ~Constant(){};
  double sample(rng &r) const {
    return _p[0];
  }
  vector<double> sample(rng &r,int n) const {
    return vector<double>(n,_p[0]);
  }
  virtual bool isconstant() const {
    return true;
  }
  Constant *clone() const {
    Constant *tmp=new Constant(_p[0]);
    return tmp;
  }
  double log_pdf(const double &x) const {
    assert(fabs(x-_p[0])<1E-8);
    return 0.0;
  }
  ostream& print(ostream &o) const {
    return o << _p[0];
  }
  double min() const {
    return _p[0];
  } 
  double max() const {
    return _p[0];
  }
  std::vector<double> parameters() const {
    return std::vector<double>(1,_p[0]);
  }
};
/** the normal distribution                                       */
class normal:
  public ctsdistribution {
public:
  normal(double mu, double sigma)
    :ctsdistribution(mu,sigma)  {
  }
  ~normal(){};
  double sample(rng &r) const {
    return _p[0]+r.normal()*_p[1];
  }
  vector<double> sample(rng &r,int n) const {
    std::vector<double> xx(n);
    for (int i=0;i<n;i++) xx[i]= _p[0]+r.normal()*_p[1];
    return xx;
  }
  void sample(rng &r,vector<double> &x) const {
    for (size_t i=0;i<x.size();i++) x[i]= _p[0]+r.normal()*_p[1];
  }
  normal *clone() const  {
    normal *tmp=new normal(mu(),sigma());
    return tmp;
  }
  double mu() const { return _p[0];}
  double sigma() const{ return _p[1];}
  double log_pdf(const double &x) const {
    // note that a specifically written version of this for
    // logarithms may be preferable
    return R::dnorm(x,_p[0],_p[1],1);

  }
  bool proper() const {
    return true;
  }
  ostream& print(ostream &o) const {
    return o << "normal(" << mu() <<"," <<sigma() <<")";
  }
  /// return the parameters as a vector
 std::vector<double> parameters() const {
    return std::vector<double>(1,_p[0]);
  }
};
/** the normal distribution                                       */
class lognormal:
  public ctsdistribution {
public:
  lognormal(double mu, double sigma)
    :ctsdistribution(mu,sigma)  {
  }
  ~lognormal(){};
  double sample(rng &r) const {
    return r.rlnorm(_p[0],_p[1]);
  }

  lognormal *clone() const  {
    return new lognormal(_p[0],_p[1]);
  }

  double log_pdf(const double &x) const {
    return R::dlnorm(x,_p[0],_p[1],1);

  }
  bool proper() const {
    return true;
  }
  ostream& print(ostream &o) const {
    return o << "lognormal(" << _p[0] <<"," << _p[1] <<")";
  }
  std::vector<double> parameters() const {
    return std::vector<double>(_p,_p+2);    
  }
};
/** The Gamma distribution
 * This is parametrised here in terns of shape and scale
 * The alternative parameterisation is shape and rate = 1/scale
 *
 * Hence the mean is a*b, and the variance a*b^2
 */
class Gamma
  : public ctsdistribution {
public:
  Gamma(double a, double s)
    :ctsdistribution(a,s) {
  }
  ~Gamma(){};
  double log_pdf(const double &x) const {
    return R::dgamma(x,_p[0],_p[1],1);

  }
  double sample(rng &r) const {
    return r.rGamma(_p[0],_p[1]);
  }
  vector<double> sample(rng &r,int n) const {
    vector<double> xx(n);
    for (int i=0;i<n;i++) xx[i]=r.rGamma(_p[0],_p[1]);
    return xx;
  }
  Gamma *clone() const {
    Gamma *tmp=new Gamma(_p[0],_p[1]);
    return tmp;
  }
  bool proper() const {
    return true;
  }
  double shape() const {
    return _p[0];
  }
  double a() const {
    return shape();
  }
  double scale() const {
    return _p[1];
  }
  double rate() const  {
    return 1./_p[1];
  }
  double mean() const {
    return shape()*scale();
  }
  double var() const {
    return shape()*scale()*scale();
  }
  ostream& print(ostream &o) const {
    return o << "Gamma(" << shape() <<"," <<scale() <<")";
  }
 std::vector<double> parameters() const {
    return std::vector<double>(_p,_p+2);
  }
private:
  Gamma();
};
/** The uniform distribution                                      */
class uniform
  : public ctsdistribution {
public:
  uniform(double a, double b)
    :ctsdistribution(a,b) {
  }
  ~uniform(){};
  double log_pdf(const double &x) const {
    if (x<_p[0])
      throw std::domain_error("x too low in log_pdf for uniform distribution");
    if (x>_p[1]) {
      std::ostringstream oss;
      oss << "x = " << x << " for distribution: ";
       print(oss);
       oss <<"x too high in log_pdf for uniform distribution";
        throw std::domain_error(oss.str().c_str());
    }
    return -log(_p[1]-_p[0]);
  }
  double sample(rng &r) const {
    return _p[0]+r()*(_p[1]-_p[0]);
  }
  uniform *clone() const {
    uniform *tmp=new uniform(_p[0],_p[1]);
    return tmp;
  }
  bool proper() const {
    return true;
  }
  double max() const {
    return _p[1];
  }
    double min() const {
    return _p[0];
  }
  ostream& print(ostream &o) const {
    return o << "uniform(" << min() <<"," <<max() <<")";
  }
 std::vector<double> parameters() const {
    return std::vector<double>(_p,_p+2);
  }
private:
  uniform();
};
/** The Beta distribution                                      */
class Beta
  : public ctsdistribution {
public:
  Beta(double a, double b)
    :ctsdistribution(a,b) {
  }
  ~Beta(){};
  double log_pdf(const double &x) const {
    return (a()-1.)*log(x)+(b()-1)*log(1.-x)
      -lgamma(a())-lgamma(b())+lgamma(a()+b());
  }
  double sample(rng &r) const {
    double x=r.rGamma(a(),1.0);
    double y=r.rGamma(b(),1.0);
    return x/(x+y);
  }
  Beta *clone() const {
    Beta *tmp=new Beta(_p[0],_p[1]);
    return tmp;
  }
  bool proper() const {
    return true;
  }
  double a()const {
    return _p[0];
  }
  double b() const {
    return _p[1];
  }
  double max() const {
    return 1.0;
  }
  double min() const {
    return 0.0;
  }
 double log_posterior(int x1, int x2) const {
    return -lgamma(x1+x2+a()+b())
      +lgamma(x2+b()) + lgamma(x1+a())
      -lgamma(a())- lgamma(b()) + lgamma(a()+b())
      + lgamma(x1+x2+1.)-lgamma(x1+1.) - lgamma(x2+1.);
  }
  ostream& print(ostream &o) const {
    return o << "Beta(" << a() <<"," << b() <<")";
  }
 std::vector<double> parameters() const {
   return std::vector<double>(_p,_p+2);
  }
private:
  Beta();
};
/** The Poisson distribution                                      */
class Poisson
  : public ctsdistribution {
public:
  Poisson(double mu)
    :ctsdistribution(mu) {
  }
  ~Poisson(){};
  double log_pdf(const double &x) const {
    return -lambda() + log(lambda())*x - lgamma(x+1.0);
  }
  double sample(rng &r) const {
    return r.rpoisson(lambda());
  }
  Poisson *clone() const  {
    Poisson *tmp=new Poisson(_p[0]);
    return tmp;
  }
  bool proper() const {
    return true;
  }
  double lambda() const {
    return _p[0];
  }
   ostream& print(ostream &o) const {
    return o << "Poisson(" << lambda() <<")";
  }
 std::vector<double> parameters() const {
    return std::vector<double>(1,_p[0]);
  }
private:
  Poisson();
};
/** The Geometric distribution - note x=1,2,3,...              */
class Geometric
  : public ctsdistribution {
public:
  Geometric(double p)
    :ctsdistribution(p) {
  }
  ~Geometric(){};
  double log_pdf(const double &x) const {
    return log(p()) + (x-1.)*log(1.-p());
  }
  double sample(rng &r) const {
    return r.rgeometric(p());
  }
  Geometric *clone()const  {
    Geometric *tmp=new Geometric(_p[0]);
    return tmp;
  }
  bool proper() const {
    return true;
  }
  double p()const {
    return _p[0];
  }
 ostream& print(ostream &o) const {
    return o << "Geometric(" << p()  <<")";
  }
 std::vector<double> parameters() const {
    return std::vector<double>(1,_p[0]);
  }
private:
  Geometric();
};

/**************************************************************************/
/*** Multivariate Distributions                                           */
/** Class for a multivariate continuous distribution on the real line     */
/**************************************************************************/
/** Abstract class for a multivariate continuous distribution on the real plane     */
class mvctsdistribution:
  public distribution<vector<double>,vector<double> > {
public:
 mvctsdistribution() { // for improper distribution
  }
  mvctsdistribution(const std::vector<double> &par) {
    _p=par;
  }
  virtual std::vector<double> sample(rng &r) const =0;

  ~mvctsdistribution() {
  }
  virtual double log_posterior(int x1, int x2) const {
    return 0.0;
  }
  virtual double log_pdf(const std::vector<double> &x) const =0;
  
 


};
/// multivariate constant distribution
class mvconstant :
 public mvctsdistribution {
public:
  mvconstant(const std::vector<double> &a)    {
      _p = a;
  }
  std::vector<double> sample(rng &r) const {
    return _p;
  }

  mvconstant *clone() const {
    mvconstant *tmp=new mvconstant(_p);
    return tmp;
  }

  double log_pdf(const std::vector<double> &x)  const {
    assert(x.size()==_p.size());
    return 0.0;
  }

  ostream &print(ostream &o) const {
    o << "mvconstant(";
    std::copy(_p.begin(),_p.end(),std::ostream_iterator<double>(o,","));
    return o << ")";
  }
};

class dirichlet
  : public mvctsdistribution {
public:
  dirichlet(const std::vector<double> &a)    {
      _p = a;
  }

  size_t size() const {
    return _p.size();
  }

  std::vector<double> sample(rng &r) const {
    return r.rdirichlet(_p);
  }

  dirichlet *clone() const {
    dirichlet *tmp=new dirichlet(_p);
    return tmp;
  }

  double log_pdf(const std::vector<double> &x)  const {
    // toremove
    assert(x.size()==_p.size());
  double lx=0.0,sum=0.0,temp=0.0;
    for (size_t i=0;i<x.size();i++) {
      lx += (_p[i]-1.0)*log(x[i]);
      sum += _p[i];
      temp += lgamma(_p[i]);
    }
    lx += lgamma(sum) - temp;
    return lx;

  }

  ostream &print(ostream &o) const {
    o << "dirichlet(";
    std::copy(_p.begin(),_p.end()-1,std::ostream_iterator<double>(o,","));
    return o <<  _p.back() << ")";
  }
};


ctsdistribution *ctsdistributionparse(std::istream &is);
mvctsdistribution *mvctsdistributionparse(std::istream &is);

ctsdistribution *fromstring(const std::string &s);
mvctsdistribution *fromstringmv(const std::string &s);

template <typename ITOR>
double ldnorm(ITOR start, ITOR lst, double mean, double sd) {
  ITOR curr=start;
  double temp=0.0;
  while (curr!=lst) {
#ifdef USR_R
    temp += pnorm(*(curr),mean,sd,1);
#else
    temp += log(gsl_ran_gaussian_pdf(*(curr)-mean,sd));
#endif
    curr++;
  }
  return temp;
}


/** The Normal Gamma distribution
 *
 */
class NormalGamma
  : public mvctsdistribution {
public:
  NormalGamma(const std::vector<double> &par)
    :mvctsdistribution(par) {
  }
 NormalGamma(double mu, double lambda, double alpha, double beta)
   :mvctsdistribution(std::vector<double>(4)) {
    _p[0]=mu;
    _p[1]=lambda;
    _p[2]=alpha;
    _p[3]=beta;
  }

  ~NormalGamma(){};
  double log_pdf(const std::vector<double> &x) const {
    return R::dgamma(x[1],alpha(),rate(),1)
      +R::dnorm(x[0],mu(),1./sqrt(lambda()*x[1]),1);
  }
  std::vector<double> sample(rng &r) const {
    std::vector<double> xx(2);
    xx[1] = r.rGamma(alpha(),rate());        // the precision
    xx[0] =  mu()+r.normal()* sqrt(1.0/(lambda()*xx[1]));
    return xx;
  }
 
  NormalGamma *clone() const {
    NormalGamma *tmp=new NormalGamma(_p[0],_p[1],_p[2],_p[3]);
    return tmp;
  }
  bool proper() const {
    return true;
  }
  double mu() const  {
    return _p[0];
  }
  double lambda() const {
    return _p[1];
  }
  double alpha() const {
    return _p[2];
  }
  double beta() const {
    return _p[3];
  }
  double rate() const {
    return 1./_p[3];
  }
  double scale() const {
    return _p[3];
  }
  double mean() const {
    return mu();
  }
  double var() const {
    return (lambda()+1.0)*beta()/(lambda()*(alpha()-1.));
  }
  ostream& print(ostream &o) const {
    return o << "NormalGamma(" << mu() << ","<<lambda()
	     << "," <<alpha() <<"," <<beta() <<")";
  }
  //std::vector<double> parameters() const {
  //  return std::vector<double>(&_p[0],&_p[4]);
  // }
private:
  NormalGamma();
};




double samplenewlogscale(const ctsdistribution *c,double &x, double tune, rng &r);
double samplenewlogscale(std::pair<double,ctsdistribution *> &val, double tune, rng &r);
#endif
