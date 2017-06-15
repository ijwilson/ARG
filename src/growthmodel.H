/** @file */
#ifndef GROWTHMODEL_H__
#define GROWTHMODEL_H__

#include <utility>     // for pair
#include "ijw_rand.h"

using std::vector;

namespace GenTL {
  /** abstract base class                                                */
  class growthmodel {
  public:
    virtual double nextcoal(const int n,const double t,rng &r)=0;
    virtual double nextcoalrel(const int n,const double t,const double rel, rng &r)=0;
    /** what is the next coalescence for multiple populations */
    std::pair<double,int> nextcoalv(const vector<int> &n,double t,rng &r);
    /** what is the next coalescence for multiple populations with different population sizes */
    std::pair<double,int> nextcoalvrel(const vector<int> &n, const vector<double> rel
				       ,const double t,rng &r);
    virtual ~growthmodel(){};
    virtual growthmodel *clone()=0;
    virtual void print(std::ostream &o) const=0; 
  };
  /** constant population size                                           */
  class constant_size: public growthmodel {
  public:
    constant_size(){};
    ~constant_size(){};
    constant_size *clone()  {
      constant_size *gm=new constant_size();
      return gm;
    }	    
    double nextcoal(const int n,const double t,rng &r){
      return r.sexp()*2.0/(double(n*(n-1)));
    }
    double nextcoalrel(const int n,const double t,const double rel,rng &r){
      return rel*r.sexp()*2.0/(double(n*(n-1)));
    }
    void print(std::ostream &o) const {
      o << "constant" << std::endl; // with relative size "<<relsize_<<std::endl;
    };
  }; 
  /** class to generate coalescence times for an exponentially growing population  */
  class exponentialgrowth: public growthmodel {
  public:
    exponentialgrowth(double gr):
      growthrate_(gr) {
    };
    exponentialgrowth(const vector<double> &x)
      :growthrate_(x[0]){
    };

    double rate() const {
      return growthrate_;
    }
    
    ~exponentialgrowth(){};
    exponentialgrowth *clone()  {
      exponentialgrowth *gm=new exponentialgrowth(growthrate_);
      return gm;
    }	
    double nextcoal(const int n,const double last,rng &r){
      if (fabs(growthrate_)>1E-3) {
	return log(
		   exp(growthrate_*last)-log(r())*2.*growthrate_/((double)n*((double)n-1.0)) 
		   )/growthrate_-last;
      } else return r.sexp()*2.0/(double(n*(n-1)));
    }
    double nextcoalrel(const int n,const double last,const double rel, rng &r){
      if (fabs(growthrate_)>1E-3) {
	return log(
		   exp(growthrate_*last)-log(r())*2.*growthrate_*rel/((double)n*((double)n-1.0)) 
		   )/growthrate_-last;
      } else return rel*r.sexp()*2.0/(double(n*(n-1)));
    }

    void print(std::ostream &o) const {
      o << "exponential("<<growthrate_<<")" 
	<< std::endl;
    };
  private:
    double growthrate_;
  };
  /** exponentially growing from a base                                                */
  class expfrombase: public growthmodel {
  public:
    expfrombase(double gr,double t0)
      :t_0_(t0),growthrate_(gr) {
      if (t_0_<0.0) 
	throw std::domain_error("Error, time since growth starts must be > 0");
    };
    expfrombase(const vector<double> &x)
      :t_0_(x[1]),growthrate_(x[0]) {
    };
    ~expfrombase(){};
    expfrombase *clone()  {
      expfrombase *gm=new expfrombase(growthrate_,t_0_);
      return gm;
    }	
    double nextcoal(const int n,const double last,rng &r) {
      if (last>=t_0_) {
	return -log(r())*2.0/(double(n)*double(n-1));
      } else {
	double t = log(
		       exp(growthrate_*last)-log(r()*exp(growthrate_*t_0_)*2.*growthrate_
						 /(double(n)*double(n-1.0))))/growthrate_-last;
	if (t+last> t_0_) 
	  return (t_0_-last) -log(r())*2.0/(double(n)*double(n-1));
	else return t;
      }
    }	
    double nextcoalrel(const int n,const double last,const double rel,rng &r) {
      if (last>=t_0_) {
	return -rel*log(r())*2.0*rel/(double(n)*double(n-1));
      } else {
	double t = log(
		       exp(growthrate_*last)-log(r()*exp(growthrate_*t_0_)*2.*growthrate_*rel
						 /(double(n)*double(n-1.0))))/growthrate_-last;
	if (t+last> t_0_) 
	  return (t_0_-last) -log(r())*2.0*rel/(double(n)*double(n-1));
	else return t;
      }
    }		
    void print(std::ostream &o) const {
      o << "expfrombase(" 
	<<growthrate_ << "," << t_0_<<")"<<std::endl; 
    };
  private:
    double t_0_,growthrate_;
  };
  /** a simple three parameter bottleneck (t_0,a,t_1)                        */
  /*  Size 1 until time t_0 in the past, then changes to size a, for a further
   * time t_1, before returning to size 1                                   
   *       
   *
   *           |                                    |  
   *           |                                    |  
   *           |                                    |  
   *                            |   |
   *                            |   |    t1
   *                            | a |
   *           |                                    |  
   *           |                                    |  
   *           |                                    |  
   *    t_0    |                                    |  
   *           |                                    |  
   *           |                                    |                          */
  class bottleneck: public growthmodel {
  public:
    bottleneck(double t0,double a,double t1)
      :t_0_(t0),t_1_(t1),a_(a) {
    };
      bottleneck(const vector<double> &x)
      :t_0_(x[0]),t_1_(x[2]),a_(x[1]) {
    };
    ~bottleneck(){};
    growthmodel *clone()  {
      growthmodel *gm=new bottleneck(t_0_,a_,t_1_);
      return gm;
    }	
    double nextcoal(const int n,const double last,rng &r) {
      if (last >= t_0_+t_1_) {                             // back to the original size
	return -log(r())*2.0/(double(n)*double(n-1));
      } else if (last >= t_0_) {                            // in the bottleneck
	double t= -a_*log(r())*2.0/(double(n)*double(n-1));
	if (t>t_1_) return t_0_+t_1_-last+nextcoal(n,t_1_,r);
	else return t;
      } else {                                             // before the bottleneck
	double t= -log(r())*2.0/(double(n)*double(n-1));
	if (last+t>t_0_) return t_0_-last+ nextcoal(n,t_0_,r);
	else return t;
      }
    }
    double nextcoalrel(const int n,const double last,const double rel,rng &r) {
      if (last >= t_0_+t_1_) {                             // back to the original size
	return -log(r())*2.0*rel/(double(n)*double(n-1));
      } else if (last >= t_0_) {                            // in the bottleneck
	double t= -a_*log(r())*2.0*rel/(double(n)*double(n-1));
	if (t>t_1_) return t_0_+t_1_-last+nextcoalrel(n,t_1_,rel,r);
	else return t;
      } else {                                             // before the bottleneck
	double t= -log(r())*2.0*rel/(double(n)*double(n-1));
	if (last+t>t_0_) return t_0_-last+ nextcoalrel(n,t_0_,rel,r);
	else return t;
      }
    }	
    void print(std::ostream &o) const {
      o << "bottleneck(" 
	<< t_0_ << "," << a_<<"," << t_1_ <<")" << std::endl; 
    };
  private:
      double t_0_,t_1_,a_;
  };
  /** A general population size change model */
  /* starts with size a_0 and keeps this size for time t_0, 
   *     then changes to a_1 and stays this size for time t_1, 
   * and so on until it changes to size a_k at time t_{k-1}
   */
  class piecewise: public growthmodel {
  public:
    /** constructor from arrays (useful from R)        */
    piecewise(double *siz, double *time, int nsize) {
      for (int i=0;i<nsize-1;i++) {
        t.push_back(time[i]);
        a.push_back(siz[i]);
      }
      a.push_back(siz[nsize-1]);
    }
    /** constructor from a vector (for reading from a command line)  */
    piecewise(const vector<double> &x) {
      if (x.size()%2!=1) 
	throw std::domain_error("error, expected an odd length vector for piecewise");
      vector<double>::const_iterator i=x.begin();
      for(;;) {
        a.push_back(*i++);
        if (i==x.end()) break;
        t.push_back(*i++);
      }
    }
    ~piecewise(){};
    /** copy a pointer to the new piecewise object                     */
    piecewise *clone()  {
      piecewise *gm=new piecewise(*this);
      return gm;
    }	
    /** What is the additional time until the next coalescent event ?   */
    double nextcoal(const int n,const double last,rng &r) {
      // first find the first interval to start in
      size_t interval;
      for (interval=0;interval<t.size();interval++) 
	if (last<t[interval]) break;
      // if it is past the final time then just return
      if (interval==t.size()) 
	return -a.back()*log(r())*2.0/(double(n)*double(n-1));
      //how much longer
      double tim= -a[interval]*log(r())*2.0/(double(n)*double(n-1));
      // if this is past the time then use the memoryless property and restart
      if (last+tim>t[interval]) 
	return t[interval]-last+nextcoal(n,t[interval],r);
      else return tim;
    }
    double nextcoalrel(const int n,const double last,const double rel,rng &r) {
      // first find the first interval to start in
      size_t interval;
      for (interval=0;interval<t.size();interval++) 
	if (last<t[interval]) break;
      // if it is past the final time then just return
      if (interval==t.size()) 
	return -a.back()*log(r())*2.0*rel/(double(n)*double(n-1));
      //how much longer
      double tim= -a[interval]*log(r())*2.0*rel/(double(n)*double(n-1));
      // if this is past the time then use the memoryless property and restart
      if (last+tim>t[interval]) 
	return t[interval]-last+nextcoalrel(n,t[interval],rel,r);
      else return tim;
    } 
    void print(std::ostream &o) const {
      o << "piecewise(";
      for (unsigned int i=0;i<t.size();i++) {
	o << a[i] << "," << t[i] <<",";
      }
      o << a.back() <<")";
    };
  private:
    vector<double> t;
    vector<double> a;
  };
}
GenTL::growthmodel *gmread(const std::string &s);
std::ostream &operator<<(std::ostream &o, const GenTL::growthmodel &g);

#endif
