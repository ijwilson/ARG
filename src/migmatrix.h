/// @file          
// @file Time-stamp: <2012-05-30 17:10:25 nijw>            
#ifndef MIGMATRIX_H_
#define MIGMATRIX_H_ 

#include <utility> // for pair
#include "ijw_rand.h"

namespace GenTL {
  /** \brief An abstract class to hold "migration matrices"
   *
   *  This class holds the generator matrix for the rates
   *   of migration between "island" populations        */
  class mig_matrix {
  public:
    /** virtual dtor                        */
    virtual ~mig_matrix() {
    };
    /** the rate of going from->to          */    
    virtual double operator()(int from, int to) const  =0;
    /** a "virtual" constructor             */
    virtual mig_matrix *clone() const =0;
    /** print                               */
    virtual std::ostream &print(std::ostream &o) const =0;
    virtual int npops() const =0;
    /** The rate of leaving population from */
    virtual double RateOut(int from) const =0;
    /** default constructor                 */
    mig_matrix() {
    };
    virtual int sample(int from, rng &r) const =0;
  };
  /** \brief an island migration matrix.
   * The rate of going from island i to j is 1/(k-1)
   * where k is the number of islands
   */
  class island: public mig_matrix {
  public:
    island(size_t k, double m):rate_(m),k_(k){};
    
    ~island() {
    };
    
    island *clone() const {
      island *tmp = new island(k_,rate_);
      return tmp;
    };
    
    double operator()(int from, int to) const {
      assert(from<k_ && from >=0);
      assert(to < k_ && to >= 0);
      return rate_/double(k_ -1);
    };
    
    std::ostream &print(std::ostream &o) const {
      o << "Island("<< k_ << "," << rate_ << ")";
      return o;
    };
    int npops() const {return k_;};
    //
    double RateOut(int from) const {
      return rate_;
    };
    /** sample a new location for the individual */
    int sample(int from,rng &r) const {
      int a=r.rint(k_-1);
      if (a>=from) a++;
      return a;
    }
    
  private:
    /// the migration rate
    double rate_;
    /// the number of islands
    int k_;
  };
  /** \brief A one dimensional stepping stone model
   * A simple 1-d stepping stone model.  This is a simpler
   * way to set this up than through a migration matrix
   */ 
  class stepping_stone1D: public mig_matrix {
  public:
    ~stepping_stone1D(){};
    double operator()(size_t from, size_t to) const {
      assert(from>=0 && from < length_);
      assert(to >= 0 && to <= length_);
      if (abs(int(from)-int(to))==1) return rate_/2.;
      else return 0.0;
    };
    std::ostream &print(std::ostream &o) const {
      o << "Steppingstone1D("<<length_<<"," << rate_ << ")";
      return o;
    }
    //
    int npops() const  {return length_;}
    //
    double RateOut(int from) const  {
      return rate_;
    };
    int sample(unsigned int from,rng &r) const {
      if (from==0) return 1;
      if (from==length_-1) return length_-2;
      if (r()<0.5) return from-1;
      else return from+1;
    }
  private:
    size_t length_;
    double rate_;
  };
  /** \brief A two dimensional stepping stone model
   * A two dimensional stepping stone model with relative 
   * dispersal rates of m_northsouth and m_eastwest */
  class stepping_stone2D: public mig_matrix {
  public:
    /** are sites North South, East West or Distant Neighbours */
    enum proximity {
      NS,EW,DN
    };
    /** constructor - note the pairs are ew, ns */
    stepping_stone2D(std::pair<int,int> dim, std::pair<double,double> rate):
      ns_(dim.second),ew_(dim.first)
      ,rate_ns_(rate.first),rate_ew_(rate.second) {};

    ~stepping_stone2D(){};

    stepping_stone2D *clone() const {
      std::pair<int,int> n(ew_,ns_);
      std::pair<double,double> rat(rate_ew_,rate_ns_);
      stepping_stone2D *s=new  stepping_stone2D(n,rat);
      return s;
    }
  
    double operator()(int from, int to) const {
      assert(from>=0 && from < ns_*ew_);
      assert(to>=- 0 && to < ns_*ew_);
      proximity p=neighbour(from,to);
      if (p==EW) return rate_ew_;
      else if (p==NS) return rate_ns_;
      else return 0.0;
    };

    std::ostream &print(std::ostream &o) const {
      o << "Steppingstone2D({" << ew_ <<"," << ns_ 
	<< "},{" << rate_ew_ << "," << rate_ns_ <<"})";
      return o;
    };

    // npops
    int npops() const {
      return ns_*ew_;
    }
    //
    double RateOut(int from) const  {
      return rate_ns_+rate_ew_;
    };
    /** sample where it goes to              */
    int sample(int from,rng &r) const {
      for (;;) { // do something really simple !!!???!!!
	int ret;
	if (r()<rate_ns_/(rate_ns_+rate_ew_)) {
	  if (r()<0.5) ret=from-ns_;
	  else ret=from+ns_;
	} else {
	  if (r()<0.5) ret=from+1;
	  else ret=from-1;
	}
	proximity p=neighbour(from,ret);
	if (p != DN) return ret;
      }
    }
    
  private:
    /** the positions are labelled from 0 (top left) to ns_*ew_-1 (bottom right)
     * by row, i.e
     *  0  1  2  3  4  5  6  7
     *  8  9 10 11 12 13 14 15
     * 16 17 18 19 20 21 22 23
     * for a 3 by 8 grid */
    /** are two sites on the same column? */
    bool incol(int a, int b) const {
      if (a%ew_==b%ew_) return true;
      return false;
    };
    /** Are two sites in the same row? */
    bool inrow(int a, int b) const {
      if (a/ew_==b/ew_) return true;
      return false;
    };
    /** what is the proximity of a pair of sites */
    proximity neighbour(int from_a, int to_b) const {
      if (inrow(from_a,to_b)) {
	if (abs(from_a-to_b)==1) return EW;
	else return DN;
      }
      if (incol(from_a,to_b)) {
	if (abs(from_a-to_b)==ew_) return NS;
      }
      return DN;
    }
    int ns_,ew_;
    double rate_ns_, rate_ew_;
  };


  /** \brief A class to hold the "migration" generator
   *  This class holds the generator matrix for the rates
   * of migration between "island" populations        */
  class general_mig_matrix: public mig_matrix {
  public:
    /// constructor for a symetric matrix
    general_mig_matrix(unsigned int siz,  std::vector<double> &m )
      // use the Meyes Effective STL trick to "convert" vector to array
      :m_(siz) {//,siz,&m[0]) {
      int start=0;
      for (size_t i=0;i<siz;i++) {
	m_[i].resize(siz);
	for (size_t j=0;j<siz;j++,start++) {
	  m_[i][j]=m[start];
	}
      }
      assert(siz*siz==m.size());
    }; 
    /// constructor from another migmatrix
    general_mig_matrix(const mig_matrix &mm);
 
    general_mig_matrix(const std::vector<std::vector<double> > &a) 
      //  :m_(a.copy()) {
      {
	m_=a;
      }; 
    // simple dtor
    ~general_mig_matrix(){};
    double operator()(int from, int to) const {
      return m_[from][to];
      }
    // print
    std::ostream &print(std::ostream &o) const ;
    // clone
    general_mig_matrix *clone() const {
      general_mig_matrix *tmp = new general_mig_matrix(m_);
      return tmp;
    };
    // npops
    int npops() const {
      //return m_.dim1();
      return m_.size();
    }
    //
    double RateOut(int from) const {
      // return std::accumulate(m_[from],m_[from]+m_.dim1(),0.0);
      return std::accumulate(m_[from].begin(),m_[from].end(),0.0);
    };
  //
    int sample(int from, rng &r) const {
      // return gen_from_pb(m_[from],m_.dim1(),r);
      return gen_from_p(m_[from],r);
    }
    /** remove a population      */
    void remove_pop(int remove) {
      for (size_t i=0;i<m_.size();i++) {
	m_[remove][i]=0.0;
	m_[i][remove]=0.0;
      }
    }
  private:
    //   const TNT::Array2D<double>  m_;
    std::vector<std::vector<double> > m_;
    general_mig_matrix();
    general_mig_matrix &operator=(const general_mig_matrix &a);
    general_mig_matrix(const general_mig_matrix &a);  
  };

}

GenTL::mig_matrix *mmread(const std::string &s);
std::ostream &operator<<(std::ostream &o, const GenTL::mig_matrix &m);
#endif
