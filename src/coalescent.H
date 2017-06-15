/** @file  */
#ifndef COALESCENT_H_
#define COALESCENT_H_

#include "growthmodel.h"
#include "migmatrix.h"
#include "ijw_rand.h"
#include "gsl_sort.h"
#include "newio.h"

/** \brief a class for keeping track of coalescent times and event
 *
 * A class to keep track of coalesence times and events.  This shall be the 
 * base of a number of coalescent classes incorporating structure and growth
 */
namespace GenTL {
  template <class T> class tree;
  class coalescent {
  public:
    /** construct with n chromosomes  */
    coalescent(int a,growthmodel *gm,rng &r)
      :nleft_(a),lrand(r){
      if (gm==0) g=new constant_size();
      else g=gm->clone();
    };
    /** return the number of chromosomes left            */
    int nleft() const {
      return nleft_;
    };
    /** reset the number of chromosomes                  */
    void nleft(int nl) {
      nleft_=nl;
    }
    /// destructor - remember to remove growthmodel
    virtual ~coalescent() {delete g;};
    /** update after a coalescence  */
    virtual void UpdateCoalescence(std::pair<int,int> &lines, int remove) {
      assert(remove==1||remove==2);
      nleft_ -=remove;
    }
    /** simplified removal and addition */
    void operator++(){
      nleft_ -=1;
    };
    void operator--(){
      nleft_-=1;
    };
    /** update after a recombination  */
    virtual void UpdateRecombination(int recomb) {
      nleft_ +=1;
    }
    /** return the next time - for a simple coalescent this is always a coalescence */
    virtual double next(const double &time) {
      return g->nextcoal(nleft_,time,lrand);
    };
    /** which pair of lines coalesces */
    virtual std::pair<int, int> coalchromo() {
      std::pair <int,int> p=lrand.sample2intsorted(0,nleft_-1);
      return p;
    };
    virtual std::ostream &print(std::ostream &o){
      g->print(o);
      o << " " << nleft_ << " lines left" << std::endl;
      return o;
    };
    virtual std::ostream &printlocations(std::ostream &o){
        for (int i=0;i<nleft_;i++) o << "0 ";
      o << std::endl;				   
      return o;
    };

  protected:
    /// the number of chromosomes left
    int nleft_;
    /// the random number generator used
    rng &lrand;
    /// the growthmodel
    growthmodel *g;
  private:
    /// keep these out of the ways - get an error if called
    coalescent();
    coalescent(const coalescent &a);
    coalescent &operator=(const coalescent &a);
  };
  /*************************************************************/
  /**********   subdivcoal class         ***********************/
  /** Subdivided coalescent.
   * This class is pure virtual class representing structured 
   * coalescents - cannot be any concrete representations      */
  /*************************************************************/
  class subdivcoal: public coalescent {
  public:
    /** First the constructor                                                     */
    subdivcoal(const std::vector<int> &where, const std::vector<double> &relsize
	       , growthmodel *gm,rng &r):
      coalescent(where.size(),gm,r),npops_(relsize.size()),where_(where)
      ,left_(npops_),relsize_(relsize) {
      std::vector<int>::iterator i=where_.begin();
      while (i != where_.end()) left_[*i++] +=1;
    }
    /** destructor                                                     */
    ~subdivcoal(){
    };
    virtual std::ostream &printlocations(std::ostream &o){
      std::copy(where_.begin(),where_.end(),std::ostream_iterator<int>(o," "));
      o << std::endl;				   
      return o;
    }
    /** print information about the subdivided coalescent for debugging */
    virtual std::ostream &print(std::ostream &o);
    /** update after a coalescent                                     */
    void UpdateCoalescence(std::pair<int,int> &lines, int remove);
    /** update after a recombination   */
    void UpdateRecombination(int recomb);
    /** return the next time (may alter the internal structure  */
    double next(const double &time)=0;
    /** which lines coalesce?  */
    std::pair<int, int> coalchromo()=0;

  protected:  
/** Pick the coalescence event                     */
    std::pair<int,int> PickCoalescence();
    int whichcoal;
    size_t npops_;
    std::vector<int> where_;
    std::vector<int> left_;
    std::vector<double> relsize_;
  };
  /*************************************************************/
  /**********   structured coalescent class ********************/
  /*************************************************************/
  /** \brief a class for the simple structured coalescent
   * This takes an array of sample sizes and migration rates
   * A correct and useful implementation is going to take
   * a little but more thought (and perhaps design).
   */
  class structured_coalescent: public subdivcoal {
  public:
    structured_coalescent(const std::vector<int> &where
			  ,growthmodel *gm, mig_matrix &m
			  ,rng &r);
    /** migrate individual first to population second 
     *  Function to correct the rates and information after 
     * a migration event                                   
     * The first parameter gives the position of the line that
     * is moving and the second says where     */
    void migrate(const std::pair<int,int> &lineto);
    /** print information about the structured coalescent for debugging */
    std::ostream &print(std::ostream &o);
    /** update after a coalescent                                     */
    void UpdateCoalescence(std::pair<int,int> &lines, int remove);
    /** update after a recombination   */
    void UpdateRecombination(int recomb);
    /** return the next time, and if it is a coalescence  */
    double next(const double &time);
    /** which lines coalesce?  */
    std::pair<int, int> coalchromo();
  protected:
    /** Pick the migration event             */
    std::pair<int,int> PickMigration();
    mig_matrix &mm_;
    double mig_rate;
    std::vector<double> rates_out;
  };
  /*************************************************************/
  /**********   splitting  coalescent class ********************/
  /*************************************************************/
  /// Class for coalescen within splitting populations
  class splitting_coalescent: public subdivcoal {
  public:
    splitting_coalescent(const std::vector<int> &where
			 , growthmodel *gm, GenTL::tree<double> *ptree
			 ,rng &r);
    /** print information about the splitting coalescent for debugging  */
    std::ostream &print(std::ostream &o);
    /** return the next time, and if it is a coalescence  */
    double next(const double &time);
    /** which lines coalesce?  */
    std::pair<int, int> coalchromo();
  private:
    /** fix the tree after a splitting event           */
    void fixpoptree();
    /** Data                                           */
    int nextjoin;
    const GenTL::tree<double> *pt_;
    gsl_index indx;
  };

 /*************************************************************/
  /**********   structured coalescent class ********************/
  /*************************************************************/
  /** \brief a class for the simple structured coalescent
   * This takes an array of sample sizes and migration rates
   * A correct and useful implementation is going to take
   * a little but more thought (and perhaps design).
   */
  class splitting_structured_coalescent: public structured_coalescent {
  public:
    splitting_structured_coalescent(const std::vector<int> &where
				    , growthmodel *gm,  general_mig_matrix &m,
				    GenTL::tree<double> *ptree,rng &r);
    /** print information about the structured coalescent for debugging */
    std::ostream &print(std::ostream &o);
    /** update after a coalescent                                       */
    double next(const double &time);
    /** which lines coalesce?  */
    std::pair<int, int> coalchromo();
  private:
    /** fix the tree after a splitting event                            */
    void fixpoptree();
    /** Data                                                            */
    general_mig_matrix &gmm;
    int nextjoin;
    const GenTL::tree<double> *pt_;
    gsl_index indx;
  };

}


GenTL::coalescent *create_coalescent(int ss,const std::string &growth_model, \
const std::string &mig_model, const std::string &pop_tree,rng &rd);


#endif


