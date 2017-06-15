/** @file */
#ifndef RECOMBNIATION__HRATE__H__
#define RECOMBNIATION__HRATE__H__

#include <stdexcept>
#include <sstream>

#include "ijw_rand.H"
#include "utilityfunctionals.H"
#include "rnodet.H"

typedef std::deque<double> RATEARRAY;

namespace GenTL {
/** \brief A class to keep track of recombination rates
 *
 * This keeps track of recombination and works out the rates
 */
template <class T>
class rec_rates {
 public:
   /// the only constructor.
   rec_rates(int n, int loci, const std::vector<double> &recombination_rates, rng &r);

    double next(double t) {	
	assert(total_>=0.0||fabs(total_)<1E-9);
	if (total_<1E-9) return 1E10;
#ifdef CHECK
	std::cerr << "returning from rec_rates.next" << std::endl;
#endif 
	return 0.5*rand.sexp()/total_;  // note the change - added factor of 0.5
    }
    double rate(){return total_;}
    std::pair<int,int> recombination(lines_of_descent<T> &curr);

    void update_coal(const std::pair <int,int> &wh, const activetype &active, int removed);
    void update_recombination(const int where, lines_of_descent<T> &lin);
    int nleft(int i){return nleft_[i];}; /// active sites at each locus
    int loci_left(){return loci_left_;} /// number of active loci left
    void onelesslocus(){loci_left_-=1;}; /// 
    bool oneless(int locus){
	nleft_[locus] -=1;
	if (nleft_[locus]==1) return true;
	return false;
    };
 
    std::list<position> remove(const std::list<position> &a) {
      std::list<position> tmp;
      std::list<position>::const_iterator i=a.begin();
      while (i!=a.end()) {
	nleft_[*i] -=1;
	if (nleft_[*i]==1) {
	  tmp.push_back(*i);
	  loci_left_-=1;
	}
	i++;
      }
      return tmp;
    }
 
 ///  simple checks
    bool check(const lines_of_descent<T> &curr);
    void print_rates(std::ostream &o) {
	for (size_t i=0;i<cr.size();i++) o << cr[i] << " ";
	o << std::endl;
	o << "leftloc: ";
	for (size_t i=0;i<nleft_.size();i++) o << nleft_[i] <<" ";
	o << std::endl;
    }
  private:
    RATEARRAY cr; /// Recombination rates for current 
    std::vector<double> cumulative_rates; /// cumulative recombination rates
    std::vector<double> rates; ///  recombination rates  
    std::vector<int> nleft_; 
    int loci_left_;
    double total_;
    int n_;
    rng &rand;
// private functions
    double ind_recomb_rate(const activetype &active);
// not defined - private for safety
    rec_rates();
    rec_rates(const rec_rates &a);
    rec_rates &operator=(const rec_rates &rhs);
};

/** \brief Constructor for relative rates 
 *
 */  
template<class T>
rec_rates<T>::rec_rates(int n, int loci
			, const std::vector<double> &recombination_rates, rng &r):
  cumulative_rates(loci-1)
  ,rates(recombination_rates),nleft_(loci,n),loci_left_(loci),total_(0.0),n_(n),rand(r)
{
    /// now calculate the cumulative recombination rates. 
    std::transform(rates.begin(),rates.end()
	      ,cumulative_rates.begin(),cumsum<double>());
    total_=double(n)*cumulative_rates[loci-2];
    /// allocates cumulative_rates[loci-2] (overall rates) to cr.
    cr.assign(n,cumulative_rates[loci-2]);
}
//
/// the recombination rates for a single chromosome - used
/// for updating the recombination rate structures.
template<class T>
double rec_rates<T>::ind_recomb_rate(const activetype &active)
{
  assert(!active.empty());
    if (active.singleton()) return 0.0;

    // std::cout << "active size " << active.size() << std::endl;
    int first=active.first();//*(active.begin());
    int last=active.last();//*(--active.end());  
    assert(last>first);
    if (first==0) return cumulative_rates[last-1];
    else return cumulative_rates[last-1]-cumulative_rates[first-1];
}
//
//
/// update after a coalescence 
/// need to update the relative rates abd  
template<class T>
  void  rec_rates<T>::update_coal(const std::pair <int,int> &wh
				, const activetype &active, int removed)
{
    assert(n_=cr.size());
//  the lines at wh.first and wh.second coalesce.  
//  we move the last line to wh.second and change wh.first to the new line
//  get the new rate for the active sites at the coalescence position
    if (removed==2) {
	total_-= cr[wh.first]+cr[wh.second];
	if (fabs(total_)<1E-10) total_=0.0;
	cr[wh.first]=cr.back();
	cr.pop_back();
	cr[wh.second]=cr.back();
	cr.pop_back();
	n_-=2;
    } else {
	double dnew=ind_recomb_rate(active);
	total_+=dnew-cr[wh.first]-cr[wh.second];
	if (fabs(total_)<1E-10) total_=0.0;

/// put the cr values in their new positions
	cr[wh.first]=dnew;
	cr[wh.second]=cr.back();
	cr.pop_back();
	n_--;
    }
}
//
//
/// update the relativerates after recombination
template<class T>
void rec_rates<T>::update_recombination(const int where, lines_of_descent<T> &lin)
{
    //  std::cout << "lin size = "<<lin.size() << std::endl;
    //std::cout << "cr size = " << cr.size() << std::endl;
    double dnew1=ind_recomb_rate(lin[where]->active);
    double dnew2=ind_recomb_rate(lin.back()->active);

    total_+=dnew1+dnew2-cr[where];
    if (fabs(total_)<1E-10) total_=0.0;

    assert(dnew1+dnew2-cr[where]<0.0);
    
    cr[where]=dnew1;
    cr.push_back(dnew2);
    n_++;
}
//
//
/// return a pair with the number of the recombining line and
/// the position that recombines
template<class T>
  std::pair<int,int> rec_rates<T>::recombination(lines_of_descent<T> &curr)
{
  std::pair <int,int> p;
  assert(n_==int(curr.size()));    
  
  p.first=gen_from_p(cr,rand);  // the section that recombines (splits)
  assert(!curr[p.first]->active.empty());
  int a=curr[p.first]->active.first();
  int b = curr[p.first]->active.last();
  // std::cout << "between " << a <<" and " << b-1 << std::endl;
  assert(a<b);
  if (a-b==1) p.second=a;
  else p.second= gen_from_p(rates,a,b-1,rand); 
  assert(p.second<b);
  //  else p.second = gen_from_p(rates.begin()+a,rates.begin()+b-1,rand);
  return p;
}
//
//
//
template <class T>
bool rec_rates<T>::check(const lines_of_descent<T> &curr) {
  assert(total_>=0.0);
  double tmptotal=0.0;
  for (size_t i=0;i<cr.size();i++) {
	tmptotal+=cr[i];
	if (fabs(cr[i]-ind_recomb_rate(curr[i]->active)) > 1E-10) {
      std::ostringstream oss;
      oss << "error in " << i << std::endl;
      oss << "have cr= " <<cr[i]<< " with calculated value=" 
          << ind_recomb_rate(curr[i]->active)   << std::endl;
      curr[i]->print(oss,"");	
      throw std::runtime_error(oss.str().c_str());
	}
  }
  assert(fabs(total_-tmptotal)<1E-8);
  assert(n_=curr.size());
  return true;
}
}
#endif
