#ifndef RECOMBNODE_H__
#define RECOMBNODE_H__

#include <list>
#include <vector>

namespace GenTL {

  const position maxloci=100000000;
  template <class T> class rnode;
  template<class T> class rec_rates;

  //////////////////////////////////////////////////////////////////////
  /** \brief recombination node
   *
   * A class to build simple recombination trees
   */
  template<class T>
  class recombnode { 
    typedef recombnode* recptr;
  public:
    /// data
    position where;           //< after which position did the recombination occur 
    activetype active;        //< which positions are active 
    double rtime;             //< what was the time of the event that led to this recombnode
    rnode<T> *coal;           //< pointer to the coalescent node from which this derives.
    /// constructor for a leaf
    recombnode(int loci, rnode<T> *here)
      : where(maxloci),rtime(0.0),coal(here) {
      active.insertrange(0,loci);
    }
    // constructor after a recombination
    recombnode(recptr old, const int wh, double ti)
    : where(wh),rtime(ti),coal(old->coal) {
      assert(where>=old->first());
      assert(where< old->last());    
      active.split(wh,old->active);
    }
    // constructor after a gene conversion
    //    recombnode(recptr old, const int wh, double ti)
    //   : where(wh),rtime(ti),coal(old->coal) {
    /// assert(where>=old->first());
    //   assert(where< old->last());    
    //   active.split(wh,old->active);
    //  }
    /** constructor after a coalescence */
    recombnode(const recptr l,const recptr r, const double ti
	       ,rec_rates<T> &rr, std::vector<rnode<T> *> &root
	       , rnode<T> *here):where(maxloci), rtime(ti),coal(here)  {
      std::list<int> rts = rr.remove(active.combine(l->active,r->active));
      if (rts.size()>0) { 
	std::list<int>::iterator i=rts.begin();
	while (i!=rts.end()) {
	  root[*i]=here;
	  active.erase(*i);
	  i++;
	}
      }
    }
    /** the first active position                                             */
    position first() { 
      return active.first();
    }
    /** which positions are present in a recombnode?                          */
    std::list<position> present() {
      std::list<position> lp;
      active.present(lp);
      return lp;
    }
    /** the last active postion                                               */
    position last() {
	assert(!active.empty());
	return active.last();
    } 
    /** print the active sites at a recombnode                                */
    std::ostream &print(std::ostream &o, const char *message) {
      o << message;
      active.print(o);
      return o;
    };
    /** Is the position loc "active" at a particular recombnode               */ 
    bool IsActive(position loc) {
      return active.present(loc);
    }
    /** Is this recombnode a leaf node ?                                      */
    bool isleaf() const  {
      if (coal->left_==0) return true;
      return false;
    }
    /** findnode - find a node that is active at pos 
     * find a node that is either an internal node at pos or a leaf 
     * below the current recombnode                                           */ 
    rnode<T> *findnode(int pos) {
      recombnode<T> *tmp=this;
      for (;;) {
	assert(IsActive(pos));
	if (tmp->isleaf()) return tmp->coal;

	assert((tmp->coal->left_->IsActive(pos)||tmp->coal->right_->IsActive(pos)));

	if (tmp->coal->left_->IsActive(pos)) {
	  if (tmp->coal->right_->IsActive(pos))
	    return tmp->coal;
	  else tmp=tmp->coal->left_;
	} else {
	  tmp=tmp->coal->right_;
	}
      } 		
    }
  };
}

#endif
