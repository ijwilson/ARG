#ifndef RECOMBLIST_H__
#define RECOMBLIST_H__  

#include <list>
#include <set>

namespace GenTL {

//  const position maxloci=100000000;
  /** forward declarations                             */
  template<class T> class rec_rates;
  template <class T> class rnode;
  
  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////
  /** recomblist is a list of recombnodes - another helper class              */
  template<class T>
    class recomblist {
      friend class rnode<T>;
    private:
      std::list<recombnode<T> *> recombdata_;
    public:
      typedef typename std::list<recombnode<T> *>::iterator itor;
      /** constructor for this recomblist                                 */
      recomblist(){};
      recomblist(recombnode<T> *a) {
	recombdata_.push_back(a);
      }
      /** Destructor for this recomblist                                  */
      ~recomblist() {
	for_each(recombdata_.begin(),recombdata_.end(),
		 DeleteObject());
      } 
      /** iterator that points to the first recombinant node in the list     */
      itor begin() {
	return recombdata_.begin();
      }
      /** iterator that points to one past the last recombinant node in the list     */
      itor end() {
	return recombdata_.end();
      }
      void splice(itor pos,  recomblist &rlist) {
	//itor iif=find(data_.begin(),data_.end(),val);
	//assert (iif!=data_.end());
	recombdata_.splice(pos,rlist.recombdata_);
      }
      itor erase (itor pos) {
	return recombdata_.erase(pos);
      }
      recombnode<T> *front() {
	return recombdata_.front();
      }
      recombnode<T> *back() {
	return recombdata_.back();
      }
      void push_back(recombnode<T> * x) {
	return recombdata_.push_back(x);
      }
      itor insert(itor pos, recombnode<T>  *x ) {
	return recombdata_.insert(pos,x);
      }
      /** print the active sites in the list of sites        */
      std::ostream &print_active(std::ostream &o, const char *message) {
	if (recombdata_.empty()) {
	  o << "no active sites ";
	  return o;
	}
	itor i=recombdata_.begin();
	int count=0;
	while (i != recombdata_.end()) {
	  o << "address:" << *i << " time:" << (*i)->rtime 
	    << "  recombination " << (*i)->where << " stretch "<< ++count <<  " of "<< recombdata_.size() << ": "; 
	  (*i)->print(o,"");
	  i++;
	}
      return o;
      };
      /** returns the recombinations before a coalescence      */
      std::set<position> recombpos(){
	std::set<int> r;
	itor i=recombdata_.begin();
	while (i != recombdata_.end()) {
	  position ps=(*i)->where;
	  r.insert(ps);
	  i++; 
	}
	return r;
      }
      void fixlist(recombnode<T> *d) {
      // move all the nodes in the list to the list d->coal 
      // and remove the node pointed to by *d
      // what are the active sites in the descendent node 
      std::list<position> actived=d->active.activesites();
      // now loop through all the segments here
      assert(!recombdata_.empty());
      recombdata_.back()->where=d->where;
      for (itor jj=recombdata_.begin();jj!=recombdata_.end();) {
	// change the coalescent node to d->coal
	(*jj)->coal=d->coal;
	std::list<position> lp = (*jj)->active.activesites();
	std::list<position>::iterator li=lp.begin(),lf;
	while (li!=lp.end()) {
	  lf=find(actived.begin(),actived.end(),*li);
	  if (lf==actived.end()) {
	    (*jj)->active.erase(*li);
	  }
	  li++;
	}
	// erase the empty recombnodes (must be at the front or back)
  	if ((*jj)->active.empty()) { 
 	  delete *jj; 
 	  jj=recombdata_.erase(jj); 
 	} 
 	else jj++; 
      }     
    }
    void fixlist(recombnode<T> *d1, recombnode<T> *d2) {
      std::list<position> actived1=d1->active.activesites();
      std::list<position> actived2=d2->active.activesites();
      // now loop through all the segments here
      for (itor jj=recombdata_.begin();jj!=recombdata_.end();) {
	// change the coalescent node to d->coal
	std::list<position> lp = (*jj)->active.activesites();
	std::list<position>::iterator li=lp.begin(),lf1,lf2;
	while (li!=lp.end()) {
	  lf1=find(actived1.begin(),actived1.end(),*li);
	  lf2=find(actived2.begin(),actived2.end(),*li);
	  if (lf1==actived1.end()&&lf2==actived2.end()) {
	    (*jj)->active.erase(*li);
	  }
	  li++;
	}
	if ((*jj)->active.empty()) {
	  delete *jj;
	  jj=recombdata_.erase(jj);
	} else jj++;
      }
    }	
  };
}
#endif
