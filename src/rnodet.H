/** @file   */         
/*Time-stamp: <2012-05-25 00:40:40 nijw>  */ 
#ifndef __RNODET__H_IJW
#define __RNODET__H_IJW

#include <vector>
#include <cassert>
#include <iosfwd>
// utility headers
#include "utilityfunctionals.H"
//#include "util.H"
// Genetics Headers
#include "activesites.H"
#include "linesOfDescent.H"
#include "recombnode.H"
#include "recomblist.H"


namespace GenTL {

  ////////////////////////////////////////////////////////////////////
  /////RNODE//////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  /// tempated recombinant node
  template <class T>
  class rnode {
    friend class recomblist<T>;
    friend class lines_of_descent<T>;
  public:
    /** Constructor after a coalescent event                                            */
    rnode(const std::pair<int, int> &wh, const double t, rec_rates<T> &rr
	  ,std::vector<rnode<T> *> &root, const lines_of_descent<T> &lin)
      :left_(lin[wh.first]),right_(lin[wh.second]),time_(t) {
      l.push_back(new recombnode<T>(lin[wh.first],lin[wh.second],t,rr,root,this));
    };
    /** constructor for leaves.          */
    rnode(int loci):left_(0),right_(0),time_(0.0) {
      l.push_back(new recombnode<T>(loci,this));
    }; 
    /** a function to return the left and right descendents at position pos 
     * of there is none then throws an exception 
     * it returns 0 if left is a leaf                                            */
    rnode *left(position pos) {
      assert(left_->IsActive(pos));
      return left_->findnode(pos);
    }
    /** right descendent                                     */
    rnode *right(position pos) {
      assert(right_->IsActive(pos)); 
      return right_->findnode(pos);
    };
    
    /** return the recombnode after a recombination         */
    recombnode<T> *recombination(recombnode<T> *rn, const position where, const double time) {
      typename recomblist<T>::itor i; 
      for (i=l.begin();i!=l.end();i++) { 
	if (where < (*i)->last()) break;  	
      } 
      assert(i!=l.end());
      assert(*i==rn);
      // create the new recombination node
      // and insert it into the list before i;
      i=l.insert(i,new recombnode<T>(*i,where,time));   
      return *i;
    }
    /** return the recombnode after a gene conversion                                    */
    //todo not completed yet!
    recombnode<T> *geneconversion(recombnode<T> *rn, const position where,const int length
				  , const double time) {
      typename recomblist<T>::itor i,j; 
      for (i=l.begin();i!=l.end();i++) { 
	if (where < (*i)->last()) break;  	
      }
      for (j=j;j!=l.end();j++) { 
	if (where+length < (*j)->last()) break;  	
      } 
	
      assert(i!=l.end());
      assert(*i==rn);
      // create the new recombination node
      // and insert it into the list before i;
      i=l.insert(i,new recombnode<T>(*i,where,time));   
      return *i;
    }
    /** What is the first recombnode  - overloaded for constant and
     * non-constant recombnodes                                                           */
    recombnode<T> *first() {
      return l.front();
    }
    recombnode<T> *first() const {
      return l.front();
    }
    /** What is the last recombination node - overloaded for constant and
     * non-constant recombnodes                                                           */
    recombnode<T> *last()  const {
      return l.back();
    }
    recombnode<T> *last()  {
      return l.back();
    }
    /** extract the data (but not change it).                                             */
    const T &data() const {
      return data_;
    }; 
    /** extract the data                                                                  */
    T &data() { 
      return data_;
    }; 
    /** which positions are present                                                       */
    std::list<position> present() {
      std::list<position> lp;
      typename recomblist<T>::itor i = l.begin(); 
      while (i!=l.end()) {
	(*i++)->active.present(lp);
      }
      return lp;
    }
    /** what are the min and max active sites here                              */
    std::pair<position,position> minmaxactive() {
      std::pair<position,position> a;
      a.first=l.begin()->active.first();
      a.second=l.back().active.last();
      return a;
    }
    /** print the sites that are active for this line                                     */
    void print_active(std::ostream &o, const char *message) {
      o << message;
      l.print_active(o,"");
    }
    /** check that each of the recombnodes points to the right place                      */
    void check() {
      typename recomblist<T>::itor i=l.begin();
      while (i != l.end()) {
	assert((*i)->coal==this);
	i++;
      }
    }
    /** checkagain */
    bool checknode() {
      typedef std::list<position> lipos;
      lipos a=present();
     
      if (!isleaf()) {
	lipos::iterator ii=a.begin();
	while(ii!=a.end()) {
	  bool l=left_->IsActive(*ii);
	  bool r=right_->IsActive(*ii);
	  if (!(l||r)) {
        std::ostringstream oss;
	    oss << "problem with node in position " << *ii << std::endl;
	    oss << "left: ";
	    left_->print(oss,"");
	    oss << "right: ";
	    right_->print(oss,"");
	    this->print_active(oss,"");
        throw std::runtime_error(oss.str().c_str());
	  }
	  ii++;
	}
      }
      return true;
    }
    
    /** Is a position active ?                                                            */
    bool IsActive(position loc) {
      typename recomblist<T>::itor i=l.begin();
      while (i != l.end()) {
	if ((*i)->IsActive(loc)) return true;
	i++;
      }
      return false;
    }
    /** Is a recombination position present in this node            */
    bool IsRecombinationActive(position where) {
      if (IsActive(where)&&IsActive(where+1)) return true;
      return false;
    }
    /** length calculation                                                                */
    double length(position pos) {
      double tmp(0.0);
      rnode<T> *l=left(pos);
      if (l->isleaf()) {
	tmp+=time();
      } else {
	tmp += time()-l->time();
	tmp += l->length(pos);
      }
      rnode<T> *r=right(pos);
      if (r->isleaf()) { 
	tmp+=time();
      } else {
	tmp += time()-r->time();
	tmp += r->length(pos);
      }
      return tmp;
    }
    /** length calculation                                                                */
    int nleaves(position pos) {
      if (isleaf()) return 1;
      return left(pos)->nleaves(pos) 
	+ right(pos)->nleaves(pos);
    }
    /**           is this a leaf               */
    bool isleaf() {
      if (left_==0) return true;
      return false;
    }
    /** what are the positions of all recombinations before coalescences  */
    std::set<position> recombinations() {
      return l.recombpos();
    }  
    /** are any of of the samples that are numbered in wh present in this node                  
     * sample0 is the address of the first sample                         */
    //    note that the function below has not been tested
    bool present(const std::vector<size_t> &which, std::vector<rnode<T> *> samp) {
      if (isleaf()) {
	size_t num=std::distance(samp.begin(),find(samp.begin(),samp.end(),this));
	std::vector<size_t>::const_iterator it=std::find(which.begin(),which.end(),num);

	if (it!=which.end()) return true;
	return false;
      }
      bool lft,rght;
      if (left_->coal) 
	lft=left_->coal->present(which,samp);
      else lft=false;
      if (lft) return true;
      if (right_->coal) rght=right_->coal->present(which,samp);
      else rght=false;
      return rght;
      //    if (left_->coal->present(which,sample0)) return true;
      //      else return right_->coal->present(which,sample0); 
    }


    bool present(rnode<T> *which,position pos) {
      if (isleaf()) {
	if (which==this) return true;
	return false;
      }
      return (left(pos)->present(which,pos)||right(pos)->present(which,pos));
    }
  
    double time() const {
      return time_;
    }
    
    /** Public data              */
   recombnode<T> *left_,*right_;    //< the recombnodes that are the left and right descendants
   recomblist<T> l;  //< a list of recombnodes that happen through the recombination of the original node


  private:
    double time_;                     //< the time of the coalescence
    T data_;  
    /** Functions that are not defined           */
    rnode(const rnode &a);
  };


}
#endif
