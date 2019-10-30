/** @file  
 * Ancestral Recombination Graph       */
// Time-stamp: <2012-06-08 16:04:27 nijw>
#ifndef RTREE__H_IJW
#define RTREE__H_IJW

#include <stdexcept>
//#include "tnt/tnt.h"
// utilities
#include "utilityfunctionals.h" // for DeleteObject
// Genetics Headers
#include "gentldeclr.h"
#include "recombinationrates.h"
#include "rnodet.h"
#include "nodet.h"
#include "coalescent.h"

namespace GenTL {
/** The default collector class - just counts the recombinations  *
 * Other collector classes can easily be written - they just have 
 * have the operator()(time, event,which) method which is used within 
 * rtree::simulate                                                    */
struct recombinations {
  /** default constructor                                             */
  recombinations():rec_(0){};
  /** function overload operator - needed for any COLLECTOR structure */
  void operator()(double time, EventType e,int which=0) {
    if (e==Recombination) rec_ +=1;
  };
  /** function operator that just returns the result                  */
  long operator()() {
    return result();
  };
  /** Return the number of recombinations                             */
  long result() {
    return rec_;
  };
  /** Data -- keeps the number of recombinations                      */
  long rec_;
};
struct recombinationLocations {
  /** default constructor                                             */
  recombinationLocations() {};
  /** function overload operator - needed for any COLLECTOR structure */
  void operator()(double time, EventType e,int which=0) {
    if (e==Recombination) where_.push_back(which);
  };
  /** function operator that just returns the result                  */
  std::list<int> &operator()() {
    return where_;
  };
  /** Return the list of recombination locations                      */
  std::list<int> &result() {
    return where_;
  };
  /** return count */
  int count() { 
    return where_.size();
  }
  /** Data -- keeps the number of recombinations                      */
  std::list<int> where_;
};

/** ARG (ancestral recombination graph) class. 
 *  The general class for 
 simulating recombination trees                                    */
template<class T,typename COLLECTOR> 
class ARG {
  typename std::vector<rnode<T> *>::iterator samp_itor;
public:
  /** constructor which creates its own coalescent object             */
  ARG(int ss, int nl,const std::vector<double> &recomb_rates,growthmodel *g,rng &r,COLLECTOR &clt)
    :sample(ss)
  ,root(nl)
  ,collect(clt) {
    for (int i=0;i<ss;i++) sample[i]=new rnode<T>(nl);
    coalescent cc(ss,g,r);
    simulate(cc,recomb_rates,r);
  };
  /** Constructor which requires a coalescent object (which includes the sample size          */
  ARG(int nl,const std::vector<double> &recomb_rates, coalescent &c,rng &r,COLLECTOR &clt)
    :root(nl), sample(c.nleft()) ,collect(clt) {
    for (int i=0;i<c.nleft();i++) sample[i]=new rnode<T>(nl);
    simulate(c,recomb_rates,r);
  };
  /** dtor                 
   * This particular formulation is that of Meyers - effective STL                             */
  ~ARG() {
  for_each(ancestors.begin(),ancestors.end(),
           DeleteObject());
    for_each(sample.begin(),sample.end(),
             DeleteObject());
    Rprintf("runing destructor on ARG\n");
  }
  /** prune the ancestral recombination graph  */
  void prune(const std::vector<size_t> &wh) {     
    // which recombnodes (ancestors to the nodes listed in wh)  are still in use?  
    std::list<recombnode<T> *> used;
    for (size_t i=0;i<wh.size();i++) {
      if (wh[i]<0 or wh[i]>=sample.size()) 
        throw::std::domain_error("Problem with pruning values");
      used.insert(used.end(),sample[wh[i]]->l.begin(),sample[wh[i]]->l.end());
    }
    typename std::list<rnode<T>* >::iterator ii=ancestors.begin();
    while (ii!=ancestors.end()) {
      
      typename recomblist<T>::itor  
      left=std::find(used.begin(),used.end(),(*ii)->left_),
        right=std::find(used.begin(),used.end(),(*ii)->right_);
      
      if (left!=used.end()) {
        if (right!=used.end()) {
          // both present 
          (*ii)->l.fixlist(*left,*right); 
          used.insert(left,(*ii)->l.begin(),(*ii)->l.end());
          used.erase(left);
          used.erase(right);
          ii++;
        } else {
          // left present but not right
          rnode<T> *newd=(*ii)->left_->coal;
          (*ii)->l.fixlist(*left); // remove all those sections that are not needed 
          used.insert(used.end(),(*ii)->l.begin(),(*ii)->l.end());
          used.erase(left);
          
          typename recomblist<T>::itor iif=find(newd->l.begin(),newd->l.end(),(*ii)->left_);
          
          assert (iif!=newd->l.end());
          newd->l.splice(iif,(*ii)->l);
          newd->l.erase(iif);
          delete(*ii);
          ii=ancestors.erase(ii);
        }
      } else if (right!=used.end()) {
        // right present but not left
        rnode<T> *newd=(*ii)->right_->coal;
        
        (*ii)->l.fixlist(*right);  // remove all the recombnodes that are not needed
        used.insert(used.end(),(*ii)->l.begin(),(*ii)->l.end());
        used.erase(right);
        
        typename recomblist<T>::itor iif=find(newd->l.begin(),newd->l.end(),(*ii)->right_);
        //	 position recombpos=(*ii)->right_.recombination;
        assert (iif!=newd->l.end());
        newd->l.splice(iif,(*ii)->l);
        newd->l.erase(iif);
        delete(*ii);
        ii=ancestors.erase(ii);
      } else { // neither present
        delete(*ii);
        ii=ancestors.erase(ii);
      }
    }
    remakeroots(wh);
    checkroots(wh);
    // Now remove the samples and delete those nodes that are no longer needed
    std::vector<rnode<T> *> newsample;
    for (size_t i=0;i<wh.size();i++) {
      newsample.push_back(sample[wh[i]]);
      sample[wh[i]] = 0;// make sure that we don't erase the samples that should be kept
    }
    for_each(sample.begin(),sample.end(), DeleteObject());
    sample=newsample;
  }
  
  /** remake the roots after a pruning operation                     */
  void remakeroots(const std::vector<size_t> &wh) {
    for (size_t pos=0;pos<root.size();pos++) {
      typename std::list<rnode<T>* >::reverse_iterator ri=ancestors.rbegin();
      for (;;)  {
        if ((*ri)->left_->IsActive(pos)&&(*ri)->right_->IsActive(pos)) {
          root[pos]=*ri;
          break;
        }
        ri++;
        if (ri==ancestors.rend()) { 
          throw std::runtime_error("should enver get here in remakeroots\n");
        }
      }
    }
  }
 /* TNT::Array2D<position> sharedsections(const std::vector<int> &wh, position pos)
  {
    TNT::Array2D<position> minmax(wh.size(),2);
    // extract the tree at a position
    node<rnode<T> *> *myroot=extract_node(root[pos],pos,sample);
    
    for (size_t i=0;i<wh.size();i++) {
      // find the node at the correct depth for wh (wasteful...)
      
      node<rnode<T> *> *current=myroot->find_node(wh[i]);
      if (current==0) {
        std::ostringstream oss;
        oss << "label " << wh[i] <<" not found";
        throw std::runtime_error(oss.str().c_str());
      }
      node<rnode<T> *> *anc=current->up;
      
      position lf=anc->data()->left_->active.first();
      position lm=anc->data()->left_->active.last();
      position rf=anc->data()->right_->active.first();
      position rm=anc->data()->right_->active.last();
      
      minmax[i][0]=std::max(lf,rf);
      minmax[i][1]=std::min(lm,rm);
      
    }
    recursively_destroy_node(myroot);
    return minmax;
  }
  TNT::Array2D<position> sharedsectionK(const std::vector<int> &wh, position pos, int depth)
  {
    TNT::Array2D<position> minmax(wh.size(),2);
    // extract the tree at a position
    node<rnode<T> *> *myroot=extract_node(root[pos],pos,sample);
    
    for (size_t i=0;i<wh.size();i++) {
      // find node with label wh[i]
      node<rnode<T> *> *current=myroot->find_node(wh[i]);
      
      if (current==0) {
        std::ostringstream oss;
        oss << "label " << wh[i] <<" not found";
        throw std::runtime_error(oss.str().c_str());
      }
      // now go up until you have sufficient ancestors
      // this needs to be fixed sometime to account for the fact that you 
      // may not be looking at all the noeds on a tree and only want those in wh!
      node<rnode<T> *> *anc=current->up;
      //	while (anc->count_leaves()<depth) anc=anc->up;
      //	if (false) {
      for(;;) {
        std::vector<int> leaves;
        anc->leaf_labels(leaves);
        std::vector<int>::iterator ii=leaves.begin();
        int count=0;
        while (ii!=leaves.end()) {
          if (std::find(wh.begin(),wh.end(),*ii)!=wh.end()) count++;
          ii++;
        }
        if (count>=depth) break;
        else anc=anc->up;
      }
      
      typedef std::vector<node <rnode<T> *> *> nodevec;
      
      nodevec a;
      anc->get_internal(a);
      std::vector<position> l,r;
      typename nodevec::iterator iii=a.begin();
      while (iii!=a.end()) {
        l.push_back((*iii)->data()->left_->active.first());
        r.push_back((*iii)->data()->left_->active.last());
        l.push_back((*iii)->data()->right_->active.first());
        r.push_back((*iii)->data()->right_->active.last());
        iii++;
      }
      minmax[i][0]=*(std::max_element(l.begin(),l.end()));
      minmax[i][1]=*(std::min_element(r.begin(),r.end()));
      
    }
    recursively_destroy_node(myroot);
    return minmax;
  }
  */
  /** What is the time of the root at position loc                                             */
  double roottime(position loc) {
    return root[loc]->time();
  };
  /** What is the length of the tree at position loc                                           */
  double length(position loc) {
    return root[loc]->length(loc);
  } 
  /** How many coalescent nodes are there in this tree                                         */
  long CoalescentNodes() const { 
    return ancestors.size();
  };
  /** is this node a root?                                              */
  bool isroot(rnode<T> *any) {
    if (find(root.begin(),root.end(),any)!=root.end()) return true;
    return false;
  }
  /** root for which positions?                                        */
  std::vector<int> rootwhichpositions(rnode<T> *n) {
    std::vector<int> pos;
    for (size_t i=0;i<root.size();i++) {
      if (root[i]==n) pos.push_back(i);
    }
    return pos;
  }
  //    T data(int ii) {
  //     return t.sample[ii]->data();
  //  }
  /** randomise the node ordering                                                              */
  void randomise(rng &ran) {
    typename std::list<rnode<T> *>::iterator i=ancestors.begin();
    while (i != ancestors.end()) {
      if (ran()<0.5) {
        recombnode<T> *tmp=(*i)->left_;
        (*i)->left_=(*i)->right_;
        (*i)->right_=tmp;
      } 
      i++;
    }
  }
  /** check that the roots have descendents to the left and right at each position             */
  bool checkroots(const std::vector<size_t> &wh) {
    for (size_t i=0;i<root.size();i++) {
      assert((root[i]->left_->IsActive(i)&&root[i]->right_->IsActive(i)));
      bool left=root[i]->left_->coal->present(wh,sample);
      bool right=root[i]->right_->coal->present(wh,sample);
      if (!(left&&right)) {
        std::ostringstream oss;
        oss << " problem with root of site " << i << std::endl;
        oss << left << "-" << right  << std::endl; 
        throw std::runtime_error(oss.str().c_str());
      }
    }
    return true;
  }
  /** What is the TMRCA for a pair of individuals s1 and s2 at position pos                   */
  double TMRCA(int s1, int s2, position pos) {
    rnode<T> *samp1=sample[s1];
    rnode<T> *samp2=sample[s2];
    rnode<T> *curr=root[pos];
    while (!curr->isleaf()) {
      if (curr->left(pos)->present(samp1,pos)) {
        if (curr->left(pos)->present(samp2,pos)) {
          curr=curr->left(pos);
        } else {
          if (curr->right(pos)->present(samp2,pos)) return curr->time();
          else  throw std::runtime_error("Error A in TMRCA\n");
        }
      } else {
        if (!curr->right(pos)->present(samp1,pos)) 
          throw std::runtime_error("Error B in TMRCA\n");
        if (curr->left(pos)->present(samp2,pos)) {
          return curr->time();
        } else {
          if (curr->right(pos)->present(samp2,pos))  curr=curr->right(pos);
          else  throw std::runtime_error("Error C in TMRCA\n");
        }
      }
    }
    throw std::runtime_error("Should never get here in TMRCA\n");
    return -99.99;
  }
  /** return a list of rnodet pointers that give the path to the sample samp
   //  //     * which should be the number of one of the leaves.  The leaf itself is not returned */
   //     std::list<rnode<T> *> path_to_node(rnode<T> *samp, position pos) {
   //       rnode<T> *curr=root(pos);
   //       std::list<rnode<T> *> ll;
   //       while (!(*curr)->isleaf()) {
   // 	if (curr->present(samp,pos)) curr=curr->left(pos);
   // 	else curr=curr->right(pos);
   // 	ll.push_back(curr);
   //       }
   //       return ll;
   //     }  
   /** Simulate the tree                                                                       */
   void simulate(coalescent &c, const std::vector<double> recomb_rates,rng &r);
  void simulate(coalescent &c, const std::vector<double> recomb_rates
                  , double f, double lambda, rng &r);
  
  std::vector<rnode<T>* > root;    //< a vector of roots
  std::vector<rnode<T> *> sample;   //< the sampled leaves
  std::list<rnode<T>* > ancestors; //< ancestral coalescences
  /** How many loci?                                                                          */
  int nloc() const  {
    return sample[0]->last()->active.last()+1;
  }
  int sites() const  {
    return root.size();
  }
  /** get the data   */
  T &data(int i) {
    return sample[i]->data();
  }
private:
  COLLECTOR &collect;
  ARG();
};

/** A utility function to do the simulations for all constructors                             */
template <class T,typename COLLECTOR>
void ARG<T,COLLECTOR>::simulate(coalescent &c, const std::vector<double> recomb_rates,rng &r)
{
#ifdef CHECK
  std::cerr << "Before simulation have " << c.nleft() << " lines" << std::endl;
#endif
  double t=0.0;
  lines_of_descent<T> lines(sample);
  rec_rates<T> rr(c.nleft(),recomb_rates.size()+1,recomb_rates,r);
  for (;;) {
#ifdef CHECK
    std::cerr << "================================================" << std::endl 
              << "Time is now " << t << std::endl;
    lines.print(std::cout);
    std::cerr << "total rate "<<rr.rate()<<std::endl;
    rr.print_rates(std::cerr);
    rr.check(lines);
    c.print(std::cerr);
#endif    
    double NextRecombination=rr.next(t);
    double NextEvent=c.next(t);
#ifdef CHECK
    std::cerr << "NextRecombination = " << NextRecombination 
              << "Next Event " << NextEvent << std::endl;
#endif
    if (NextEvent<NextRecombination) {
      t += NextEvent;
      std::pair<int,int> chromo = c.coalchromo();
#ifdef CHECK
      std::cout << "coalescing "<< chromo.first 
                << " and " << chromo.second << std::endl;
#endif
      if (chromo.first!=chromo.second) { // a coalescence 
        collect(t,Coalescence);
        ancestors.push_back(new rnode<T>(chromo,t,rr,root,lines));
        // Now fix lines 
        int removed=lines.coal(ancestors.back(),chromo);
        c.UpdateCoalescence(chromo,removed);    
        rr.update_coal(chromo,lines[chromo.first]->active,removed); 
      } else collect(t,Migration);  // other
    } else {//	  recombination
      t += NextRecombination;
      std::pair <int,int> wh=rr.recombination(lines);
      collect(t,Recombination,wh.second);
#ifdef CHECK
      std::cout << "recombination for line " << wh.first
                << " at position " << wh.second << std::endl;
#endif
      recombnode <T> *nw=
      lines[wh.first]->coal->recombination(lines[wh.first],wh.second,t);	    
      
      lines.push_back(nw);
      c.UpdateRecombination(wh.first);
      rr.update_recombination(wh.first,lines);
    }
    if (c.nleft()<=1) return;
  }
}
/** A utility function to do the simulations for all constructors                             
 * This time with gene conversion rather than just crossing over 
 * here f is the ratio of gene conversion to crossing over and
 * lambda is the mean tract length - where the tract length has a 
 * geometric distribution                                          */
//todo   not completed yet!
#ifdef FALSE
template <class T,typename COLLECTOR>
void ARG<T,COLLECTOR>::simulate(coalescent &c, const std::vector<double> recomb_rates
                                  , double f, double lambda, rng &r)
{
#ifdef CHECK
  std::cerr << "Before simulation have " << c.nleft() << " lines" << std::endl;
#endif
  double t=0.0;
  lines_of_descent<T> lines(sample);
  rec_rates<T> rr(c.nleft(),recomb_rates.size()+1,recomb_rates,r);
  for (;;) {
#ifdef CHECK
    std::cerr << "================================================" << std::endl 
              << "Time is now " << t << std::endl;
    lines.print(std::cout);
    std::cerr << "total rate "<<rr.rate()<<std::endl;
    rr.print_rates(std::cerr);
    rr.check(lines);
    c.print(std::cerr);
#endif    
    double NextRecombination=rr.next(t);
    double NextEvent=c.next(t);
    
#ifdef CHECK
    std::cerr << "NextRecombination = " << NextRecombination 
              << "Next Event " << NextEvent << std::endl;
#endif
    if (NextEvent<NextRecombination) {
      t += NextEvent;
      std::pair<int,int> chromo = c.coalchromo();
#ifdef CHECK
      std::cout << "coalescing "<< chromo.first 
                << " and " << chromo.second << std::endl;
#endif
      if (chromo.first!=chromo.second) { // a coalescence 
        collect(t,Coalescence);
        ancestors.push_back(new rnode<T>(chromo,t,rr,root,lines));
        // Now fix lines 
        int removed=lines.coal(ancestors.back(),chromo);
        c.UpdateCoalescence(chromo,removed);    
        rr.update_coal(chromo,lines[chromo.first]->active,removed); 
      } else collect(t,Migration);  // other
    } else {//	  recombination
      t += NextRecombination;
      // which node recombines and at what position
      std::pair <int,int> wh=rr.recombination(lines);
      
      bool GeneConversion=(r()<f);
      int tracklength;
      if (GeneConversion) {
        int trackLength = r.rgeometric(1./lambda);
        // check that wh+rtracklength is still within active sites
        // if not then this a recombination
        GeneConversion=false;
      }
      if (!GeneConversion) {
        recombnode <T> *nw=
          lines[wh.first]->coal->recombination(lines[wh.first],wh.second,t);	    
        
        lines.push_back(nw);
        c.UpdateRecombination(wh.first);
        rr.update_recombination(wh.first,lines);
      } else {
        // recombine the first node twice
        // and then coalesce the two ends.
        recombnode <T> *nw=
          lines[wh.first]->coal->recombination(lines[wh.first],wh.second,t);	    
        lines.push_back(nw);
        rr.update_recombination(wh.first,lines);
        *nw=
          lines[wh.first]->coal->recombination(lines[wh.first],wh.second+tracklength,t);	    
        lines.push_back(nw);
        rr.update_recombination(wh.first,lines);
        ancestors.push_back(new rnode<T>(chromo,t,rr,root,lines));
        // Now fix lines 
        int removed=lines.coal(ancestors.back(),chromo);
        rr.update_coal(chromo,lines[chromo.first]->active,removed); 
      }
      
      
      collect(t,Recombination,wh.second);   //ijw to fix
#ifdef CHECK
      std::cout << "recombination for line " << wh.first
                << " at position " << wh.second << std::endl;
#endif
    }
    if (c.nleft()<=1) return;
  }
}
#endif
/** recursive function to extract the underlying tree at a given position */
template<typename T>
node<rnode<T> *> *extract_node(rnode<T> *here,int position
                                 ,const  std::vector<rnode<T> *> &samp)
{
  static int label =1000000;
  
  if (here->isleaf()) {
    int lab =  distance(samp.begin()
                          ,std::find(samp.begin(),samp.end(),here));
    
    return new node<rnode<T> *>(lab,here);
  }
  
  rnode<T> *Left = here->left(position);    
  rnode<T> *Right = here->right(position);
  
  assert(here->left(position)!=0);
  assert(here->right(position)!=0);
  
  
  return new  node<rnode<T> *>(here->time()
                                 ,extract_node(Left,position,samp)
                                 ,extract_node(Right,position,samp)
                                 ,++label,here);
}
/** extract a tree and convert to something else recursively                         */
template<typename T,typename D>
node<D> *extract_node(rnode<T> *here,int position
                        ,const  std::vector<rnode<T> *> &samp,D mydef)
{
  static int label =1000000;
  
  if (here->isleaf()) {
    int lab =  distance(samp.begin()
                          ,std::find(samp.begin(),samp.end(),here));
    
    return new node<D>(lab,mydef);
  }
  
  rnode<T> *Left = here->left(position);    
  rnode<T> *Right = here->right(position);
  assert(here->left(position)!=0);
  assert(here->right(position)!=0);
  
  return new  node<D>(here->time()
                        ,extract_node(Left,position,samp,mydef)
                        ,extract_node(Right,position,samp,mydef)
                        ,++label,mydef);
}
}
#endif   

