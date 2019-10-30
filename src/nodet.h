/** @file  */ 
#ifndef __NODE__H__T_IJW
#define __NODE__H__T_IJW

#include <stack>
#include "gsl_distributions.h"
#include "coalescent.h"


namespace GenTL {

  /** \brief template node class.  
   *
   * A class for buildins simple trees.  the class for extracting trees
   * from rnodes is built directly from this class
   */

  template <class T>
  class node {
    //  friend class tree;
    //    friend class iterator;
    typedef node* nodeptr; /// pointer to a node
  public:
    node()                                        /// default constructor
      :left(0),right(0),time(0.0),data_(){ };                            
    node(int lab,T data)                          /// leaf constructor
      :left(0),right(0),time(0.0),label(lab),data_(data){}; 
    node(int lab)                          /// leaf constructor
      :left(0),right(0),time(0.0),label(lab){}; 
    node(T data)                                  /// leaf constructor
      :left(0),right(0),time(0.0),data_(data){ };                   
    node (double t, nodeptr l, nodeptr r)
      :left(l),right(r),time(t) {
      left->up=this;
      right->up=this;
    };
    node (double t, nodeptr l, nodeptr r,int lab,T &a)
      :left(l),right(r),time(t),label(lab),data_(a){
      left->up=this;
      right->up=this;
    };
    node (double t, nodeptr l, nodeptr r,T &data)
      :left(l),right(r),time(t),data_(data){
      left->up=this;
      right->up=this;
    };
    
    void assign(double t, nodeptr l, nodeptr r,T &a){
      left=l;
      right=r;
      time=t; // assumes data already correct
      left->up=this;
      right->up=this;
      data_=a;
    }
    void recursiveprint(std::ostream &o,bool quote=true, bool internalnodes=false) const; 
    void recursiveprintnodata(std::ostream &o,bool printtimes=true,bool internalnodes=true) const;
    ///recursively print the tree in Newick format.
    int classifyancestralnode(const std::vector<int> &location);
    int count_leaves() const {
      if (left==0) return 1;
      else return left->count_leaves()+right->count_leaves();
    }
    /** get a vector of labels below a node       */
    void leaf_labels(std::vector<int> &leaves) {
       if (left==0) {
	 leaves.push_back(label);
       } else {
	 left->leaf_labels(leaves);
	 right->leaf_labels(leaves);
       }
       return;
    }
    /** get the list of leaves below a node       */
    void get_leaves(std::vector<node<T> *> &leaves) {
       if (left==0) {
	 leaves.push_back(this);
       } else {
	 left->get_leaves(leaves);
	 right->get_leaves(leaves);
       }
       return;
    }
    /** change everything at position pos to value   */
    void change_below(int pos, char value);


    /** get the list of internal node below (and including) a node       */
    void get_internal(std::vector<node<T> *> &a) {
       if (left!=0) {
	 a.push_back(this);
	 left->get_internal(a);
	 right->get_internal(a);
       } 
       return;
    }
    /** calculate the length of tree below a node */
    double length() const {
      if (left==0) return 0.0;
      double t=2.*time-left->time-right->time;
      return t+right->length()+left->length();
    }
    /** is this node a leaf ?                                            */
    bool isleaf() const {
      return (left==0);
    } 
    /** Find a node with this label                  */
    node<T> *find_node(int lab) {
      if (isleaf()) {
	if (label==lab) return this;
	return 0;
      }
      node <T> *ret=left->find_node(lab);
      if (ret==0) ret=right->find_node(lab);
      return ret;
    };
    void write_leaves(int pos);


    double recursemutate(double mtime ,const int pos);
    void recursemutate(const std::vector<double> &theta, rng &r);
    void recursemutate(size_t start,const std::vector<double> &theta, rng &r);
    
    void nodefreq(int *nf, int npops);

    void print(std::ostream &o){ o << data_;}
    nodeptr left;  /// pointer to "left" node 
    nodeptr right; /// pointer to "right" node
    nodeptr up;  /// pointer to node above 
    std::vector<int> countleft;
    std::vector<int> countright;
    double time;   /// coalescence time
    int label;     /// label for the node
    const T &data() const {return data_;}; ///extract the data (but not change it).
    T &data() { return data_;}; //extract the data
  
    //  friend node<T> readnode(std::istream &in);
  protected:
    T data_;
  };

  /** forward declaration of iterator                       */
  template <class T> class iterator; 
  /** tree class                                            */
  template <class T> 
  class tree {
  public:
    /** default constructor                        */
    tree() {
    };
    /** constructor for n haplotypes               */
    tree(int n):sample_(n),ancestors_(n-1) {
      root = &(ancestors_[n-2]);
    };
    /** constructor for a coalescent               */
    tree(coalescent *c):
      sample_(c->nleft())
      ,ancestors_(c->nleft()-1) {
      root=&(ancestors_[c->nleft()-2]);
      simulate(c);
    }; 
    /** constructor for a branching process               */
    tree(int n, ctsdistribution &nexteventprior, rng &r):
      sample_(n),ancestors_(n-1) {
      root=&(ancestors_[n-2]);
      bpsimulate(n,nexteventprior,r);
    };

    /** declaration for the simulation             */		
    void simulate(coalescent *c);
    /** declaration for the branching process simulation             */
    void bpsimulate(int n, ctsdistribution &nexteventprior, rng &r);
    /** iterate through the tree                   */
    iterator<T> begin() {
      iterator<T> i(*this);
      return i;
    };
    /** remakes a tree from the root */
    node<T> *root;
 
    void getnodefreq(std::vector<int> &location, int npops);
    /** Count to the left and right from without locations    */
    void getnodefreq();

 
    node<T> *readnode(std::istream &in, std::vector<node<T> > &anc
		      , std::vector<node<T> > &sample);
    
    iterator<T> end();

    node<T> *first() {
      return &sample_[0];
    }

    const node<T> *ancestorptr(int i) const {
      return &ancestors_[i];
    }

    node<T> *ancestorptr(int i)  {
      return &ancestors_[i];
    }

    const node<T> &ancestor(int i) const {
      return &ancestors_[i];
    }

    node<T> sample(int i) const {
      return sample_[i];
    }

    const T &data(int i) const {
      return sample_[i].data();
    }

    T &data(int i) {
      return sample_[i].data();
    }

    double anctime(int i) const {
      return ancestors_[i].time;
    }

    /** How many leaves dows the tree have ?                                  */
    size_t nleaves() const {
      return sample_.size();
    }
    /** mutate - template  specialisation for different types of node         */
    void mutate(int nloci, rng &r);
    void mutate(const std::vector<int> &nmuts, rng &r);
    void mutate(const std::vector<double> &theta, rng &r);
    void mutate(size_t start, const std::vector<double> &theta, rng &r);
    std::vector<int> haplogroups(int pos,double minsplit);

 



    /** friend to the input operator                  */
    template<class V>
    friend  std::istream &operator>>(std::istream &in,  tree<V>  &a);
  protected:
    std::vector<node<T> > sample_;
    std::vector<node<T> > ancestors_;
  private:
    tree(const tree &a);
  };
  /** iterator class
   *
   */
  template <class T> class iterator {
  public:
    iterator(const tree<T> &t) {
      add(t.root);
    };
    node<T> *operator++();

    node<T> *operator *() {	 
      return where_.top();
    };
    bool operator != (const node<T> *tmp) {
      if (where_.empty())  return (tmp!=0);
      else return (where_.top()!=tmp);
    };
               
  private:
    void add(node<T> *here){
      for (;;) {
	if (here==0) break;
	where_.push(here);
	here=here->left;
      } 
    };
    /** we keep a stack of the  nodes.  The current position
     * of the iterator is given by where_.top() */
    std::stack<node<T> *> where_;
  };

  /**
   * output a tree on ostream 
   */
  template <class T>
  std::ostream &operator<<(std::ostream &o,const  node <T> &root) {
    root.recursiveprint(o);
    return o;
  }
  /**
   * output a tree on ostream 
   */
  template <class T>
  std::ostream &operator<<(std::ostream &o,const  tree<T> &t) {
    t.root->recursiveprint(o,true,true);
    o << ";";
    return o;
  }

  /** recursive read nodes         */
  template <class T>
  node<T> *readnode(std::istream &in, std::vector<node<T> > &anc, std::vector<node<T> > &sample
		    , bool quote=true, bool internalnodes=false,bool uselabel=true)
  {
    int ch=skipblank(in);
    if (ch=='(') {
        node<T> *l,*r;
        l=readnode<T>(in,anc,sample,quote,internalnodes,uselabel);
        checkreadcharacter(in,':',"after left node");
    
        double ltime,rtime;
        in >> ltime;
    
        checkreadcharacter(in,',',"comma after time");  
        r=readnode<T>(in,anc,sample,quote,internalnodes,uselabel);
    
        checkreadcharacter(in,':',"reading right node");
        in >> rtime;
    
        if (fabs(ltime+l->time-rtime- r->time)>1E-5) 
        throw std::domain_error("error with times - I expected a ... tree");     
    
        checkreadcharacter(in,')',"error after reading right node");
        int lab(0);
        T d(0);
    
        if (internalnodes) {
            ch=skipblank(in);
            if (quote) {
                if (ch!='\'') {

                    in.putback(ch);
                } else {
                    if (uselabel) in >> lab;
                    in >> d;
                    std::cerr << "read " << lab << " and data " << d << std::endl; 
                    checkreadcharacter(in,'\'',"error - expected quote after reading internal node");  
                }
            } else {
                in.putback(ch);
                if (uselabel) in >> lab;
                in >> d;
            }
        }
        anc.push_back(node<T>(ltime+l->time,l,r,lab,d)); 
        return &anc[anc.size()-1];
    } else {
        int lab(0);
        T d(0);
        if (quote) {
            if (ch!='\'') {

                in.putback(ch);
            } else {
                if (uselabel) in >> lab;
                in >> d;
                std::cerr << "read " << lab << " and data " << d << std::endl; 
                checkreadcharacter(in,'\'',"error - expected quote");  
            }
        } else {
            in.putback(ch);
            if (uselabel) in >> lab; 
            in >> d;
        }
        sample.push_back(node<T>(lab,d));
        return &sample[sample.size()-1];
    }
  }
  /**
   * input a tree from stream
   */
  template <class T>
  std::istream &operator>>(std::istream &in,  tree <T> &a) {
    a.sample_.resize(0);
    a.sample_.reserve(1000);
    a.ancestors_.reserve(1000);
    a.ancestors_.resize(0);
    a.root=readnode<T>(in,a.ancestors_,a.sample_,true,true);
    return in;
  }

  template <class T>
  node<T> *iterator<T>::operator++() {
    if (where_.empty()) return 0;
    node<T> *curr=where_.top(); 
    where_.pop();
    if (curr->right !=0) {
      add(curr->right);
    }
    return curr;
  };
  template <class T>
    void node<T>::recursiveprint(std::ostream &o, bool quote, bool internalnodes) const 
  {
    if (left!=0) {
      o << "(";
      left->recursiveprint(o,quote,internalnodes);
      o << ":" << time-left->time << ",";
      right->recursiveprint(o,quote,internalnodes);
      o << ":" << time-right->time << ")";
      if (internalnodes) {
	if (quote) o << "'";
	o << label << " " << data_;
	if (quote) o << "'";
      }
    } else {
      if (quote) o << "'";
      o << label << " " << data_;
      if (quote) o << "'";
    }
  }

 template <class T>
   void node<T>::recursiveprintnodata(std::ostream &o,bool printtime,bool internalnodes) const  {
    if (left!=0) {
      o << "(";
      left->recursiveprintnodata(o,printtime,internalnodes);
      if (printtime) o << ":" << time-left->time;
      o << ",";
      right->recursiveprintnodata(o,printtime,internalnodes);
      if (printtime) o << ":" << time-right->time;
      o << ")";
      if (internalnodes) o << label;
    } else {
      o << label;
    }
  }
  /** What sort of nodes is this node ancestral to?
   * types of location i only ( returns i) 
   * or to a mixture of different nodes - returns -1        */
  template <class T>
  int node<T>::classifyancestralnode(const std::vector<int> &location)
  {
    int l,r;
    if (left!=0) {
      int l= left->classifyancestralnode(location);
      if (l==-1) return -1;
      int r=right->classifyancestralnode(location);
      if (l==r) return l;
      else return -1;
    }
    else {
      return location[label];
    }
  }
  /** A utility function to do the simulations for all constructors                             */
  template <class T>
  void tree<T>::simulate(coalescent *c)
  {
    double t=0.0;
    std::vector<node<T> *> lines(c->nleft());
    for (int i=0;i<c->nleft();i++) lines[i]=&sample_[i];
    int pos=0;    
    for (;;) {
      double NextEvent=c->next(t);
      if (NextEvent>1E5) {
        std::ostringstream oss;
        oss << "Problem with coalescent: " << NextEvent << std::endl;
        throw std::runtime_error(oss.str().c_str());
      }
      t += NextEvent;
     
      std::pair<int,int> chromo = c->coalchromo();
      if (chromo.first!=chromo.second) { // a coalescence 
	ancestors_[pos].left=lines[chromo.first];
	ancestors_[pos].right=lines[chromo.second];
	ancestors_[pos].time=t;
	lines[chromo.second]=&ancestors_[pos];
	lines[chromo.first]=lines[c->nleft()-1];
	pos++;
	c->UpdateCoalescence(chromo,1);
	if (c->nleft()==1) break;
      }
    }
    return;
  }
  template <class T> 
  void tree<T>::getnodefreq(std::vector<int> &location, int npops)
  {
    for (int i=0;i<ancestors_.size();i++) {
      ancestors_[i].countleft.resize(npops,0);
      ancestors_[i].countright.resize(npops,0);
      if (ancestors_[i].left->left==0) { // a leaf
	ancestors_[i].countleft[location[&sample_[0]-ancestors_[i].left]] =1;
      } else {
	for (int j=0;j<npops;j++) ancestors_[i].countright[j]=
	  ancestors_[i].right->countleft[j]+ancestors_[i].right->countright[j];
      }
      if (ancestors_[i].right->left==0) { // a leaf
	ancestors_[i].countright[location[&sample_[0]-ancestors_[i].right]] =1;
      } else {
	for (int j=0;j<npops;j++) ancestors_[i].countleft[j]=
	  ancestors_[i].left->countleft[j]+ancestors_[i].left->countright[j];
      }
    }
  }
 template <class T> 
  void tree<T>::getnodefreq()
  {
    for (size_t i=0;i<ancestors_.size();i++) {
      ancestors_[i].countleft.resize(1,0);
      ancestors_[i].countright.resize(1,0);
      if (ancestors_[i].left->left==0) { // a leaf
	ancestors_[i].countleft[0] =1;
      } else {
	ancestors_[i].countright[0]=
	  ancestors_[i].right->countleft[0]+ancestors_[i].right->countright[0];
      }
      if (ancestors_[i].right->left==0) { // a leaf
	ancestors_[i].countright[0]=1;
      } else {
	 ancestors_[i].countleft[0]=
	  ancestors_[i].left->countleft[0]+ancestors_[i].left->countright[0];
      }
    }
  }
  template <typename T>
    void recursively_destroy_node(node <T> *here) {
    if (!here->isleaf()) {
      recursively_destroy_node(here->left);
      recursively_destroy_node(here->right);
    }
    delete here;
  }

  /// rotate the labels on a tree so that the smallest label is to the left
 template <typename T>
   int rotate_tree_labels(node <T> *here) {
   if (here->isleaf()) return here->label;
   int lft=rotate_tree_labels(here->left);
   int rgt=rotate_tree_labels(here->right);
   if (lft>rgt) {
     node<T> *tmp=here->left;
     here->left=here->right;
     here->right=tmp;
     return rgt;
   }
   return lft;
 }

}
#endif
