/// @file 
/// templated functions for mutations
// @file      Time-stamp: <2012-01-27 17:20:19 ijw>
#ifndef MUTATION_H_IJW__
#define MUTATION_H_IJW__

#include <stdexcept>
// utility
#include "ijw_rand.h"
#include "rnodet.h"

namespace GenTL {
  /** recursively write leaves   */
  template<typename T>
  void write_leaves(rnode<T > *here, position pos)
  {
    if (here->isleaf()) { // at a leaf
      here->data()[pos]=1;
    } else {
      write_leaves(here->left(pos),pos);
      write_leaves(here->right(pos),pos);
    }
    return;
  }
  /** Write leaves when position is used to find LEft and Right, but the 
   * mutation written is at position mutno                               */
  template<typename T>
  void write_leaves(rnode<T > *here, position pos, int mutno)
  {
    if (here->isleaf()) { // at a leaf
      here->data()[mutno]=1;
    } else {
      write_leaves(here->left(pos),pos,mutno);
      write_leaves(here->right(pos),pos,mutno);
    }
    return;
  }
  /** Recursively find the position of the mutation              
   * once successful this returns  - the time at which the mutation 
   * occurred.  In this way I can keep all th4e old code, but still allow
   * the function to return the time at which the mutation occurred                                                    
   */
  template<typename T>
  double recursemutateINF(rnode<T > *here,double mtime
			  ,const position pos)
  {
    assert(!here->isleaf());
      
    rnode <T > *l(here->left(pos));
    mtime -= here->time()-l->time();
    
    if (mtime<0.0) { // mutation happens here
      mtime = -mtime-here->time();    // -1 time the time of the mutation
      write_leaves(l,pos);
      return mtime;
    } else { // no mutation here so try to the left
      if (!l->isleaf()) {
	mtime = recursemutateINF(l,mtime,pos);
	if (mtime<0.0) return mtime;
      }
    }
    // no mutation so try to the right    
    rnode <T > *r(here->right(pos));
    mtime -= here->time()-r->time();
    
    if (mtime<0.0) {
      mtime = -mtime-here->time();  // - the time of the mutation
      write_leaves(r,pos);
      return mtime;
    } else {
      if (!r->isleaf()) 
	mtime = recursemutateINF(r,mtime,pos);
    }
    return mtime;
  }
  /** First mutation function - infinite site     */
  template<typename T>
  double mutateINF(rnode<T > *root, position pos, rng &r)
  {
    double len=root->length(pos);
    double mut=r()*len;
    
    double mutated=recursemutateINF(root,mut,pos);
    if (mutated>0.0) 
      throw std::domain_error("SNPmutation failed in mutateINF");
    return -mutated;
  }
 /** Mutate above all the nodes below this root
     Assumes that there is enough memory allocated  
     Also assumes that the node *here has already been done */
  template<typename T>
  void recursemutateAboveNode(rnode<T > *here,int pos,int &count)
  {
    rnode <T > *Left(here->left(pos));
    rnode <T > *Right(here->right(pos));
    if (!Left->isleaf()) {
      write_leaves(Left,pos,count);
      count +=1;
      recursemutateAboveNode(Left,pos,count);
    }
    if (!Right->isleaf()) {
       write_leaves(Right,pos,count);
       count +=1;
       recursemutateAboveNode(Right,pos,count);
    }
  }
  /** Function to return the n-3 different arrangments of mutations
   * for a tree - we only do one of the descendents of the root as 
   * the two resulting sets are completely equivalent                           */
 template<typename T>
 int mutateNodes(rnode<T > *root,int pos)
  {
    int count=0;
    rnode <T > *Left(root->left(pos));
    rnode <T > *Right(root->right(pos));
    if (!Left->isleaf()&&!Right->isleaf()) {  // none are leaves
      write_leaves(Left,pos,0);
      count+=1; 
      recursemutateAboveNode(Left,pos,count);
      recursemutateAboveNode(Right,pos,count);
    } else if (Left->isleaf()) {
      recursemutateAboveNode(Right,pos,count);
    } else {
      recursemutateAboveNode(Left,pos,count);
    }
    return count;
  }
 /** Symmetric STR mutation function for recombination node at position pos                   
     we are putting the mutated STR at position 0*/
  template<typename T>
  void mutateSTR(rnode<T > *here,int currCount, position pos, double theta, rng &r)
  {
    assert(!here->isleaf());
    rnode <T > *Left(here->left(pos));
    rnode <T > *Right(here->right(pos));
 
    int nmut = r.rpoisson(theta*(here->time()-Left->time())/2.);
    int countLeft=currCount;
    if (nmut>0)  countLeft += 2*r.rBinom(0.5,nmut) - nmut;

     if (!Left->isleaf()) 
      mutateSTR(Left,countLeft,pos,theta,r);
    else 
      Left->data()[0]=countLeft;

    nmut = r.rpoisson(theta*(here->time()-Right->time())/2.);
     int countRight=currCount;
    if (nmut>0) countRight += 2*r.rBinom(0.5,nmut)-nmut;
 
    if (!Right->isleaf()) 
      mutateSTR(Right,countRight,pos,theta,r);
    else 
      Right->data()[0]=countRight;
  }

  /** A mutation function that only considers 2 classes (so just considers transitions) 
      T can be any type but we just assign 0 and 1  returns a bool that determines whether a
      site has mutated (but could still not be variable */
  template<typename T>
  int mutateBinary(rnode<T > *here,int curr, position pos, double theta, rng &r)
  {
    int muts=0;
    assert(!here->isleaf());
    rnode <T > *Left(here->left(pos));
    rnode <T > *Right(here->right(pos));

    int nmut = r.rpoisson(theta*(here->time()-Left->time())/2.);
    int bLeft=(nmut%2==1)?1-curr:curr;
    muts+=nmut;
    if (!Left->isleaf()) {
      muts += mutateBinary(Left,bLeft,pos,theta,r);
    } else  Left->data()[pos]=bLeft;

    nmut = r.rpoisson(theta*(here->time()-Right->time())/2.);
    int bRight=(nmut%2==1)?1-curr:curr;
    muts+=nmut;
    if (!Right->isleaf()) {
      muts+=mutateBinary(Right,bRight,pos,theta,r);
    }
    else Right->data()[pos]=bRight;
    
    return muts;
  }

  /** A mutation function that only considers 2 classes (so just considers transitions) 
      T can be any type but we just assign 0 and 1  returns a bool that determines whether a
      site has mutated (but could still not be variable 
      This is the same as mutateBinary but now all the results are put at position 0
*/
  template<typename T>
  int mutateBinary2(rnode<T > *here,int curr, position pos, double theta, rng &r)
  {
    int muts=0;
    assert(!here->isleaf());
    rnode <T > *Left(here->left(pos));
    rnode <T > *Right(here->right(pos));

    int nmut = r.rpoisson(theta*(here->time()-Left->time())/2.);
    int bLeft=(nmut%2==1)?1-curr:curr;
    muts+=nmut;
    if (!Left->isleaf()) {
      muts += mutateBinary2(Left,bLeft,pos,theta,r);
    } else  Left->data()[0]=bLeft;

    nmut = r.rpoisson(theta*(here->time()-Right->time())/2.);
    int bRight=(nmut%2==1)?1-curr:curr;
    muts+=nmut;
    if (!Right->isleaf()) {
      muts+=mutateBinary2(Right,bRight,pos,theta,r);
    }
    else Right->data()[0]=bRight;
    
    return muts;
  }
}

#endif

