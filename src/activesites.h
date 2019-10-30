/** @file */
#ifndef ACTIVE_SITES_H__
#define ACTIVE_SITES_H__

#include <set>
#include <iterator>
#include <list>
#include <iostream>
#include <algorithm>
#include <cassert>

/// makes it plain where we are talking about position on a sequence
typedef int position;

#define ACTIVERANGE

namespace GenTL {
#ifdef ACTIVERANGE


  /** A structure for representing ranges          */
  template <typename VAL>
    struct range { 
      range(const VAL &m,const  VAL &x):mn(m),mx(x) {
      };
      /** is rhs contained in the range ?          */
      bool contains(const range &rhs) {
	if (mn<=rhs.mn and mx >= rhs.mx) return true;
	return false;
      };
      /** is range contained in the rhs ?          */
      bool iscontained(const range &rhs) {
	if (mn>=rhs.mn and mx <= rhs.mx) return true;
	return false;
      };
      /** print a representation of a range        */
      std::ostream &print(std::ostream &o) const {
	o <<  mn <<"-"<<mx;
	return o;
      }
      /** The minimum and maximum values   */
      VAL mn,mx;
    };
  /** Are two ranges the same ?                    */
  template <typename VAL>
    bool operator==(const range<VAL> &lhs,const range<VAL> &rhs)  {
    if (lhs.mn==rhs.mn and lhs.mx==rhs.mx) return true;
    return false;
  }
  
  /** do rhs and lhs overlap ?                   */
  template <typename VAL>
    bool overlaps(const range<VAL> &a,const range<VAL> &b)  {
    if (a.mn<=b.mn) {
      if(b.mn<=a.mx) return true;
      return false;
    } else {
      if (a.mn <= b.mx) return true;
      return false;
    }
  }
  /** ostream<< operator for the range class     */
  template <typename VAL>
    std::ostream &operator<<(std::ostream &o, const range<VAL> &a) {
    o <<  a.mn <<"-"<<a.mx;
    return o;
  }
  /** The main activetype class                  */
  class activetype {
  public:
    typedef std::list<range<position> >::iterator itor;
    typedef std::list<range<position> >::const_iterator const_itor;
    /** default constructor   */
    activetype():lr() {
    };
    /** is it empty?          */
    bool empty() const {
      return lr.empty();
    }
    /** used in the constructor to initialise */
    void insertrange(position first, position end) {
      lr.insert(lr.begin(),range<position>(first,end-1));
    }
    /** is the position x present?            */
    bool present(const position &x) {
      itor i=begin();
      while (i!=end()) {
	if ((*i).mn<=x) {
	  if (x<=(*i).mx) return true;
	  i++;
	} else return false;
      }
      return false;
    };
    /** return all the positions present  */
    void present(std::list<position> &lp) {
      itor i =begin();
      while (i!=end()) {
	position s=(*i).mn;
	position e=(*i).mx+1;
	while(s!=e) lp.push_back(s++);	
	i++;      
      }
    }
   /** return all the positions present  */
    std::list<position> activesites() {
      std::list<position> lp;
      itor i =begin();
      while (i!=end()) {
	position s=(*i).mn;
	position e=(*i).mx+1;
	while(s!=e) lp.push_back(s++);	
	i++;      
      }
      return lp;
    }
    /** The end iterator                 */
    const_itor end() const {
      return lr.end();
    }
    itor end() {return lr.end();}

    const_itor begin() const {
      return lr.begin();
    }

    itor begin() {
      return lr.begin();
    }
 
    position first() const {
      return lr.front().mn;
    };
    /** The last position   */
    position last() const {
      return lr.back().mx;
    };

    void insert(std::list<position> &intersect, const range<position>  &addr) {
      
      if (addr.mn>lr.back().mx) {
	lr.push_back(addr);
	return;
      }
	
      itor i=begin();
      while (i != end()) {
	if (addr.mx<(*i).mn) {
	  lr.insert(i,addr);
	  return;
	}  
	if (overlaps(*i,addr)) {
	  position start=std::max(int((*i).mn),int(addr.mn));
	  position last =std::min(int((*i).mx),int(addr.mx));
	  while (start!=last+1) 
	    intersect.push_back(start++);
	  if ((*i).mn>addr.mn) (*i).mn=addr.mn;
	  if ((*i).mx<addr.mx)  (*i).mx=addr.mx;
	  return;
	}
	i++;
      }
      assert(false); // should never get here
    }
    /** fix the intersections  */
    void fix(std::list<position> &insect) {
      if (lr.size()<2) return;      
      itor i=begin();
      itor j=i;j++;      
      while (j!=end()) {
	assert(i!=j);
	if ((*i).mx>=(*j).mn) {
	  position start=(*j).mn;
	  position last =std::min(int((*i).mx),int((*j).mx));
	  while (start!=last+1) 
	    insect.push_back(start++);
	  if ((*j).mx>(*i).mx)
	    (*i).mx=(*j).mx;
	  lr.erase(j);
	  if (lr.size()==1) return;
	  j=i;
	  j++;
	} else {
	  ++i;
	  ++j;
	}
      }
    }
    /** Combine two sets of active sites             
     * returns the intersection of then (i.e. those sites
     * where we loose an ancestral site               */
    std::list<position> combine(const activetype &l, 
				const activetype &r) {
   
      std::list<position> intersect;
      // first put a copy of l in the activetype 
      copy(l.begin(),l.end(),back_inserter(lr));
#ifdef ACTIVECHECK
      std::cout << "have:";
      print(std::cout);
      std::cout << "\n adding:";
      r.print(std::cout);
      std::cout << std::endl;
#endif
      // now try to put r into the activetype
      const_itor right=r.begin();
      while (right != r.end()) insert(intersect,*right++);
#ifdef ACTIVECHECK 
        std::cout << "before fixingm, activetype \n";
       print(std::cout);
       std::cout << std::endl;
#endif
      fix(intersect);
#ifdef ACTIVECHECK 
      std::cout << "new activetype \n";
       print(std::cout);
       std::cout << "\n intersect: ";
       list<position>::iterator j=intersect.begin();
       while (j!=intersect.end()) std::cout << *j++ <<" ";
       std::cout << std::endl;
#endif
      return intersect;
    }
    /** is there only a singleton                    */
    bool  singleton() const {
      if (lr.front().mn==lr.back().mx) return true;
      return false;
    }
    /** Erase the positon                            */
    void erase(position pos) {
      itor i=begin();
      while (i != end()) {
	assert((*i).mn<=pos);  // necessary condition
	if ((*i).mn==pos) {
	  if ((*i).mx==pos) {
	    lr.erase(i);
	  } else  (*i).mn+=1;
	  break;
	} 
	if ((*i).mx >=pos) {
	  if ((*i).mx==pos) {
	    (*i).mx-=1;
	  } else {
	    lr.insert(i,range<position> ((*i).mn,pos-1));
	    (*i).mn=pos+1;
	  }
	  break;
	}
	i++;
      }  
      assert(i!=end());
    } 
    /** Split the actives sites at position wh      */
    void split(position wh, activetype &old) {
      itor i=old.begin();
      while (i != old.end()) { // need to split i and 
        if ((*i).mn <= wh && (*i).mx > wh) {
          lr.insert(i,range<position> ((*i).mn,wh));
          (*i).mn=wh+1;
          lr.splice(begin(),lr,old.begin(),i);
          break;
        } else if ((*i).mn > wh) {  // split at i
          lr.splice(begin(),lr,old.begin(),i);
          break;
        }
        i++;
      }
    }
    /** Print the active sites                       */    
    std::ostream &print(std::ostream &o) const {
      const_itor i=lr.begin();
      while (i!=lr.end()) {
	o << *i << " ";
	i++;
      }
      o  << std::endl;
      return o;
    };
  private:
    /** A list of ranges that gives the active sites  */ 
    std::list<range<position> > lr;
  };



/** A class that acts like an iterator 
 * I don't want to create a fully fledged iterator
 * here yet as I don't need all the functionality  
 * Note that we also need an iterator through the 
 * recombnodes to find all the positions at a node */ 
/* class posItor { */
/*   friend class activetype; */
/*  public:  */
/*   posItor(const activetype &a) { */
/*     _current=a.begin(); */
/*     _curpos=(*_current).mn; */
/*     return _curpos; */
/*   } */

/*   operator==(const activetype::itor &ai) { */
/*     if (_current==end */
  
/*   operator++() { */
/*     if (_curpos<(*_current).mx) { */
/*     _curpos++; */
/*     } else if (_curpos==(*current_.mx) */
/*     return _curpos; */
/*   }; */
  
/*   operator++(int) { */
/*      position tmp=_curpos; */
/*      operator++(); */
/*      return tmp; */
/*   }; */
  
/*  private: */
/*   activetype::itor _current; */
/*   position _curpos; */
/* }; */



#endif

#ifdef ACTIVESET

typedef int position;
/** class activetype holds the set of active sites at a node
 * initially inherited from the base set class for simplicity
 * but this should allow larger changes in time             */
 class activetype {
 public:
   typedef std::set<position>::iterator iterator;
   /** default constructor                            */
   activetype():sp() {
   };
   /** is it empty?                                   */
   bool empty() const {
     return sp.empty();
   }
   /** insert a range of values (for initialisation   */
   void insertrange(position first, position end) {
     for (position i=first;i<end;i++) sp.insert(i);
   }
   /** is the position ky present?                    */
   bool present(const position &ky) const {
     return bool(sp.count(ky));
   }
   /** The first active position                      */
   position first() const {
     return *sp.begin();
   }
   /** And the last active position                   */
   position last() const {
     return *(--sp.end());
   }
   /** Combine two sets of active sites             
    * returns the intersection of then (i.e. those sites
    * where we loose an ancestral site               */
   std::list<position> combine(const activetype &l, 
			       const activetype &r) {
     set_union(l.sp.begin(),l.sp.end() 
	       ,r.sp.begin(),r.sp.end() 
	       ,inserter(sp,sp.begin()));
     std::list<position>  tmp;
     set_intersection(l.sp.begin(),l.sp.end() 
		      ,r.sp.begin(),r.sp.end() 
		      ,inserter(tmp,tmp.begin()));
     return tmp;
   }
   /** is there only a singleton                    */
   bool singleton() const {
     return (sp.size()==1);
   }
   /** Erase the positon                            */
    void erase(position pos) {
      assert(sp.count(pos)>0);
      sp.erase(pos);
    }
    /** Split the actives sites at position wh      */
    void split(position wh, activetype &old) {
        iterator split=old.sp.lower_bound(wh);
	if (*split!=wh) { // not present so more tricky
	  sp.insert(old.sp.begin(),split);
	  old.sp.erase(old.sp.begin(),split);
	} else {
	  sp.insert(old.sp.begin(),++split);
	  old.sp.erase(old.sp.begin(),split);
	}
    }
    /** Print the active sites             */    
    std::ostream &print(std::ostream &o) const {
      iterator i=sp.begin();
      while (i!=sp.end()) o << *(i++) <<"-";
      return o << std::endl;
    };
 private:
    /** keep the information as a set              */
    set<position> sp;
};
#endif



}



#endif
