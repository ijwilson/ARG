/** @file
 * A set of functional classes that can be useful  
*/
#ifndef UTILITYFUNCTIONALS_H_
#define UTILITYFUNCTIONALS_H_
#include <R.h>

#include <vector>
#include <algorithm>
#include <numeric>   // for accumulate
#include <set>
#include <cmath>
#include <iosfwd>
#include <iostream>
#include <iterator>  // for ostream_iterator
#include <sstream>

template <typename ITOR> std::vector<int> tabulate(ITOR first, ITOR last);
template <typename T> std::ostream &print(std::ostream &o
		    ,const std::vector<std::vector<T> > &a);
template <typename T> T rescale(  T *a, int count);
template <typename ITOR> std::vector<std::vector<int> > tabulate(ITOR first, ITOR last, ITOR f2, ITOR l2);
template <typename T> void keep_first_copy(std::vector<T> &a);



/// A class to act as a cumulative sum functional
template <typename T>
class cumsum {
public:
    cumsum():sum(0){}
    T& operator()(const T &item) {
	sum += item;
	return sum;
    }
private:
    T sum;
};

/// functional to help delete vectors of pointers
/// from Meyers, Effective STL p38
/// foreach
struct DeleteObject {
  template <typename T>
  void operator()( const T* ptr) const
  {
    delete ptr;
  }
};

/** A class for set comparators   */
//template<class T>
//struct less {
//  bool operator()(T x, T y) { return x<y; }
//};
template<class T>
struct moreorequal {
  bool operator()(T x, T y) { return !(x<y); }
};

/// rescale a list of doubles by dividing by the sum
/** returns tha log of the sum.  Used when working
 * with small numbers, for example for HMMs
 */
template <typename T>
T rescale( T *a, int count)
{
  T sm = std::accumulate(a,a+count,0.0);
  sm /= (double)count;
  for (int i=0;i<count;i++) a[i]/=sm;
  return log(sm);
}

// only written for character variables
template <typename ITOR>
std::vector<int> tabulate(ITOR first, ITOR last)
{
  std::vector<int> b(*(std::max_element(first,last))+1);
  ITOR i=first;
  while (i!=last) {
     b[*i]+=1;
     i++;
  }
  return b;
}

/** cross tabulation using iterators     */
template <typename ITOR>
std::vector<std::vector<int> >tabulate(ITOR first, ITOR last, ITOR f2, ITOR l2)
{
  int cols = *(std::max_element(f2,l2));
  int rows =*(std::max_element(first,last));
  std::vector<std::vector<int> > b(rows+1);
  for (size_t i=0;i<b.size();i++) b[i].assign(cols+1,0);
 

  ITOR i=first;
  ITOR j=f2;
  while (i!=last) {
     b[*i][*j]+=1;
     i++;j++;
  }
  return b;
}

template <typename T>
std::ostream &print(std::ostream &o
		    ,const std::vector<std::vector<T> > &a)
{
  for (int i=0;i<a.size();i++) {
    std::copy(a[i].begin(),a[i].end()
	      ,std::ostream_iterator<T >(o," "));
    o << std::endl;
  }
  return o;
}


/// class that simulates a factorial function- but includes a lookup table.
class lfactorl {
public:
  lfactorl() {
    for (int i=0;i<32;i++) a[i]=0.0;
  }
  double operator()(int n) {
    if (n<=1) return 0.0;
    if (n<=31) return (a[n]>0.0) ? a[n] : (a[n]=lgamma((double)(n+1)));
    else return lgamma((double)(n+1));
  };
 private:
  double a[32];
};

/// Class to keep success rates of MCMC update steps
class success_rate {
 public:
  /// constructor - set everything to zero
  success_rate() {
    success=trials=0.0;
  }
  // Add a value (true if a success)
  void operator+=(bool att) {
    if (att) success+=1.0;
    trials+=1.0;
  }
  void operator++() {
    success+=1.0;
    trials +=1.0;
  }
  void operator--() {
    trials+=1.0;
  }
  double operator()() const {
    if (fabs(trials)<0.0001) return 0.0;
    return success/trials;
  }
  void clear() {
    success=0.0;
    trials=0.0;
  }
 private:
  double success;
  double trials;
};


/** A template function that removes all but the last copy of
    a value in a vector - there must be a better way as this is
    order I am not happy with this as the original iterator based one
    was better but failed when iterators fell off the end
*/
template <typename T>
void keep_first_copy(std::vector<T> &a)
{
  // first remove adjacent values - cause problems for the rest
  // of the algorithma and this is more efficient
  a.erase(std::unique(a.begin(),a.end()),a.end());
  std::set<T> s;
  typename std::vector<T>::iterator i=a.begin();
  while (i!=a.end()) {
    if (s.count(*i)==0) {
	  s.insert(*i);
	  i++;
    } else {
      i=a.erase(i);
    }
  }
}
/// A very simple warning class
class warning {
public:
  /** Only consists of a constructor that prints a message */
  warning (const char *text) {
    oss  << text << std::endl;
    print();
  };
  template <typename T>
  warning (const char *text,const T &p1,const char *end=0) {
    oss << text << p1 << end << std::endl;
    print();
  } 
  template<typename T,typename R>
  warning(const char *text, const T p1,const char *text2, const R p2, const char *text3=0) {
    oss << text << p1 << text2 << p2 << text3 << std::endl;
    print();
  }
  template<typename T,typename R,typename S>
  warning(const char *text, const T p1,const char *text2, const R p2, const char *text3
          , const S p3, const char *text4=0) {
    oss << text << p1 << text2 << p2 << text3 << p3 << text4 << std::endl;
    print();
  }
private:
  std::ostringstream oss;
  void print() {
    Rprintf("%s",oss.str().c_str());

  }
};

#endif
