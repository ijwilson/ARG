/**    @file */
#ifndef GSL_SORT_H__
#define GSL_SORT_H__
//#include <vector> - we can't need vector if we have to pass one!


#include "Rmath.h"


#undef sexp

class gsl_index {
public:
  gsl_index() {
  }
  gsl_index(const std::vector<double> &x) {
    sort(x);
  }
  ~gsl_index() {
    if (p) delete[] p;
  }
  void sort(const std::vector<double> &x) {
     
    p=new int[x.size()];
    std::vector<double> y(x);
    rsort_with_index(&y[0],p,x.size());
  }
  size_t operator[](int i ) const {
    return p[i];
  }
  private:
  int *p;
};


#endif
