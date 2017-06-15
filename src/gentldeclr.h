/// @file
/// declarations used in the ARG files
#ifndef GENTL_DECLR_H_
#define GENTL_DECLR_H_

namespace GenTL {
  struct recombinations;
  typedef int position;
  template<class T,typename COLLECTOR=recombinations> class rtree;
  template<class T> class rnode;
  enum EventType{Coalescence=0, Recombination, Migration};
}

#endif
