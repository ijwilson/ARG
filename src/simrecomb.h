#include <Rcpp.h>
#include "ARG.h"

class simrecomb {
public: 
  typedef GenTL::rnode<std::vector<int> > nodetype;
  typedef GenTL::ARG<std::vector<int> ,GenTL::recombinationLocations > ARGtype;
  std::vector<int>  mutatetheta(double *theta, rng &r);
  int mutatetheta(double theta, int pos,rng &r);
  simrecomb(int ss, int sites,double rec,rng &r):
    g(new GenTL::constant_size())
    ,c(new GenTL::coalescent(ss,g,r))
    ,tr(sites,std::vector<double>(sites-1,rec),*c
          ,r,recomb) {
  };
  
  simrecomb(int ss, int sites,const std::vector<double> &rec
              , std::string &gm,rng &r):
    c(new GenTL::coalescent(ss,gmread(gm),r)),
    tr(sites,rec,*c,r,recomb) {
  };
  
  simrecomb(const std::vector<double> &rec
              , const std::vector<int> &location
              ,std::string &gm
              , std::string &migm
              , rng &r):
    g(gmread(gm))
    ,mm(mmread(migm))
    ,c(new GenTL::structured_coalescent(location,g,*mm,r))
    ,tr(rec.size()+1,rec,*c,r,recomb) {
  };
  std::vector<int> mutate(int var,const std::string &a,rng &r);
  void STRmutate(int pos, double theta, rng &r);
  
  void mutateAboveNode(int pos);
  
  GenTL::growthmodel *g;
  GenTL::mig_matrix *mm;
  GenTL::coalescent *c;
  GenTL::recombinationLocations recomb;
  ARGtype tr;
};

