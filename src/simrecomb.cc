// all the functions for R - I shall assume that all functions here 
// do *not* change labels from 0 to n-1 -- rather that shall
// be done in R

#include <deque>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>                // for ostringstream
#include <fstream>
#include "ARG.h"
#include "ascertain.h"
#include "growthmodel.h"
#include "migmatrix.h"
#include "mutation.h"
#include "simrecomb.h"


std::vector<int>  simrecomb::mutate(int var, const std::string &a,rng &r) 
{
  size_t sites=tr.root.size();
  // assign data - if needed
  for (size_t i=0;i<tr.sample.size();i++) {
    tr.sample[i]->data().assign(sites,0);
  } 
  // calculate lengths
  std::vector<double> len(sites);
  for (size_t i=0;i<sites;i++) len[i]=tr.length(i);
  // now mutate 
  return ascertainedMutation(a,tr,var,len,r);
}
void simrecomb::STRmutate(int pos, double theta, rng &r)
{
  if ( tr.sample[0]->data().size()==0) {
    for (size_t i=0;i<tr.sample.size();i++) {
      tr.sample[i]->data().assign(1,0);
    }
  }
  
  GenTL::mutateSTR(tr.root[pos],10,pos,theta,r);
  
  return;
}
int simrecomb::mutatetheta(double theta, int pos,rng &r)
{
  // assumes that data is already assigned
  return GenTL::mutateBinary2(tr.root[pos],0,pos,theta,r);
}

std::vector<int>  simrecomb::mutatetheta(double *theta, rng &r) 
{
  size_t sites=tr.root.size();
  for (size_t i=0;i<tr.sample.size();i++) {
    tr.sample[i]->data().assign(sites,0);
  }
  std::vector<int> var(sites);
  for (size_t j=0;j<sites;j++) {
    var[j] = GenTL::mutateBinary(tr.root[j],0,j,theta[j],r);
  }
  return var;
}

void simrecomb::mutateAboveNode(int pos) 
{
  // reassign the data
  size_t ss=tr.sample.size();
  for (size_t i=0;i<ss;i++) {
    tr.sample[i]->data().assign(ss-3,0);
  }  
  size_t count=GenTL::mutateNodes(tr.root[pos],pos);
  if (count != ss-3){
    std::ostringstream oss; 
    oss << "Count = " << count << " ss = " << ss << " - Problem with mutateAboveNode \n";
    throw std::runtime_error(oss.str().c_str());
  }
}
