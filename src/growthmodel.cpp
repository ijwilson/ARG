#include "growthmodel.h"
#include "newio.h"
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

GenTL::growthmodel *gmread(const std::string &s)
{
  std::vector<std::string> tok;
  Tokenize(s,tok,"(), ");
    if (tok[0]=="constant") {
	if (tok.size()==1) 
	  return new GenTL::constant_size();
	else {
	  throw std::domain_error("Need no parameters for const size");
	} 
    } else if (tok[0]=="exponential" or tok[0]=="Exponential") {
	if (tok.size()!=2) 
	    throw std::domain_error("Need one parameter for exponential growth");
	std::vector<double>x=stringtodouble(tok,1);
	return new  GenTL::exponentialgrowth(x);
    } else if (tok[0]=="expfrombase") {
	if (tok.size()!=3) 
	    throw std::domain_error("Need two parameters for expfrombase");
	std::vector<double>x=stringtodouble(tok,1);
	return new  GenTL::expfrombase(x);
    } else if (tok[0]=="bottleneck") {
	if (tok.size()!=4) 
	    throw std::domain_error("Need three parameters for a bottleneck");
	std::vector<double>x=stringtodouble(tok,1);
	return new GenTL::bottleneck(x);
    } else if (tok[0]=="piecewise") {
      if (tok.size()%2!=0) 
	  throw std::domain_error("Need an odd number of parameters for a piecewise model");
      std::vector<double>x=stringtodouble(tok,1);
      return new GenTL::piecewise(x);
    }  else {
      std::ostringstream o;
      o << "'" << tok[0] << "' unknown growth type in gmread" << std::endl;
      throw std::domain_error(o.str().c_str());
    }
} 

std::ostream &operator<<(std::ostream &o, const  GenTL::growthmodel &g)
{
    g.print(o);
    return o;
}
/** what is the next coalescent                              */
std::pair<double,int> 
 GenTL::growthmodel::nextcoalv(const std::vector<int> &n
		       ,const double t,rng &r)
{
  int which=-1;
  double mn=1E10;

  for (size_t i=0;i<n.size();i++) {
    if (n[i]>0) {
      double cl=nextcoal(n[i],t,r);
      if (cl<mn) {
	which=i;
	mn=cl;
      }   
    } 
  }
  return std::pair<double,int>(mn,which);
}
/**                               */
std::pair<double,int> 
GenTL::growthmodel::nextcoalvrel(const std::vector<int> &n, const std::vector<double> rel
		       ,const double t,rng &r)
{
  int which=-1;
  double mn=1E10;

  for (size_t i=0;i<n.size();i++) {
    if (n[i]>0) {
      double cl=nextcoalrel(n[i],t,rel[i],r);
      if (cl<mn) {
	which=i;
	mn=cl;
      }   
    } 
  }
  return std::pair<double,int>(mn,which);
}
