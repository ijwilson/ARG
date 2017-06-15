#include "migmatrix.H"
#include "newio.H"
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

#include "util.H"

// bool IsSquare(int n, int &siz)
// {
//   // take the square root
//   double b=sqrt(double(n));
//   // and round to the nearest integer
//   siz=int(b+0.5);
//   // does this equal the original value?
//   if (siz*siz==n) return true;
//   return false;
// }



/** \brief returns a mig_matrix from a string */


namespace GenTL {
  std::ostream &general_mig_matrix::print(std::ostream &o) const
  {
    o << "General(" << std::endl;
    
    for (size_t i=0;i<m_.size();i++) {  
      o << "{";
      for (size_t j=0;j<m_[0].size()-1;j++) {
	o << m_[i][j]<<",";
      }
      o << m_[i][m_[0].size()-1] << "}" << std::endl;
    }
    return o << ")";   
  }    
}

/** read a migration matrix from a string */
GenTL::mig_matrix *mmread(const std::string &s) {
  if (s=="") return 0;
  std::vector<std::string> tok;
  Tokenize(s,tok,"(), ");
  if (tok[0]=="Island"||tok[0]=="island") {
    if (tok.size()!=3) 
      throw std::domain_error("Need two parameters, k and m for the island model");
    std::vector<double> x=stringtodouble(tok,1);
    return new GenTL::island(size_t(x[0]),x[1]);
  }  
  if (tok[0]=="General"||tok[0]=="general") {
    int siz;
    if (IsSquare(int(tok.size()-1),siz)) {
    } else {
      std::ostringstream os;
      os << "Length of parameters, "<<  tok.size()-1 
	 << " does not give a square matrix in General";
      throw std::domain_error(os.str().c_str());  
    };
    std::vector<double> x=stringtodouble(tok,1);
    return new GenTL::general_mig_matrix(siz,x);
  } else {
      throw std::domain_error("No valid migration model specified");
  }
}


/** Print the mig_matrix type */
std::ostream &operator<<(std::ostream &o, const GenTL::mig_matrix &m)
{
  return m.print(o);
}


GenTL::general_mig_matrix::general_mig_matrix(const GenTL::mig_matrix &mm)
  :m_(mm.npops()) {
  for (int i=0;i<mm.npops();i++) {
    m_[i].resize(mm.npops());
    for (int j=0;j<mm.npops();j++) {
      m_[i][j]=mm(i,j);
    }
  }
}; 
