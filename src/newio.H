/**     @file                                    */
#ifndef NEWIO_H___
#define NEWIO_H___

#ifdef USE_R
#include <R.h>
#endif

#include <cstdio> // for EOF
#include <cstdlib>

#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <bitset>
#include <limits>
#include <stdexcept>


#ifdef USETNT
#include <tnt/tnt.h>
#endif

using std::istream;
using std::ostream;
using std::string;

/** Several templates to make it easier to join character tests to classes */
template <typename T> 
std::string jointext(const char *t, const T &x, const char *end="") {
  std::ostringstream oss;
  oss << t << " " << x << end;
  return oss.str();
}
template <typename T, typename S> 
std::string jointext(const char *t, const T &x,const char *t2, const S &x2, const char *end="") {
  std::ostringstream oss;
  oss << t << " " << x << " " << t2 << " " << x2 << end;
  return oss.str();
}
template <typename T,typename S, typename U> 
std::string jointext(const char *t, const T &x,const char *t2, const S &x2
                     ,const char *t3, const U &x3, const char *end="") {
  std::ostringstream oss;
  oss << t << " " << x<< " " << t2 << " " << x2 << " " << t3 << " " << x3<< end;
  return oss.str();
}

/** Test if a file can be opened for reading                                */
bool CanReadFile(const char *filename);


template <typename T> void fromstring(const std::string &s, std::vector<T> &a, const char sep=',');
template <class T> std::vector<T> readvector(std::istream &in,T dummy);
template <class T> void printvectorb(std::ostream &o,const std::vector<T> &u,int len=5,int width=14);
std::vector<double> doublevector_scan(std::istream &in, const char *namestring);
template <typename T>
const std::string stemLeafFilename(const T &stem, const char *leaf,  char sep='.');

/** A class to make logging to file of errors and messages easier    
 * Allows you to write both to an output stream (cerr by deafault) and 
 * to a file */
#ifndef USE_R
class logger {
public:
  logger(const char *filename, ostream &screenstream=std::cerr)
    :logfile(filename),screen(screenstream) {
    if (!logfile) {
      std::ostringstream oss;
      oss << "Error opening logfile " << filename << std::endl;
      throw std::runtime_error(oss.str().c_str());
    }
  }
  ~logger() {
    logfile.close();
  }
  void operator()(const std::exception &e) throw() {
     logfile << "error: " << e.what();
     screen  << "error: " << e.what();
#ifdef USE_R
     Rprintf("error %s:", e.what());
#else
     std::cerr << "error: " << e.what();
#endif 
    logfile.flush();
  }
  void operator()(const char *message,bool usescreen=true) {
    logfile << message;
    if (usescreen) screen  << message;
    logfile.flush();
  }
private:
  std::ofstream logfile;
  ostream &screen;
};
#endif
/** 
 * My own exception class - need to work on this           
 */
class ioerror: public std::exception {
 public:
  template<typename T> 
  ioerror(const char *mess, const T &add): message(jointext(mess,add).c_str()) {
#ifdef USE_R
    Rprintf(message.c_str());
#else
    std::cerr << message;
#endif
   
  }

  ioerror(const char *mess)
    :message(mess)  {
#ifdef USE_R
    Rprintf(message.c_str());
#else
    std::cerr << message;
#endif
    
  };

  ~ioerror() throw() {
  }
  
  const char *what() const throw()  {
     return message.c_str();
  }
private:
   const std::string message;
};
/** 
 * My own exception class - need to work on this           
 */
class ijwerror: public std::exception {
 public:
  template<typename T> 
  ijwerror(const char *mess, const T &add): message(jointext(mess,add).c_str()) {
#ifdef USE_R
    Rprintf(message.c_str());
#else
    std::cerr << message;
#endif
  }

  ijwerror(const char *mess)
    :message(mess)  {
#ifdef USE_R
    Rprintf(message.c_str());
#else
    std::cerr << message;
#endif
  };

  ~ijwerror() throw() {
  }
  
  const char *what() const throw()  {
     return message.c_str();
  }
private:
   const std::string message;
};


/**********************************************************************************************/
/**  Read a vector of comma separated Ts from a string   
     cannot return a because otherwise template would not work*/
template <typename T>
void fromstring(const std::string &s, std::vector<T> &a, const char sep) 
{
  std::istringstream iss(s);
  a.resize(0);

  T itmp;
  char ch;
  for (;;) {
    iss >> itmp; 
    if (!iss) break;
    a.push_back(itmp);
    iss >> ch;
    if (!iss) break;
    if (ch != sep) {
      std::ostringstream oss;
      oss << "problem in fromstring,\nExpected separator " << sep 
          << " and got " << ch << std::endl;
      throw ioerror(oss.str().c_str());
    }
  }
}  

template <typename T,typename L>
const std::string stemLeafFilename(const T &stem, const L &leaf)
{
  std::ostringstream oss;
  oss << stem << "." << leaf;
  return oss.str().c_str();
}
template <typename T,typename L,typename R>
const std::string stemLeafFilename(const T &stem, const L &leaf,const R &ex)
{
  std::ostringstream oss;
  oss << stem << "." << leaf << ex;
  return oss.str().c_str();
}

// specialisation for a const char *leaf
template <typename T>
const std::string stemLeafFilename(const T &stem, const char *leaf,  char sep)
{
  std::ostringstream oss;
  oss << stem << sep << leaf;
  return oss.str().c_str();
}

#include "myenums.h"

std::vector <std::string> readstringline(std::istream &in
			       ,int suggestedlength=10, char until='\n');
std::vector<std::vector<std::string> > readstringlines(std::istream &in);
std::vector <int> readintvector(std::istream &in,bool usechar=false, char until=EOF);
std::vector <std::vector<int> > readcharintmat(std::istream &in,int &cols, bool usechar=false);

int skipblank(std::istream &in, char until='\n');
void checkreadcharacter(std::istream &in ,const char wh, const char *message="error");
int skipline(std::istream &in);
void Tokenize(const std::string& str,std::vector<std::string>& tokens
	      ,const std::string& delimiters = " ");
bool skipto(std::istream &in, const char *pattern, const char comment='#');
bool findparametername(std::istream &in, const char *pattern);

std::vector<double> stringtodouble(std::vector<std::string> s,size_t from=0);
//template <class T>
//std::ostream &printvector(std::ostream &o,std::vector<T> &u,bool newline=false,const char *between=" ");
//template <class T>
//std::ostream &printvector(std::ostream &o,const std::vector<T> &u,bool newline=false,const char *between=" ");
//template <class T>
//void printvectorb(std::ostream &o,const std::vector<T> &u,int len,int width);
/*******************
 *
 *
 *
 *******************/
/** Little function to output an STL   */
template <class T>
 std::ostream &printvector(std::ostream &o,const std::vector<T> &u,bool newline=false,const char *between=" ") {
   std::copy(u.begin(),u.end(),std::ostream_iterator<T>(o,between));
   if (newline) o << std::endl;
   return o;
}
// template <class T>
// std::ostream &printvector(const std::vector<T> &u,std::ostream &o, bool newline,const char *between) {
//   std::copy(u.begin(),u.end(),std::ostream_iterator<T>(o,between));
//   if (newline) o << std::endl;
//   return o;
// }
// template <class T>
// std::ostream &printvector(std::ostream &o,const std::vector<T> &u,bool newline,const char *between) {
//   std::copy(u.begin(),u.end(),std::ostream_iterator<T>(o,between));
//   if (newline) o << std::endl;
//   return o;
// }
template <class T>
void printvec(const std::vector<T> &u) {
  std::ostringstream oss;
   std::copy(u.begin(),u.end(),std::ostream_iterator<T>(oss," "));
#ifdef USE_R
   Rprintf(oss.str().c_str());
#else   
   std::cerr << oss.str() << std::endl;
#endif

}



template <class T>
void printvectorb(std::ostream &o,const std::vector<T> &u,int len,int width) {

  for (size_t i=0;i<u.size();i++) {
    if (i%len==0) o << std::endl << std::setw(5) << i << ": ";
    o << std::setw(width) << u[i] << " ";
  }
  o << std::endl;
}
template <class T>
void printvectorfixed(std::ostream &o,const std::vector<T> &u,int precision=5) {
  //  o.setf(0,std::floatfield);
  for (size_t i=0;i<u.size();i++) {
    o   << std::fixed << std::setprecision(precision) << u[i] << " ";
  }
  o << std::endl;
}
#ifdef USETNT
/** read a vector from the TNT::Array format    */
template <class T>
std::vector<T> readvector(std::istream &in,T dummy) {
  TNT::Array1D<T> a;
  in >> a;

  return std::vector<T>(&a[0],&a[0]+a.dim());
}
#endif
template <typename T>
std::vector <std::vector<T> > readmat(std::istream &in,int rows, int cols, T dummy)
{
  std::vector<std::vector<T> > a(rows,std::vector<T>(cols));
  for (int i=0;i<rows;i++) {
    for (int j=0;j<cols;j++) in >> a[i][j];
  }
  return a;
}
/** returns a pointer and a length to a string representation 
 * (provided that T  has a << operator                         */
template <typename T>
std::string tostring(T *start, int length) {
  std::ostringstream oss;
  std::copy(start,start+length,std::ostream_iterator<T>(oss,""));
  return oss.str();
}

double getval(const std::string &str, const char *val,bool loud=false);
std::string getstringfromkeyboard(const char*message, const char *def_val);


int findstart(std::istream &in,const char *pattern);
int skipspace(std::istream &in);
char *readfromquotes(std::istream &in,int *len);

/***************************************************************************/
/*    Now for scanning of old style parameter files                        */
/***************************************************************************/
template <typename T>
T scan(std::istream &in, const char *namestring,const T &default_val, volume vol=quiet)
{
    if (findparametername(in,namestring)) {
      T val;
      in >> val;
      return val;
    } else {
        if (vol==loud) {
          std::ostringstream oss;          
          oss << namestring
		  << " not found in parameter file, using "
		  << default_val << std::endl;
#ifdef USE_R
          Rprintf("%s",oss.str().c_str());
#else
          std::cerr << oss.str();
#endif
        }
      return default_val;
    }
}
/** scan an integer vector after "namestring" from input file          */
template<typename S,typename T>
std::pair<S,T> pair_scan(std::istream &in, const char *namestring,const std::pair<S,T> &def, volume vol=quiet) 
{
  if (findparametername(in,namestring)) {
    S a;
    T b;
    in >> a >> b;
    return std::pair<S,T>(a,b);
  } else {
    if (vol==loud) {
    std::ostringstream oss;
          
    oss  << namestring
         << " not found in parameter file, using (" << def.first << "," << def.second << ")" << std::endl;
#ifdef USE_R
    Rprintf("%s",oss.str().c_str());
#else
    std::cerr << oss.str();
#endif
    }
    return def;
  }
}

template<typename T>
std::vector <T> readvector(std::istream &in, char until='\n') {
  std::vector <T> data;
  T x;
  for (;;) {
    int ch = skipblank(in,until);
    if ((isdigit(ch))||(ch=='-')||ch=='.') {
	in.putback(ch);
	in >> x;
    } else break;
    data.push_back(x);
  }
  return data;
}


std::vector<int> intvector_scan(std::istream &in, const char *namestring,volume vol=quiet);


char *namescan(std::istream &in, const char *namestring, const char *default_val,volume vol=quiet);
int *charintvector_scan(std::istream &in, char *namestring,int *default_val,volume vol=quiet);
std::string stringscan(std::istream &in, const char *namestring, const char *default_val,volume vol=quiet);



//
//
//
template<typename T>
std::vector <T> readvec(std::istream &in, char until, T dummy) 
{
    std::vector <T> data;
    T x;
    for (;;) {
      int ch;
      do {
	ch = in.get();
      } while (ch!=until&&isspace(ch));

      if ((isdigit(ch))||(ch=='.')||(ch=='-')) {
	in.putback(char(ch));
	in >> x;

	data.push_back(x);
      } else  break;
    }
    return data;
}

int nextname(std::istream &in, char *filename);
int get_doublesfrombrackets(std::istream &in,double *par);
bool findparametername(std::istream &in, const char *pattern);
std::vector <double> readdoublevector(std::istream &in, char until='\n');
template <typename T> 
T get_val(std::istream &in) {
  T g;

  char ch=(char)skipspace(in);
  if (ch!='<') return 0.0;
  in >> g;
  return g;
}
template<typename T>
void get_vals(istream &in,T *x, int *count) {
  *count=0;
  char ch=(char)skipspace(in);
  if (ch!='<')
    return;
  for (;;) {
    *count +=1;
    in >> x[*count];
    ch=(char)skipspace(in);
    if (ch=='>')
      break;
    }
    return;
}

int findnonindentedchar(const string &s,int start, char tok=',',char up = '(', char down=')');

void ltrim(std::string &s,char c=' '); 

std::string stringwrap(const std::string &is, const size_t width=72
		       ,const std::string startline="") ;
std::vector <std::string> readLinesUntil(std::istream &in, char FirstChar,bool IgnoreFirst=true);


std::string tempfilename(std::string path, int len);


template<typename T>
std::string xxx_to_bin(const T& value)
{
	const std::bitset<std::numeric_limits<T>::digits + 1> bs(value);
	return std::string(bs.to_string());
}

#endif
