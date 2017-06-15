#include "newio.H"
#include <stdexcept>
#include <cstring>
#include <ctime>

bool CanReadFile(const char *filename)
{
  std::ifstream inf(filename);
  inf.close();
  if (inf.fail()) return false;
  return true;
}



std::string nextstring(std::istream &in);

/** trim whitespace from thee beginning of a string   */
void ltrim(std::string &s,char c) 
{
  if (s.empty()) return;
 
  std::string::iterator p=s.begin();
  while (*p==c and p!=s.end()) {
    p++;
  }
  s.erase(s.begin(),p);
}

/** Wrap a string at word breaks of length width, and insert line breaks followed 
 the string startline at these points.  Return a string.   Note that this also uses the length of the brk to make 
 sure that evewrything would line up correctly. */

std::string stringwrap(const std::string &is, const size_t width,const std::string startline) 
{
  std::string tmp;
  std::ostringstream oss;
  char cur='\0';
  char last='\0';
  size_t i=0;

  std::string::const_iterator inp=is.begin();

  while (inp!=is.end()) {
    cur=*inp++;
    if (cur=='\n') {
      i=0;
      oss << tmp;
      tmp.clear();
    } else if (++i == width) {
      ltrim(tmp);
      oss << "\n" << startline << tmp;
      i=tmp.length()+startline.length();
      tmp.clear();
    } else if (isspace(cur) && !isspace(last)) { // a word
        oss << tmp;
        tmp.clear();
    }
    tmp += cur;
    last=cur;
  }
  return oss.str();
}

/** read a string from the input stream (ignoring comments) */
std::vector <std::string>
readstringline(std::istream &in,int suggestedlength,char until) {
    std::vector <std::string> data;
    data.reserve(suggestedlength);
    std::string x;
    int ch;

    for (;;) {
        ch = skipblank(in,until);
        if ((isdigit(ch))||isalpha(ch)||(ch=='-')||(ch=='+')||(ch=='(')) {
            in.putback(ch);
            in >> x;
            data.push_back(x);
        }
        else
            break;
    }
    return data;
}
/** read a "matrix" of strings separated by spaces from a file */
std::vector<std::vector<std::string> > readstringlines(std::istream &in) {
    std::vector<std::vector<std::string> > a;
    // check that this is not a commented out line
    int ch = skipblank(in,'\n');
    if (in.eof())
        return a;
    in.putback(ch);
    // read in a line
    std::vector<std::string> x= readstringline(in);
    int len=x.size();
    a.push_back(x);

    for(;;) {
        int ch = skipblank(in,'\n');
        if (ch==EOF)
            break;
        in.putback(ch);
        a.push_back(readstringline(in,len));
    }
    return a;
}
//
//
//
int skipblank(std::istream &in, char until) {
  static int ch;

  for (;;) {
    ch = in.get();
    if (ch=='#')
      ch = skipline(in);
    if (ch==until)
      return (ch);
    if (!isspace(ch))
      return (ch);
  }
}
//
//
//
bool skipto(std::istream &in, const char *pattern, const char comment) {
    int c;
    const char *p = pattern;

    if (in.eof())
        return false;
    while ((c = in.get()) != EOF) {
        if (*p == c)
            p++;
        else {
            if (comment == c) in.ignore(100000,'\n');
            p = pattern;
        }
        if (!*p) {
            return true;
        }
    }
    return false;
}
//
//
//
int skipline(std::istream &in)
/*  a little function to read to the end of a line
    to remove comments !!                            */
{
    char tmp[2000];
    in.getline(tmp,2000);
    if (in.eof()) return EOF;
    return '\n';
}
//
//
void checkreadcharacter(std::istream &in ,const char wh, const char *message) {
    int ch=skipblank(in);
    if (ch!=wh) {
        std::ostringstream oss;
        oss << message << " expected " << wh << ", got " << char(ch);
        throw std::domain_error(oss.str().c_str());
    }
}
//
//
//
std::vector <int> readintvector(std::istream &in,bool usechar, char until) {
  std::vector <int> data;
  int x;

  for (;;) {
    int ch = skipblank(in,until);
    if (ch=='/') 
      ch=skipblank(in,until);
    if ((isdigit(ch))||(ch=='-')) {
	in.putback(ch);
	in >> x;
    } else if (ch=='N'||(usechar&&isalpha(ch))) {
      if (!usechar) {
	int ch2 = in.get();
	if (ch2=='A') {
	  x = -1;
	} else {
	  throw std::domain_error("got 'N' in readintvector");
	}
      } else if (usechar&&isalpha(ch)) {
	x=static_cast<int>(ch);
      }
    } else break;

    data.push_back(x);
  }
  return data;
}

//
//
//
std::vector <std::vector<int> >
readcharintmat(std::istream &in,int &cols, bool usechar)
//  this reads a mixture of integer
{
    std::vector<std::vector<int> > a;
    for (;;) {
      int ch;
      // skip over any coment lines 
      while((ch= skipblank(in,'\n'))=='\n');

      if (in.eof())
	break;
      else
	in.putback(ch);
      std::vector<int> x= readintvector(in,usechar,'\n');
      a.push_back(x);
    }

    if (a.size()>0) {
      cols=a[0].size();
      for (size_t i=1;i<a.size();i++) {
	if (a[i].size()!=size_t(cols)) {
      std::ostringstream oss;
	  oss << "problem in row " << i+1 << std::endl;
	  throw std::length_error(oss.str().c_str());
	}
      }
    }
    return a;
}
// from http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
// used as
// vector<string> tokens;
//    string str("Split me up! Word1 Word2 Word3.");
//    Tokenize(str, tokens);
//
void Tokenize(const std::string& str,
              std::vector<std::string>& tokens,
              const std::string& delimiters) {
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}
//
//
//
std::vector<double> stringtodouble(std::vector<std::string> s,size_t from) {
    if (s.size()<from)
        std::domain_error("vector opf strings too short to convert");
    std::vector<double> x;
    x.reserve(s.size()-from);
    char *pEnd;
    for (size_t i=from;i<s.size();i++) {
        x.push_back(strtod(s[i].c_str(),&pEnd));
        if (*pEnd!='\0') {
            std::ostringstream oss;
            oss << "problem with string " << s[i] <<" in stringtodouble";
            throw std::domain_error(oss.str().c_str());
        }
    }
    return x;
}
/** get a value */

double getval(const std::string &str, const char *val, bool loud) {
    double v=-99.9;
    std::istringstream iss(str);
    if (!findparametername(iss,val)) {
      if (loud) {
#ifdef USE_R
        Rprintf("%g not found in string %s in findparametername\n",val,str.c_str());
#else
        std::cerr << val << " not found in string " << str << " in findparametername" << std::endl;
#endif
      }
    }  else {
        char ch;
        iss >> ch;
        if (ch != '=') {
#ifdef USE_R
         Rprintf("expected = after %s in string\n",val);
#else
          std::cerr << "expected = after " << val << " in string " << str << std::endl;
#endif
        }

        iss >> v;
    }
    return v;
}
/*************************************************************************/
/*  function to get a string from keyboard - skipping
	all blanks in the input - suitable for getting file
	names for example.  If nothing is entered then
	def_val is used as a default                   */
/*************************************************************************/
#ifndef USE_R
std::string getstringfromkeyboard(const char *message,const  char *def_val) {
    std::cerr << "please enter " << message
    << ", or press <enter> for default (" << def_val << ")\n";

    std::string tmp;
    getline(std::cin,tmp);
    if (tmp.empty()) tmp=def_val;
    return tmp;
}
#endif
//
//  search the file for the first
// time that pattern appears, 1 if a success 0 otherwise
//
int findstart(std::istream &in,const char *pattern) {
    int c;
    const char *p = pattern,comment='#';

    if (in.eof())
        in.clear();     // if at end of file clear error
    in.seekg(0);                            // get back to start of file
    while ((c = in.get()) != EOF) {
        if (*p == c)
            p++;
        else {
            if (comment == c)
                skipline(in);
            p = pattern;
        }
        if (!*p)  return 1;
    }
    return 0;
}
//
//  skip spaces including line breaks
//
int skipspace(std::istream &in) {
    static int ch;

    for (;;) {
      ch=in.get();
      if (ch=='#')
	ch = skipline(in);
      if (!isspace(ch))
	return (ch);
    }
}
/****************************************************************/
char *readfromquotes(std::istream &in,int *len) {
    int ch = skipspace(in);
    if (ch!=39) {
      std::ostringstream oss;
      oss << "in readfromquotes: error - should be a quote - ' is " << ch << "\n";
      std::range_error(oss.str().c_str());
    }
    char *data = new char[1024];
    in.getline(data,1024,39);

    *len = 1024; // maybe alter this to realloc
    return data;
}
//
//
//
bool findparametername(std::istream &in, const char *pattern) {
    int c;
    const char *p = pattern,comment='#';

    if (in.eof())
        in.clear();     // if at end of file clear error
    in.seekg(0);				  // get back to start of file
    while ((c = in.get()) != EOF) {
        if (*p == c)
            p++;
        else {
            if (comment == c)
                skipline(in);
            p = pattern;
        }
        if (!*p) {
            if ((c=in.get()) == ':')
                return true;
            p=pattern;
        }
    }
    return false;
}
#ifndef USE_R
/** scan an integer vector after "namestring" from input file          */
std::vector<int> intvector_scan(std::istream &in, const char *namestring,volume vol) {
    if (findparametername(in,namestring)) {
        return readintvector(in);
    } else {
      if (vol!=quiet)
        std::cerr << namestring
		  << " not found in parameter file, using vector of length 0\n";
        return std::vector<int>(0);
    }
}

/** scan an integer vector after "namestring" from input file          */
std::vector<double> doublevector_scan(std::istream &in, const char *namestring) {
    if (findparametername(in,namestring)) {
        return readvector<double>(in);
    } else {
        std::cerr << namestring
        << " not found in parameter file, using vector of length 0\n";
        return std::vector<double>(0);
    }
}
#endif
/***************************************************************************/
char *namescan(std::istream &in, const char *namestring, const char *default_val,volume vol) {
    char *temp=NULL;

    if (findstart(in,namestring)) {
        temp = new char[201];
        if (!temp)
            throw ioerror("error allocating string in namescan\n");
        if (!nextname(in,temp))
            throw ioerror("error reading datafile name\n");
    }
    else if (default_val!=NULL) {
      if (vol!=quiet)
#ifndef USE_R
        std::cerr << namestring << " not found: using default "
		  << default_val << std::endl;
#endif
        temp = new char[201];
        if (!temp) {
          std::ostringstream oss;
          oss << "error allocating string in namescan";
          throw std::range_error(oss.str().c_str());
        }
        strcpy(temp,default_val);
    }
    else {
        temp=NULL;
    }
    return temp;
}
/***************************************************************************/
std::string stringscan(std::istream &in, const char *namestring, const char *default_val,volume vol) {

    if (findstart(in,namestring)) {
        return nextstring(in);
    } else if (default_val!=NULL) {
      if (vol!=quiet) {
#ifndef USE_R
        std::cerr << namestring << " not found: using default "
		  << default_val << std::endl;
#endif
      }
        return std::string(default_val);
    } else {
    	return std::string(""); 
    }
}
/***************************************************************************/
int *charintvector_scan(std::istream &in, char *namestring,int *default_val,volume vol) {
    int *tmp;
    if (findstart(in,namestring)) {
        std::vector<int> a=readintvector(in);
        tmp=new int[a.size()+1];
        tmp[0]=a.size();
        for (size_t i=0;i<a.size();i++)
            tmp[i+1]=a[i];
        return tmp;
    }
    else {
        if (default_val==NULL) {
          if (vol==loud) {
#ifndef USE_R
            std::cerr << namestring << " not found in parameter file, using NULL\n";
#endif
          }     
          tmp=NULL;
        }
        else {
            tmp=new int[2];
            tmp[0]=1;
            tmp[1]=*default_val;
            if (vol==loud) {
#ifndef USR_R
                std::cerr << namestring << " not found in parameter file, using "
                << tmp[1] << std::endl;
#endif
            }
        }
    }
    return tmp;
}
/*********************************************************************/
int nextname(std::istream &in, char *filename)
/*  filename must be at least 255 characters long */
{
    int ch;
    ch= skipspace(in);
    if (ch!=':')
        throw ioerror("expected a :");
    in >> filename;
  

    //    for (i=0;i<255;i++) {
    //        if ((isspace(ch))||ch=='('||ch==EOF) {
    //
    //	  if (ch!=ungetc(ch,in))
    //	    throw "error putting back ch in nextname";
    //            filename[i]='\0';
    //            return 1;
    //        }
    //        filename[i]=(char)ch;
    //        ch = getc(in);
    //    }
    return 1;
}
/*********************************************************************/
std::string nextstring(std::istream &in)
/*  filename must be at least 255 characters long */
{
    int ch;
    ch= skipspace(in);
    if (ch!=':')
        throw ioerror("expected a :");
    std::string filename;
  in >> filename;

    return filename;
}
/*********************************************************************/
int get_doublesfrombrackets(std::istream &in,double *par) {
    char ch;
    //    float f;

    ch=(char)skipspace(in);
    if (ch!='(') {
        par[0]=par[1]=0.0;
        return 0;
    }
    in >> par[0];
    ch=(char)skipspace(in);
    if (ch!=',') {
        par[1]=0.0;
        return 1;
    }
    in >> par[1];

    ch=(char)skipspace(in);
    if (ch==')')
        return 2;
    else {
      std::ostringstream oss;

      oss << "only two parameters can be read in "
          << "getdoublesfrombrackets so far\n ";
      oss << "read (" << par[0] << "," << par[1] << ch
          << "instead of (" << par[0] << "," << par[1] << ")";
      throw std::range_error(oss.str().c_str());
    }
}

int findnonindentedchar(const string &s, int start,char tok,char up, char down)
{
  int endpos=s.length();
  int pos=start;
  int level=0;

  while (pos!=endpos) {
    if (s[pos]==tok&&level==0) break;
    else if (s[pos]==up) level++;
    else if (s[pos]==down) level--;
    pos++;
    if (level<0) throw std::domain_error("level of indentation should never be below zero");
  }
  return pos-start;
}


//
//  read lines from the input until either the end of file or the first character of the 
//  next line is FirstChar
// 
//
// returns a vector of lines
std::vector <std::string> readLinesUntil(std::istream &in, char FirstChar,bool IgnoreFirst) 
{  
  std::string tmp;
  std::vector<std::string> aa(0);

  if (IgnoreFirst) {
    if (!getline(in,tmp)) return aa;
    else aa.push_back(tmp);
  }

  while (in.peek()!=FirstChar && getline(in,tmp)) {
  
    aa.push_back(tmp);
  }
  return aa;
}

std::string tempfilename(std::string path, int len) 
{
  static bool used=false;
  std::string s; 
  if (used==false) {
    //#ifdef WIN32
    //  srand(GetTickCount());
    //#else
    srand( time(NULL));
    //#endif   
    used=true;
  }
  for (int i=0;i<len;i++) {
    int r=rand()%62+48;
    if (r<=57) {
      s.push_back(static_cast<char> (r));
    } else {
      r += 7;
      if (r<=90) {
	s.push_back(static_cast<char> (r));
      } else {
	r += 6;
	s.push_back(static_cast<char> (r));
      }
    }
  }
  return path+"/"+s;
}


