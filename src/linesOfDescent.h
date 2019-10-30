#ifndef LINES_OF_DESCENT_H__
#define LINES_OF_DESCENT_H__

#include <deque>

namespace GenTL {

  template <class T> class recombnode;
  template <class T> class rnodet;
  /** a simple class for holding the lines of descent  */ 
  template<class T> 
   class lines_of_descent: public std::deque<recombnode<T> *> {
      public:
      /** construct from the sample */
      lines_of_descent( std::vector<rnode<T> *> &sample)
	: std::deque<recombnode<T> *>(sample.size()) {
	for (size_t i=0;i<sample.size();i++) {
	  this->at(i)=sample[i]->first();
	}
      }
      /** add a recombnode         */
      void add(recombnode<T> *a) {  
	push_back(a);
      } 
      /** coalescence -- rearrange the lines-of-descent */
      int coal(rnode<T> *nw,  std::pair<int, int> &wh) {
	this->at(wh.first) = nw->l.front();
	if (this->at(wh.first)->active.empty()) {
	  assert(wh.first>wh.second);
	  this->at(wh.first)=this->back();
	  this->pop_back();
	  this->at(wh.second)=this->back();
	  this->pop_back();
	  return 2;
	}
	this->at(wh.second)= this->back();
	this->pop_back();
	return 1;
      }
      /** print out the current lines                   */
      void print(std::ostream &o) {
	for (size_t i=0;i<this->size();i++) {
	  std::ostringstream s;
	  s <<  i<<": ";
	  this->at(i)->print(o,s.str().c_str());
	}
      }
    private:
      lines_of_descent();
    };
}

#endif
