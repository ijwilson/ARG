/** @file */
#ifndef ASCERTAIN_H__IJW
#define ASCERTAIN_H__IJW

#include <vector>
#include <string>
#include "ARG.H"
#include "ijw_rand.h"
#include "mutation.H"
#include "utilityfunctionals.H" // for cumsum

namespace GenTL {
  /** A class to automate ascertaiment of variable sites from a tree         */
  template<typename T,typename COLLECTOR>
  class locusSelector {
  public:
    locusSelector(ARG<T,COLLECTOR> &t, std::vector<double> &len):
      nloci(t.nloc()),l(len),cuml(nloci),localtree(t) {
      // calculate the cumulative lengths
      std::transform(l.begin(),l.end(),cuml.begin(),cumsum<double>());
    };
    /** select a panel of k loci based on a hit in panels                   */
    std::vector<int> panelhit(int k, const std::vector<int> &panel
                              ,rng &ran,bool repeatmut=true);
   
    /** select a panel of k loci based on a double hit in panels a and b     */
    std::vector<int> doublehit(int k, const std::vector<int> &panela
                               ,const std::vector<int> &panelb
                               ,rng &ran,bool repeatmut=true);
    
    /** Select loci at random (selecting proportional to length of tree)     */
    std::vector<int> select(int k, rng &ran);
    /** Select loci with mutations older than time age                       */
    std::vector<int> older(int k, double age,rng &ran);
  private:
    int nloci;
    std::vector<double> &l;
    std::vector<double> cuml;
    ARG<T,COLLECTOR> &localtree;
  };

  template<typename T,typename COLLECTOR>
  std::vector<int> locusSelector<T,COLLECTOR>::select(int k, rng &ran) {
    std::vector<int> pos;
    pos.reserve(k);
    if (double(k)/nloci >0.1) {
      std::ostringstream oss;
      oss << "only supports up to 10% selection of loci" << std::endl 
          << "we have k = " << k << " with " << nloci << " loci " << std::endl;  
      throw std::domain_error(oss.str().c_str());
    }
    for (int count=0;count<k;) {
      int locus=gen_from_cump(cuml,ran);
      if (std::find(pos.begin(),pos.end(),locus)==pos.end()) {  // not already there
        count++;
        pos.push_back(locus);
        // reset the position
        for (size_t ii=0;ii<localtree.sample.size();ii++) 
          localtree.sample[ii]->data()[locus]=0;	
        mutateINF(localtree.root[locus],locus,ran);
      }
    }
    std::sort(pos.begin(),pos.end());
   
    return pos;
  }
  /** select mutations that are older than age time                        */
  template<typename T,typename COLLECTOR>
  std::vector<int> locusSelector<T,COLLECTOR>::older(int k, double age,rng &ran) {
    std::vector<int> pos;
    pos.reserve(k);
    if (double(k)/nloci >0.1) {
      std::ostringstream oss;
      oss << "only supports up to 10% selection of loci" << std::endl 
          << "we have k = " << k << " with " << nloci << " loci " << std::endl;  
      throw std::domain_error(oss.str().c_str());
    }
    for (int count=0;count<k;) {
      int locus=gen_from_cump(cuml,ran);
      if (std::find(pos.begin(),pos.end(),locus)==pos.end()) {  // not already there
        for (size_t ii=0;ii<localtree.sample.size();ii++) localtree.sample[ii]->data()[locus]=0;
        double muttime = mutateINF(localtree.root[locus],locus,ran);
        if (muttime>age) {  // keep
          count++;
          pos.push_back(locus);
        }
      }
    }
    std::sort(pos.begin(),pos.end());
    return pos;
  }
  /** select a panel of k loci based on a hit in panel  */
  template<typename T,typename COLLECTOR>
  std::vector<int> locusSelector<T,COLLECTOR>::panelhit(int k, const std::vector<int> &panel, 
                                                        rng &ran, bool repeatmut) {
    std::vector<int> pos;
    pos.reserve(k);
    if (double(k)/nloci >0.1) {
      std::ostringstream oss;
      oss << "only supports up to 10% selection of loci" << std::endl 
          << "we have k = " << k << " with " << nloci << " loci " << std::endl;  
      throw std::domain_error(oss.str().c_str());
    }
    for (int count=0;count<k;) {
      int locus=gen_from_cump(cuml,ran);
      if (std::find(pos.begin(),pos.end(),locus)==pos.end()) {  // not already there
        // reset values to zero
        for (size_t ii=0;ii<localtree.sample.size();ii++) 
          localtree.sample[ii]->data()[locus]=0;	
        // mutate at position locus
        mutateINF(localtree.root[locus],locus,ran);
        // what is the value in the first member of the panel
        int val=localtree.sample[panel[0]]->data()[locus];
        bool uselocus=false;
        for (size_t i=1;i<panel.size();i++) {
          if (localtree.sample[panel[i]]->data()[locus]!=val) {
            uselocus=true;
            break;
          }
        }
        if (uselocus) {
          count++;
          pos.push_back(locus);
        } else if (!repeatmut) count++;
      }
    }
    std::sort(pos.begin(),pos.end());
    
    return pos;
  }
  /** select a panel of k loci based on a double hit in panel a and b */
  template<typename T,typename COLLECTOR>
  std::vector<int> locusSelector<T,COLLECTOR>::doublehit(int k, const std::vector<int> &panela, 
                                                         const std::vector<int> &panelb, 
                                                         rng &ran,bool repeatmut) {
    std::vector<int> pos;
    pos.reserve(k);
    if (double(k)/nloci >0.1) {
      std::ostringstream oss;
      oss << "only supports up to 10% selection of loci" << std::endl 
          << "we have k = " << k << " with " << nloci << " loci " << std::endl;  
      throw std::domain_error(oss.str().c_str());
    }
    for (int count=0;count<k;) {
      int locus=gen_from_cump(cuml,ran);
      if (std::find(pos.begin(),pos.end(),locus)==pos.end()) {  // not already there
        for (size_t ii=0;ii<localtree.sample.size();ii++) 
          localtree.sample[ii]->data()[locus]=0;	
        mutateINF(localtree.root[locus],locus,ran);
        int vala=localtree.sample[panela[0]]->data()[locus];
        bool uselocus=false;
        for (size_t i=1;i<panela.size();i++) {
          if (localtree.sample[panela[i]]->data()[locus]!=vala) {
            uselocus=true;
            break;
          }
        }

        if (uselocus) {  // in the first panel
          // std::cout << "in first panel" << std::endl;
          int valb=localtree.sample[panelb[0]]->data()[locus];
          uselocus=false; 
          for (size_t i=1;i<panelb.size();i++) {
            if (localtree.sample[panelb[i]]->data()[locus]!=valb) {
              uselocus=true;
              break;
            }
          }
        }
        if (uselocus) {
          count++;
          pos.push_back(locus);
        } else if (!repeatmut) count++;
      }
    }
    std::sort(pos.begin(),pos.end());
    
    return pos;
  }

  /// ascertained at random
  template <typename T,typename COLLECTOR>
  std::vector<int> atrandom(ARG<T,COLLECTOR > &t,int k
                            , rng &ran)
  {
    std::vector<int> u=ran.integer_choose(k,t.nloc());

    for (int i=0;i<k;i++) {
      for (size_t ii=0;ii<t.sample.size();ii++) 
        t.sample[ii]->data()[u[i]]=0;
      mutateINF(t.root[u[i]],u[i],ran);   
    }
    return u;
  }


  template<typename T,typename COLLECTOR>
  std::vector<int> ascertainedMutation(std::string ascertainment
                                       , ARG<T ,COLLECTOR> &t
                                       ,int k, std::vector<double> &len , rng &ran
                                       , bool repeatmut=true) 
  {
    std::vector<std::string> va;
    Tokenize(ascertainment,va," (),");
    if (va.size()==0||va[0]=="none") {  // select at random
      return atrandom(t,k,ran);
    } else if (va[0]=="bylength") {
      // std::cout << "by length" << std::endl;
      if (va.size()!=1) 
        throw std::domain_error("bylength takes no parameters");
      GenTL::locusSelector<T,COLLECTOR> sel(t,len);
      return sel.select(k,ran);
    }
    else if (va[0]=="older") {
      // std::cout << "by length" << std::endl;
      if (va.size()!=2) 
        throw std::domain_error("older takes a single double parameter");
      double age =atof(va[1].c_str());
      GenTL::locusSelector<T,COLLECTOR> sel(t,len);
      return sel.older(k,age,ran);
    } else if (va[0]=="panel") {
      if (va.size()!=2) 
        throw std::domain_error("panel takes a single integer parameter");
      int add =atoi(va[1].c_str());
      //std::cout << "Panel size is " << add <<" diploids " << std::endl << "Panel: ";
  
      std::vector<int> panel;
      panel.reserve(2*add);;
      for (int i=0;i<2*add;i++) {
        panel.push_back(i);
      }
      //copy(panel.begin(),panel.end(),std::ostream_iterator<int>(std::cout," "));
      // std::cout << std::endl;
      GenTL::locusSelector<T,COLLECTOR> sel(t,len);
      return sel.panelhit(k,panel,ran,repeatmut);
    } else if (va[0]=="doublepanel") {
      if (va.size()!=4) 
        throw std::domain_error("doublepanel takes three integer parameters");
      int add =atoi(va[1].c_str());
      std::vector<int> panela;
      panela.reserve(2*add);;
      for (int i=0;i<2*add;i++) panela.push_back(i);
      int start =atoi(va[2].c_str());
      add =atoi(va[3].c_str());
      std::vector<int> panelb;
      panelb.reserve(2*add);;
      for (int i=2*start;i<2*start+2*add;i++) panelb.push_back(i);
      //   std::cout << "Panela: ";
      //        copy(panela.begin(),panela.end(),std::ostream_iterator<int>(std::cout," "));
      //        std::cout << std::endl;
      //        std::cout << "Panelb: ";
      //        copy(panelb.begin(),panelb.end(),std::ostream_iterator<int>(std::cout," "));
      //        std::cout << std::endl;
      GenTL::locusSelector<T,COLLECTOR> sel(t,len);
      return sel.doublehit(k,panela,panelb,ran,repeatmut);
    } else {
      std::ostringstream oss;
      oss << ascertainment << " Not recognised as an "
          << "ascertainment method";
      throw std::domain_error(oss.str().c_str());
    }
  }
}



#endif
