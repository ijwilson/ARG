//#include "GenTL/coalescent.h"
//#include "GenTL/nodet.h"
#include "coalescent.h"
#include "nodet.h"
#include "ijw_rand.h"
using namespace GenTL; 

/******************************************************************************************/
/******************************************************************************************/
/** Pick the migration event             */
/* returns a pair that is the line which moves
 * and the location to which it is going to move */
  std::pair<int,int> structured_coalescent::PickMigration() {
    int frompop=gen_from_p<std::vector<double> >(rates_out,lrand);
#ifdef CHECK
    std::cout << "rates out ";
    printvector(std::cout,rates_out);
    std::cout << "Migration from population " << frompop <<"  "
	      << left_[frompop] << " left "<< std::endl;
#endif
    int from = lrand.rint(left_[frompop]);
    int to;
    
    for (;;) {  // fix to get "positive" populations
      to=mm_.sample(frompop,lrand);
      if (relsize_[to]>0.0) break;
    }
#ifdef CHECK
	std::cout << " to population " << to;
	std::cout << " individual " << from << std::endl;
#endif
	for (size_t i=0;i<where_.size();i++) {
	  if (where_[i]==frompop) {
	    if (from==0) 
	      return  std::pair<int,int>(i,to);
	    from--;
	  }
	}
	assert(false);
	return std::pair<int,int>(0,0);
      }

double structured_coalescent::next(const double &time) {
#ifdef CHECK
  for (size_t i=0;i<rates_out.size();i++) {
    if (rates_out[i]<0.0) {
      std::ostringstream oss;
      oss << "thrown as rates_out[" << i <<"] = " << rates_out[i];
      throw std::range_error(oss.str().c_str());
    }
  }
     
  if (fabs(accumulate(rates_out.begin(),rates_out.end(),0.0)-mig_rate)>1E-5) {
    printvector(std::cerr,rates_out);
    printvector(std::cerr,left_);
    std::cerr << "mig rate = " << mig_rate << std::endl;
  }
  if (mig_rate<0) {
    std::cerr << " left\n";
    printvector(std::cerr,left_,true);
    std::ostringstream oss;
    oss << "thrown as mig = " << mig_rate;
    throw std::range_error(oss.str().c_str());
  }
#endif
  assert(mig_rate>=0.0);
  double mig=1E50;
  if (mig_rate>0.0) 
    mig=lrand.sexp()/mig_rate;

  std::pair<double,int> coal=g->nextcoalvrel(left_,relsize_,time,lrand);
#ifdef CHECK
  std::cout << "Coalescence time "<< coal.first << " in population "
	    << coal.second << " migration time " << mig << std::endl;
#endif
  if (mig<coal.first) {
    whichcoal=-1;
    return mig;
  } else {
    whichcoal=coal.second;
    return coal.first;
  }
};
/** ------------------ which lines coalesce?  ------------------    */
std::pair<int, int> structured_coalescent::coalchromo() {
#ifdef CHECK
	std::cout << " in coalchromo - whichcoal = "
		  << whichcoal << std::endl;
#endif
	if (whichcoal<0) {
	  std::pair<int,int> p= PickMigration();
	  migrate(p);
	  return std::pair<int,int>(0,0);
	} else {
	  return PickCoalescence();
	}
      };
/** ------------------      update after a coalescent   ------------------       */
void structured_coalescent::UpdateCoalescence(std::pair<int,int> &lines, int remove) {
	assert(remove==1||remove==2);
	assert(lines.first>lines.second);
#ifdef CHECK
	printvector(std::cout,where_,true);
	std::cout << "Before fixing, rates_out\n";
	printvector(std::cout,rates_out,true);
	std::cout << "left = " << std::endl;
	printvector(std::cout,left_,true);
	std::cout << "Remove = " << remove << std::endl;
#endif
	assert(where_[lines.first]==where_[lines.second]);
		  
	if (remove==2) {
	  left_[where_[lines.first]] -=1;	
	  double change = mm_.RateOut(where_[lines.first]);
	  //	  if (fabs(rates_out[where_[lines.first]])<1E-10) change=0.0;
	  rates_out[where_[lines.first]] -= change;
	  mig_rate -= change;
	  where_[lines.first]=where_.back();
	  where_.pop_back();
	} 
	left_[where_[lines.second]] -=1;
	double change = mm_.RateOut(where_[lines.second]);
	//	if (fabs(rates_out[where_[lines.first] ])<1E-10) change=0.0;
	rates_out[where_[lines.second]] -= change;
	mig_rate -= change;
	where_[lines.second]=where_.back();
	where_.pop_back();
	nleft_-=remove;
#ifdef CHECK
	std::cout << "++++++++++++++++++++++++++++++++++++++++\n" 
		  << "After change, where\n";
	printvector(std::cout,where_,true);
	std::cout << "After fixing, rates_out\n";
	printvector(std::cout,rates_out,true);
	std::cout << "left = " << std::endl;
	printvector(std::cout,left_,true);
	std::cout << "Remove = " << remove << std::endl;
#endif	
      }
/**   print the data for debugging                              */
std::ostream &structured_coalescent::print(std::ostream &o) {
  subdivcoal::print(o);
  o << "Structure ";
  mm_.print(o);
  o << std::endl;
  o << "Rates out" << std::endl;
  printvector(o,rates_out);
  return o;
};
/** update the recombination rates                                         */
void structured_coalescent::UpdateRecombination(int recomb) {
  nleft_ +=1;
  int location=where_[recomb];
  where_.push_back(location);
  left_[location] +=1;
  double change=mm_.RateOut(location);
  rates_out[location] +=change;
  mig_rate += change;
};
/** perform the migration                                                 */
void  structured_coalescent::migrate(const std::pair<int,int> &lineto) {
#ifdef CHECK
  std::cout << "Moving line " << lineto.first << " to "
	    << "population " << lineto.second << std::endl;
#endif 
  assert(lineto.first>=0 && lineto.first < nleft_);
  assert(lineto.second < mm_.npops() && lineto.second >=0);
  
  int remove=where_[lineto.first];
  assert(remove >=0 && remove < mm_.npops());
  left_[remove] -=1;
  left_[lineto.second]+=1;
  rates_out[remove] -= mm_.RateOut(remove);
  rates_out[lineto.second] += mm_.RateOut(lineto.second);
  where_[lineto.first]=lineto.second;
  mig_rate = std::accumulate(rates_out.begin(),rates_out.end(),0.0);
};
/** -------------------- constructor ----------------------------------  */
structured_coalescent::structured_coalescent(const std::vector<int> &where
					     ,growthmodel *gm, mig_matrix &m
					     ,rng &r):
  subdivcoal(where,std::vector<double>(m.npops(),1.0),gm,r)
  ,mm_(m),mig_rate(0.0),rates_out(npops_,0.0) {    
  for (size_t j=0;j<npops_;j++) {
    rates_out[j]  = double(left_[j])*mm_.RateOut(j);
    mig_rate += rates_out[j];
  }
};
/********************************************************************************************/
/********************************************************************************************/
splitting_coalescent::splitting_coalescent(const std::vector<int> &where
					   ,growthmodel *gm
					   , GenTL::tree<double> *ptree,rng &r):
  subdivcoal(where,std::vector<double>(ptree->nleaves(),0.0),gm,r),pt_(ptree)
{    
  std::vector<double> tm(ptree->nleaves()-1);
  
  if (int(ptree->nleaves()-1)!= *(std::max_element(where.begin(),where.end())))
    throw std::range_error("should have equal leaves and max pop label");

  if (ptree->nleaves()!= npops_)
    throw std::range_error("should have equal leaves and npops_");

  for (unsigned int i=0;i<npops_-1;i++) {
    std::cerr << i << " " << pt_->anctime(i) << " " << pt_->data(i) << std::endl;
    tm[i]=pt_->anctime(i);
    relsize_[i]=pt_->data(i);
  }
  relsize_[npops_-1]=pt_->data(npops_-1);
  
  indx.sort(tm);

#ifdef CHECK
  printvector(std::cerr,tm);
  for (size_t i=0;i<npops_-1;i++) std::cerr << indx[i]<< " ";
  std::cerr << std::endl;
#endif
  nextjoin=0;
};
/********************************************************************************************/
void splitting_coalescent::fixpoptree() 
{ 
    const node<double> *jnode(pt_->ancestorptr(indx[nextjoin]));
    int keep,move;
    if (jnode->label==jnode->left->label) {
	keep=jnode->left->label;
	move=jnode->right->label;
    } else if (jnode->label==jnode->right->label) {
	keep=jnode->right->label;
	move=jnode->left->label;
    } else 
	throw std::domain_error("error - internal node does not specify "
			    "which node to keep in fixpoptree");
  
  relsize_[keep]=jnode->data();
  relsize_[move]=0.0;

  for (int i=0;i<nleft_;i++) {
    if (where_[i]==move) where_[i]=keep;
  }
  left_[keep]+=left_[move];
  left_[move]=0;
  
  if (jnode==pt_->root) nextjoin=-99; // is this the root
  else nextjoin+=1;

 }


/** ----------------- which lines coalesce? --------------------------- */
std::pair<int, int>  splitting_coalescent::coalchromo() {
  if (whichcoal<0) { // no coalescence, just a fusion event
    fixpoptree();
    return std::pair<int,int>(0,0);
  } else {
    return PickCoalescence();
  }
};

/** ------------------- Pick the coalescence event --------------------  */
std::pair<int,int> subdivcoal::PickCoalescence()  {
  std::pair <int,int> pr= lrand.sample2intsorted(0,left_[whichcoal]-1);
#ifdef CHECK
  std::cout << "sampling coalescence between 0 and " 
	    << left_[whichcoal]-1 << std::endl;
  std::cout << "Using lines " << pr.second <<" and " 
	    << pr.first << std::endl;
#endif
  pr.first-=(pr.second+1);
  unsigned int i,j;
  for (i=0;i<where_.size();i++) {
    if (where_[i]==whichcoal) {
      if (pr.second==0) break;
      pr.second--;
    }
  } 
  pr.second=i;
#ifdef CHECK
  std::cout << "First is at position "<< i << std::endl;
#endif
  for (j=i+1;j<where_.size();j++) {
    if (where_[j]==whichcoal) {
      if (pr.first==0) break;
      pr.first--;
    }      
  }
  pr.first=j;
#ifdef CHECK
  std::cout << "and second at position "<< j << std::endl;
#endif
  return pr;
}

/** update after a recombination   */
void subdivcoal::UpdateRecombination(int recomb) {
  nleft_ +=1;
  int location=where_[recomb];
  where_.push_back(location);
  left_[location] +=1;
}
/** update after a coalescence   */
void subdivcoal::UpdateCoalescence(std::pair<int,int> &lines, int remove) {
  assert(remove==1||remove==2);
  assert(lines.first>lines.second);
#ifdef CHECK
  std::vector<int>::iterator i=where_.begin();
  while (i != where_.end()) {
    std::cout << *i++ << " ";
  }
  std::cout << std::endl;
#endif
  assert(where_[lines.first]==where_[lines.second]);	  
  if (remove==2) {
    left_[where_[lines.first]] -=1;
    where_[lines.first]=where_.back();
    where_.pop_back();
  } 
  left_[where_[lines.second]] -=1;
  where_[lines.second]=where_.back();
  where_.pop_back();
  nleft_ -=remove;
}
std::ostream &subdivcoal::print(std::ostream &o) {
  o << "Growth Model ";
  g->print(o);
  o << std::endl;
  o << nleft_ << " Individuals in " << npops_ <<" Populations with sizes: ";
  copy(left_.begin(),left_.end(),std::ostream_iterator<int>(o," "));
  o << std::endl;
  return o;
};

/** return the next time, and if it is a coalescence  */
double splitting_coalescent::next(const double &time) {
  std::pair<double,int> coal=g->nextcoalvrel(left_,relsize_,time,lrand);
#ifdef CHECK
  std::cout << "Relative sizes: ";
  printvector(std::cout,relsize_);
  std::cout << "Lines Left: ";
  printvector(std::cout,left_);
  std::cout << "Coalescence time "<< coal.first << " in population "
	    << coal.second << std::endl;
#endif
  if (nextjoin>=0) {
    if( time+coal.first>pt_->anctime(indx[nextjoin])) {
#ifdef CHECK
      std::cout << "joining populations at time " 
		<< pt_->anctime(indx[nextjoin]) << std::endl;
#endif
      whichcoal=-1;
      return pt_->anctime(indx[nextjoin])-time;//get time
    }
  }
  whichcoal=coal.second;
  return coal.first;
}
std::ostream &splitting_coalescent::print(std::ostream &o) {
  subdivcoal::print(o);
  o << pt_ << std::endl;
  if (nextjoin>=0) {
    const node<double> *jnode(pt_->ancestorptr(indx[nextjoin]));
    int keep=jnode->left->label;
    int move=jnode->right->label;
    o<< "Next population join between " << keep << " and " << move
     << " at time " << pt_->anctime(indx[nextjoin]) << " nextjoin: " 
     << nextjoin << " index " << indx[nextjoin];
  }  else o << "All populations have joined";
  o << std::endl;
  return o;
};
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
double  splitting_structured_coalescent::next(const double &time) {
 if (nextjoin>=0) { 
   double mt=structured_coalescent::next(time);
   if( time+mt>pt_->anctime(indx[nextjoin])) {
#ifdef CHECK
     std::cout << "joining populations at time " 
	       << pt_->anctime(indx[nextjoin]) << std::endl;
     std::cout << "nextjoin is nextjoin" << std::endl;
#endif
     whichcoal=-2;
     return pt_->anctime(indx[nextjoin])-time;
   }
   return mt;
 } else {
   whichcoal=where_[0]; // must be the only location left
   return g->nextcoalrel(left_[whichcoal],time,relsize_[whichcoal],lrand);
 }
}
/** print                                                              */
std::ostream &splitting_structured_coalescent::print(std::ostream &o) {
  structured_coalescent::print(o);
  o << *pt_ << std::endl;
  if (nextjoin>=0) {
    const node<double> *jnode(pt_->ancestorptr(indx[nextjoin]));
    int keep=jnode->left->label;
    int move=jnode->right->label;
    o<< "Next population join between " << keep << " and " << move
     << " at time " << pt_->anctime(indx[nextjoin]) << " nextjoin: " 
     << nextjoin << " index " << indx[nextjoin];
  }  else o << "All populations have joined";
  o << std::endl;
  return o;
}
/********************************************************************************************/
/** which lines coalesce?  */
std::pair<int, int> splitting_structured_coalescent::coalchromo() {
#ifdef CHECK
  std::cout << " in splitting_coalescent::coalchromo - whichcoal = "
	    << whichcoal << std::endl;
#endif
  if (whichcoal==-1) {
    std::pair<int,int> p= PickMigration();
    migrate(p);
    return std::pair<int,int>(0,0);
  } else if (whichcoal==-2) {
    fixpoptree();
#ifdef CHECK
    std::cout << "fixed poptree" << std::endl;
#endif
    return std::pair<int,int>(-1,-1);
  } else {
    return PickCoalescence();
  }
};
/********************************************************************************************/
void splitting_structured_coalescent::fixpoptree() 
{
 const  node<double> *jnode(pt_->ancestorptr(indx[nextjoin]));
  int keep,move;
  if (jnode->label==jnode->left->label) {
    keep=jnode->left->label;
    move=jnode->right->label;
  } else if (jnode->label==jnode->right->label) {
    keep=jnode->right->label;
    move=jnode->left->label;
  } else 
    throw std::domain_error("error - internal node does not specify "
			    "which node to keep in fixpoptree");
  relsize_[keep] = jnode->data();
  relsize_[move]=0.0;

  for (int i=0;i<nleft_;i++) {
    if (where_[i]==move) where_[i]=keep;
  }
  left_[keep]+=left_[move];
  left_[move]=0;

  mig_rate=0.0;

  if (jnode==pt_->root) {
    nextjoin=-99; // is this the root
    for (size_t j=0;j<npops_;j++) rates_out[move]  = 0.0;
  } else {
    nextjoin+=1;
    gmm.remove_pop(move);
    for (size_t j=0;j<npops_;j++) {
      rates_out[j]=double(left_[j])*mm_.RateOut(j);
      mig_rate += rates_out[j];
    }
  }
 
}

/** constructor                                                          */
splitting_structured_coalescent::splitting_structured_coalescent(const std::vector<int> &where
								 ,growthmodel *gm,  general_mig_matrix &m
								 ,GenTL::tree<double> *ptree,rng &r):
  structured_coalescent(where,gm,m,r),gmm(m),pt_(ptree)
{    
  std::vector<double> tm(ptree->nleaves()-1);

  for (size_t i=0;i<npops_-1;i++) {
    tm[i]=pt_->anctime(i);
    relsize_[i]=pt_->data(i);
  }
  relsize_[npops_-1]=pt_->data(npops_-1);
   
  indx.sort(tm);
#ifdef CHECK
  for (size_t i=0;i<tm.size();i++) std::cout << indx[i] <<" "<<tm[indx[i]] << std::endl;
#endif  
  nextjoin=0;        
}

coalescent *create_coalescent(int ss,const std::string &growth_model, 
			      const std::string &mig_model, const std::string &pop_tree, rng &rd)
{
  
    tree<double> pt;
      
    try {
      growthmodel *gm=gmread(growth_model);
      mig_matrix *m=mmread(mig_model);
      if (pop_tree!="") {
	std::istringstream iss(pop_tree);
	iss >> pt;
      }
      
      if (m==0&&pt.nleaves()==0) {
	return new GenTL::coalescent(ss,gm,rd);
      } else if (m!=0&&pt.nleaves()>0) {
	GenTL::general_mig_matrix gmm(*m);
	std::vector<int> a;
	a.reserve(ss*m->npops());
	for (int i=0;i<m->npops();i++) 
	  for (int j=0;j<ss;j++) a.push_back(i);
	return new GenTL::splitting_structured_coalescent(a,gm,gmm,&pt,rd);
      }  else if (m!=0) {
	std::vector<int> a;
	a.reserve(m->npops()*ss);
	for (int i=0;i<m->npops();i++) 
	  for (int j=0;j<ss;j++) a.push_back(i);
	return new GenTL::structured_coalescent(a,gm,*m,rd);
      } else {
	std::vector<int> a;
	a.reserve(ss*pt.nleaves());
      for (size_t i=0;i<pt.nleaves();i++) 
	for (int j=0;j<ss;j++) a.push_back(i);
      return new GenTL::splitting_coalescent(a,gm,&pt,rd);
      } 
    }
    catch (std::exception &e) {
      std::ostringstream oss;
      oss << "error: " << e.what() << "\n";
      Rprintf(oss.str().c_str());

      throw;
    }
    catch (...) {
      Rprintf("unknown exception in create_coalescent");

      throw;
    }
}
