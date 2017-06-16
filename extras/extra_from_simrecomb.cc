

//  a data structure for holding the ARGs
class ARGlist {
public:
  ARGlist() {
    //std::cout << "Initialising the list of ARGs" << std::endl;
    present=0;
  }
  ~ARGlist() {
    //  std::cout << "Destroying the list of ARGs" << std::endl;
    //std::cout << present << " still not destroyed\n";
    //  if (present>0) 
    //   for_each(vsr.begin(),vsr.end(),DeleteObject());
  }
  void push_back(simrecomb *sr) {
    vsr.push_back(sr);
    present +=1;
  }
  size_t size() const {
    return vsr.size();
  }
  size_t nonnull() const {
    return present;
  }
  size_t sites(int wh) const {
    return vsr[wh]->tr.root.size();
  }
  simrecomb *operator[](int wh) {
    return vsr[wh];
  }
  void destroy(int i) {
    if (i<0||static_cast<size_t>(i)>=vsr.size()) {
      Rprintf("warning, %d is not in range (0,%d)\n",i,vsr.size()-1);
      return;
    }
    if (vsr[i]==0) {
      Rprintf("warning %d has already been deleted",i);
      return;
    }
    delete vsr[i];
    present--;
  }
  
private:
  std::deque<simrecomb *> vsr;
  int present;
};

//
static ARGlist ARG;


extern "C" {
  void simplesim(int *ss, int *sites, double *rec, char **growthmod, int *count) {
    try {
      std::string gm(*growthmod);
      rng r;
      std::vector<double> rr(*sites-1);
      for (int i=0;i<*sites-1;i++) rr[i]=rec[i];
      ARG.push_back(new simrecomb(*ss,*sites,rr,gm,r));
      *count=ARG.size()-1;
    }
    catch (const std::exception &e) {  // return  and signal an error
      Rprintf("exception caught: %s\n",e.what());
      *count=-1;
    }
  }
  // add a migration matrix simulation (migration matrix passed as string)
  void mmsim(int *location, int *ss,int *sites, double *rec
               , char **growthmod, char **migmodel, int *count) {
    try {
      rng r;
      std::string gm(*growthmod);
      std::string mm(*migmodel);
      std::vector<double> rr(rec,rec+(*sites-1));
      std::vector<int> loc(location,location+*ss);
      ARG.push_back(new simrecomb(rr,loc,gm,mm,r));
      *count=ARG.size()-1;
    }
    catch (const std::exception &e) {  // return  and signal an error
      Rprintf("exception caught in mmsim: %s\n",e.what());
      *count=-1;
    }
  }
  void popsim(int *ss, int *sites, double *rec, char **growthmod, char **migmodel
                , char **poptree, int *count) {
    rng r;
    std::string gm(*growthmod);
    std::string mm(*migmodel);
    std::string pt(*poptree);
    
    std::vector<double> rr(rec,rec+(*sites-1));
    
    ARG.push_back(new simrecomb(*ss,*sites,rr,gm,r));
    *count=ARG.size()-1;
  }

  /** Get stepping stone model mutations at positions gives by pos           */
  void STRmutate(int *whichARG, double *theta, int *pos, int *npos, int *data)
  {
    rng r;
    simrecomb::ARGtype &st=ARG[*whichARG]->tr;
    int ss=st.sample.size(); 
    int start=0;
    for (int j=0;j<*npos;j++) {
      ARG[*whichARG]->STRmutate(pos[j],theta[j],r);
      for (int i=0;i<ss;i++) 
        data[start+i]=st.sample[i]->data()[0];
      start += ss;
    }
  }
  /** Get the mutations based on a mutation rate theta at each site                     
  * Assumes that theta is a vector of the same length as the 
  * number of sites.  data should be given enough space to take all the 
  * possible mutations - but it should not use them                     
  */
  void mutatetheta(int *whichARG, double *theta, int *data, int *u, int *var, int *muts) {
    rng r;
    simrecomb::ARGtype &st=ARG[*whichARG]->tr;
    int ss=st.sample.size();
    int sites = st.root.size();
    if (sites != *var) {
      Rprintf("Problem, sites is of incorrect size\nHave %d and passed %d\n"
                ,st.root.size(),*var);
      return;
    }
    std::vector<int> nmuts=ARG[*whichARG]->mutatetheta(theta,r);
    
    *var=0;
    for (size_t i=0;i<nmuts.size();i++) {
      if (nmuts[i]>0) {
        u[*var]=static_cast<int>(i);
        muts[*var]=nmuts[i];
        (*var)++;
      }
    }
    
    int count=0;
    for (int ii=0;ii<ss;ii++) 
      for (int j=0;j<*var;j++) 
        data[count++]=st.sample[ii]->data()[u[j]];
  }
  /** Get the mutations based on a mutation rate theta at each site                     
  * Assumes that theta is a vector of the same length as the 
  * number of sites.  data should be given enough space to take all the 
  * possible mutations - but it should not use them                     
  */
  void mutatethetapos(int *whichARG, int *pos,double *theta, int *data
                        , int *u, int *var, int *muts)
  {
    try {
      rng r;
      simrecomb::ARGtype &st=ARG[*whichARG]->tr;
      for (size_t i=0;i<st.sample.size();i++) {
        st.sample[i]->data().assign(1,0);
      } 
      int ss=st.sample.size();
      int mutcount=0;
      int count=0;
      for (int i=0;i<*var;i++) {
        
        int nmuts=ARG[*whichARG]->mutatetheta(theta[i],*pos,r);
        
        if (nmuts>0) {
          u[mutcount]=i+1;
          muts[mutcount++]=nmuts;
          for (int j=0;j<ss;j++) 
            data[count++]=st.sample[j]->data()[0];
        }
      }
      *var=mutcount;
      return;
    }
    catch (const std::exception &e) {  // return  and signal an error
      Rprintf("exception caught: %s\n",e.what());
      return;
    }
  }
  /** Get the mutations above the nodes in the ARG                 */
  void mutateAboveNodes(int *whichARG, int *pos, int *data) 
  {
    try {
      // int ss=ARG[*whichARG]->tr.root[*pos]->nleaves(*pos);
      int ss=ARG[*whichARG]->tr.sample.size();
      ARG[*whichARG]->mutateAboveNode(*pos);
      
      if (*pos<0 or *pos>=static_cast<int>(ARG[*whichARG]->tr.root.size())) {
        Rprintf("Problem, pos must be between 0 and sites-1 in C++, have %d\n",*pos);
        return;
      }
      
      int count=0;
      for (int i=0;i<ss;i++) 
        for (int j=0;j<ss-3;j++) 
          data[count++]=ARG[*whichARG]->tr.data(i)[j];
      
      return;
    }
    catch (const std::exception &e) {  // return  and signal an error
      Rprintf("exception caught: %s\n",e.what());
      return;
    }
  }
  /** extract the lengths that  prior to the first coalescence                     
  lengths is an array of length 2*ss                                       */
  void haplengths(int *which, int *position, double *lengths) {
    try {     
      simrecomb::ARGtype &st=ARG[*which]->tr;
      int ss=st.sample.size();
      int sites=st.root.size();
      double *mnp=lengths;
      double *mxp=lengths+ss;
      
      for (int ii=0;ii<ss;ii++) {
        std::set<GenTL::position> pp=st.sample[ii]->recombinations();
        pp.insert(-1);
        
        std::set<GenTL::position>::iterator iii=pp.lower_bound(*position);
        iii--;
        if (*iii==-1)  mnp[ii]=1;
        else mnp[ii] = (*iii)+2;
        GenTL::position mp=*(pp.upper_bound(*position));
        if (mp==GenTL::maxloci)  mxp[ii]=sites;
        else mxp[ii] =  mp+1;
      }
      return;
    }
    catch (const std::exception &e) {  // return  and signal an error
      Rprintf("exception caught: %s\n",e.what());
      return;
    }
  }
  /** extract the lengths that  prior to the first coalescence                     
  lengths is an array of length 2*ss                                       */
  /*    void sharedsection(int *whichARG, int *samples,int *nsamp, int *position, double *lengths) {
  simrecomb::ARGtype &st=ARG[*whichARG]->tr;
  //   int ss=st.sample.size();
  //  int sites=st.root.size();
  std::vector<int> a(samples,samples+(*nsamp));
  
  TNT::Array2D<int> res=st.sharedsections(a,*position);
  int len=res.dim1();
  for (int i=0;i<len;i++) {
  lengths[i]=res[i][0];
  lengths[i+len]=res[i][1];
  }
  return;
}*/
  /** extract the lengths that  prior to the first coalescence                     
  lengths is an array of length 2*ss                                       */
  /*   void sharedsectionK(int *whichARG, int *samples,int *nsamp, int *position, double *lengths, int *depth) {
  simrecomb::ARGtype &st=ARG[*whichARG]->tr;
  //      int ss=st.sample.size();
  //int sites=st.root.size();
  std::vector<int> a(samples,samples+(*nsamp));
  
  TNT::Array2D<int> res=st.sharedsectionK(a,*position,*depth);
  int len=res.dim1();
  for (int i=0;i<len;i++) {
  lengths[i]=res[i][0];
  lengths[i+len]=res[i][1];
  }
  return;
  }*/
  /** extract the tree at a given position                                  */
  void extracttree(int *whichARG, int *position, unsigned char *tree, int *nchars) {
    simrecomb::ARGtype &st=ARG[*whichARG]->tr;
    
    GenTL::node<simrecomb::nodetype *> 
      *extractroot=extract_node(st.root[*position]
      ,*position
                                  ,st.sample);
    
    rotate_tree_labels(extractroot);
    std::ostringstream oss;      
    extractroot->recursiveprintnodata(oss,true,false);
    
    std::string os=oss.str();
    if (static_cast<int>(os.size())>=*nchars) {
      Rprintf("warning -- need a longer string in extracttree\n");
      return;
    }
    for (size_t j=0;j<os.size();j++) tree[j]=os[j];
    tree[os.size()]=';';
    recursively_destroy_node(extractroot);
  }
  /** prune an ancestral recombination graph                                  */
  void pruneARG(int *which, int *samples, int *nsamps) {
    simrecomb::ARGtype &st=ARG[*which]->tr;
    
    std::vector<size_t> wh(samples,samples + (*nsamps));
    st.prune(wh);
    
  }
  /** See if we can get some haplogroups at a position   */
  void haplogroups(int *whichARG, int *position, double *minsplit, int *hg) {
    simrecomb::ARGtype &st=ARG[*whichARG]->tr;
    std::vector<char> mdef(1);
    GenTL::node<std::vector<char> > 
      *extractroot=extract_node(st.root[*position]
      ,*position
                                  ,st.sample,mdef);
    
    //  std::vector<int> hg=t.haplogroups(0,minsplit);
    /** This is unfinished - I need to be able to get 
     * haplotypes from the rot of a tree - rather than 
     * a tree class             */ 
    
    
  }

}

