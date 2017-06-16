#include <Rcpp.h>
#include "ijw_rand.h"
#include "simrecomb.h"



// [[Rcpp::export]]
SEXP simplesim(int ss, int sites, double rec, const std::string growthmod) {
  try {
    rng r;
    std::vector<double> rr(sites-1, rec);

    simrecomb *s = new simrecomb(ss, sites, rr, growthmod, r);
    Rcpp::XPtr< simrecomb > pt(s, true);          // get pointer as SEXP
    return pt;
  }
  catch (const std::exception &e) {  // return  and signal an error
    Rprintf("exception caught: %s\n", e.what());
    
  }
  return Rcpp::XPtr< simrecomb >(NULL, true);
}
// add a migration matrix simulation (migration matrix passed as string)
// [[Rcpp::export]]
SEXP mmsim(Rcpp::IntegerVector location,int sites, double rec
             , const std::string growthmod, const std::string migmodel) {
  try {
    rng r;
    std::vector<double> rr(sites-1, rec);
    location = location-1;
    std::vector<int> loc = Rcpp::as<std::vector<int> >(location);
    simrecomb *s = new simrecomb(rr, loc, growthmod, migmodel, r);
    Rcpp::XPtr< simrecomb > pt(s, true);          // get pointer as SEXP
    return pt;
  }
  catch (const std::exception &e) {  // return  and signal an error
    Rprintf("exception caught in mmsim: %s\n",e.what());
  }
}

/**  extract the height of the trees in the ARG          */
// [[Rcpp::export]]
Rcpp::NumericVector treeheight(SEXP ptr) {
    Rcpp::XPtr< simrecomb > s(ptr);
    size_t sites=s->tr.sites();
    Rcpp::NumericVector times(sites);
      
    for (size_t i=0; i<sites; i++) {
        times[i]=s->tr.root[i]->time();
    }
    return times;
}
/** extract the tree lengths from the trees in the ARG   */
// [[Rcpp::export]]
Rcpp::NumericVector treelength(SEXP ptr) {
  Rcpp::XPtr< simrecomb > s(ptr);
  size_t sites=s->tr.sites();
  Rcpp::NumericVector times(sites);
  
  for (size_t i=0;i<sites;i++) times[i]=s->tr.length(i);
  return times;
}

/** Get the mutations at var sites in the ARG          */
// [[Rcpp::export]]
Rcpp::List mutate(SEXP ptr, int var, const std::string ascert="panel(3)") {
  rng r;
  Rcpp::XPtr< simrecomb > s(ptr);
  simrecomb::ARGtype &st=s->tr;
  
  int ss=st.sample.size();
  Rprintf("sample size = %d\n",ss);
  std::vector<int> positions=s->mutate(var, ascert, r);
  
  Rcpp::IntegerMatrix haps(ss, var);
  for (int i=0; i<ss; i++) 
    for (int j=0; j< var; j++) 
      haps(i, j) = st.sample[i]->data()[positions[j]];
  
  return Rcpp::List::create(Rcpp::Named("position")=Rcpp::IntegerVector(positions.begin(), positions.end()) +1, 
                            Rcpp::Named("haplotypes")=haps);
  
}

/** Get the TMRCA between two leaves at position pos                        */
// [[Rcpp::export]]
double  TMRCA(SEXP graph, int position, Rcpp::IntegerVector samps) {
  Rcpp::XPtr< simrecomb > s(graph);
  return s->tr.TMRCA(samps[1]-1, samps[2]-1, position-1);
}

/** prune an ancestral recombination graph                                  */
// [[Rcpp::export]]
void pruneARG(SEXP ptr, Rcpp::IntegerVector samples) {
  Rcpp::XPtr< simrecomb > s(ptr);
  samples = samples -1;   // convert to 0 offset
  std::vector<size_t> wh(samples.begin(), samples.end());
  s->tr.prune(wh);
}


/*** R
library(ARG)
a <- simplesim(500, 20000, 0.01, "constant")
b <- simplesim(500, 20000, 0.01, "exponential(1)")

m <- mmsim(rep(c(1,2),c(250,250)), 20000, 0.01, "constant", "island(2,1)")

ma <- mutate(a, 1000)
mb <- mutate(b, 1000)
mm <- mutate(m, 1000)

opar=par(mfrow=c(2,1))

plot(treeheight(a), type="l", xlab="position", ylab="Tree Height", axes=FALSE, ylim=c(0,10))
lines(treeheight(b),col="blue");

plot(density(colMeans(mb$haplotypes), from=0, to=1), col="blue", main="")
lines(density(colMeans(ma$haplotypes), from=0, to=1), col=1)

axis(1)
axis(2)

par(opar)


TMRCA(a,100, 1,2)
small <- simplesim(100, 5000, 0.005, "constant")
plot(treeheight(small), type="l")
pruneARG(small, 1:10)
lines(treeheight(small), col="red")
rm(a)
rm(b)
gc()

*/
