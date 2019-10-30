#include <Rcpp.h>
#include "ijw_rand.h"
#include "simrecomb.h"



// [[Rcpp::export]]
SEXP simplesim_rcpp(int ss, int sites, double rec, const std::string growthmod) {
  try {
    std::string gm(growthmod);
    rng r;
    std::vector<double> rr(sites-1);
    for (int i=0; i<sites-1; i++) rr[i]=rec;
    
    simrecomb *s = new simrecomb(ss, sites, rr, gm, r);
    Rcpp::XPtr< simrecomb > pt(s, true);          // get pointer as SEXP
    return pt;
  }
  catch (const std::exception &e) {  // return  and signal an error
    Rprintf("exception caught: %s\n", e.what());
  }
}


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


/*** R
library(ARG)
a <- simplesim_rcpp(500, 10000, 0.01, "constant")
b <- simplesim_rcpp(500, 10000, 0.01, "exponential(1)")

plot(treeheight(a), type="l", xlab="position", ylab="Tree Height", axes=FALSE, ylim=c(0,10))
lines(treeheight(b),col="blue");

axis(1)
axis(2)
rm(a)
rm(b)
gc()

*/
