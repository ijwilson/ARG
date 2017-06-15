
#include "util.h"
#include "utilityfunctionals.h"

/************************************************************************/
/*  A function that is needed for WIN32's - where there is no
lgamma defined.  Also needed if you want agreement between
runs on unix machines and on windowsNT - loosely based on the
function in Numerical Recipies									       */
/************************************************************************/
//functionoids
///
lfactorial lfactrl;
//
//
// 

/***************************************************************/
extern "C" {
  void myerror( char *message)
    /* Error handler for my routines     
     * to make sure that old routines that include cutil.h work*/
  {
    Rprintf("Run-time error in ijw's code...\n");
    Rprintf("%s\n",message);
    Rprintf("...I hoped I wouldn't have to do this but...\n");

throw std::runtime_error(message);
  }
}

int pow2(int x)
{
  if (x>30||x<0) error("pow2 only for 0 <= i <=30\n");
  return 1 << x;
}
/** returns the sum of the powers of 2 up to n
    sumpow2(6) = 2+4+8+16+32+64 */
int sumpow2(int x)
{
  if (x>30||x<1) {
    std::ostringstream oss;
    oss << "x = " << x << std::endl;
    oss << "sumpow2 only for 1 <= i <=30\n";
    error(oss.str().c_str());
  }
  return 2*(pow2(x)-1);
}

#ifdef CHECK
double mylog(double x)
{
  if (x<=0.0) {
    warning("problem in mylog: value = ",x,"\n");
  }
  return log(x);
}
#endif

#ifdef WIN32
double lgamma(double xx)
{
	double x,y,tmp,ser;
	static  double cof[6] = {
		76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5	};
	int j;

	y = x = xx;
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser = 1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
#endif


double log_D(const double *b, int n)
 {
   double sum=0.0,temp=0.0;
   for (int i=0;i<n;i++) {
     sum += b[i];
     temp += lgamma(b[i]);
   }
   temp -= lgamma(sum);
   return temp;
 }
//
// Sorting routines
//
void hpsort(unsigned long n, double ra[])
{
	unsigned long i,ir,j,l;
	double rra;

	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
		} else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra;
	}
}

/**************************************************************************/
void hpsort(unsigned long n, int ra[])
{
	unsigned long i,ir,j,l;
	int rra;

	if (n < 2) return;
 	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
		} else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra;
	}
}
//
//  Now for an indexing function
//
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void indexx(unsigned long n, double arr[], unsigned long indx[])
{
  unsigned long i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  double a;

  istack= new int[NSTACK+1];
  for (j=1;j<=n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j];
	a=arr[indxt];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	SWAP(indx[l+1],indx[ir])
	  }
      if (arr[indx[l]] > arr[indx[ir]]) {
	SWAP(indx[l],indx[ir])
	  }
      if (arr[indx[l+1]] > arr[indx[l]]) {
	SWAP(indx[l+1],indx[l])
	  }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	SWAP(indx[i],indx[j])
	  }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK)
        throw std::range_error("NSTACK too small in indexx.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  delete[] istack;
}
void indexi(unsigned long n, int arr[], int indx[])
{
  int i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  int a;

  istack= new int[NSTACK+1];
  for (unsigned long jj=1;jj<=n;jj++) indx[jj]=jj;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j];
	a=arr[indxt];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	SWAP(indx[l+1],indx[ir])
	  }
      if (arr[indx[l]] > arr[indx[ir]]) {
	SWAP(indx[l],indx[ir])
	  }
      if (arr[indx[l+1]] > arr[indx[l]]) {
	SWAP(indx[l+1],indx[l])
	  }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	SWAP(indx[i],indx[j])
	  }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK)
        throw std::range_error("NSTACK too small in indexx.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  delete[] istack;
}

#undef M
#undef NSTACK
#undef SWAP

// quantiles for an unsorted vector *x of length n
// this function for 1 allocated arrays
// does not return a sorted array
// this method is the same as that in R
void quantiles(double *res,double *x, int nx, double *p,int np)
{
  unsigned long *index = new unsigned long[nx+1];
  indexx(nx,x,index);
  for (int i=1;i<=np;i++) {
    double r=1.0 + (nx-1.)*p[i];
    int ind=int(r);
    assert(ind>=1&&ind<=nx);
    double f=r-double(ind);
    assert(f>=0.0&&f<=1.0);
    if (ind<nx) res[i]=(1.-f)*x[index[ind]]+ f*x[index[ind+1]];
    else res[i]=x[index[ind]];
  }
  delete []index;
}

bool IsSquare(int n, int &siz)
{
  // take the square root
  double b=sqrt(double(n));
  // and round to the nearest integer
  siz=int(b+0.5);
  // does this equal the original value?
  if (siz*siz==n) return true;
  return false;
}



