/*kappa keepo
*vptree.h header with the accessors
*Doinakis Michail
*e-mail: doinakis@ece.auth.gr
*/
#ifndef VPTREE_H
#define VPTREE_H
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp; // swaps the values of a and b
extern  int * IDSin;
extern  int * IDSout;
typedef struct vptree_{
  double *vp; // the vantage point
  double md; // the median distance for the corresponding tree
  int vpid; // the index of the input vector of data points
  struct vptree_ *inner;// the inner subtrees
  struct vptree_ *outer;// the outer subtree

}vptree;

//Accessorts and other functions
vptree * buildvp(double *X , int n , int d);
vptree * Vpt(double *X,int *IDS,int n, int d);
vptree * getInner(vptree * T);
vptree * getOuter(vptree * T);
double * getVP(vptree * T);
double getMD(vptree * T);
int getIDX(vptree * T);
double quickselect(double *arr, int n, int k);
double median(double *arr, int n);






#endif
