/*
*Implementation of Vantage point tree cilk
*Doinakis Michail
*e-mail: doinakis@ece.auth.gr
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <sys/time.h>
#include <sys/times.h>

#include "../inc/vptree.h"



// where X is the n by d matrix with the points
//IDS is a matrix containing the indexes of the original set of points
// n is the number of points
// d is the number of dimensions
vptree * Vpt(double *X,int *IDS,int n, int d){
  vptree *vantageptree = malloc(sizeof(vptree)*1); //allocating the right amount of space for the tree
  if(n == 0)
    {
      return NULL; // if there are no points in the X array we return NULL to end the function
    }else if (n == 1) // base case of 1 element
    {
      vantageptree->vp = (double *) malloc(sizeof(double)*d); // allocate the space for the vantage point
      vantageptree->vpid = IDS[n-1]; // set the value of it to be the last point of the IDS matrix
      vantageptree->md = 0;
      vantageptree->inner = NULL;// since its only one it has no inner nor outer sub-tree
      vantageptree->outer = NULL;
      for(int i = 0; i < d; i++){
        vantageptree->vp[i]=*(X + (n-1)*d + i); //setting the vantage point
      }
    }else {
    vantageptree->vp = (double *) malloc(sizeof(double)*d); //allocating the space for the vantage point
    vantageptree->vpid = IDS[n-1]; // we set the vantage point to be the last point of the corresponding tree
    double *distances = (double *) malloc((n-1)*sizeof(double)); // a matrix to keep the distances we are gonna calculate
    double *distance = (double *) malloc((n-1)*sizeof(double)); // contains the same valuse as the distances matrix(used just for help)
    int innerelements = 0;       // counter for the inner sub-tree elements
    int outerelements = 0;      // same counter for the outer sub-tree elements
    for(int i = 0; i < d; i++){     // setting the last point in the X array to be the vantage point
      vantageptree->vp[i]=*(X + (n-1)*d + i);
    }
    //Calculating the distances from all the other points
    double temp=0;
    cilk_for(int i = 0; i < n-1; i++){
      for(int j = 0; j < d; j++){
        temp = temp + pow((*(X+i*d+j)-vantageptree->vp[j]),2);
      }
      distances[i] = sqrt(temp);
      distance[i] = sqrt(temp);
      temp=0;
    }
    vantageptree->md = median(distance,n-1);  // finding the median value
    free(distance);
    /* calculation of the inner and outer elements and allocation of the right space for the IDS matrixes
     and the Inner and Outer set of elements */
    for(int i=0; i<n-1;i++){
      if(distances[i] <= vantageptree->md)
      {
        innerelements++;
      }else
      {
        outerelements++;
      }
    }
    int *IDSin = (int *) malloc(innerelements*sizeof(int));
    int *IDSout= (int *) malloc(outerelements*sizeof(int));
    double *XInner = (double *) malloc(innerelements*d*sizeof(double));
    double *XOuter = (double *) malloc(outerelements*d*sizeof(double));
    /* Assigning points to each matrix weather the the distance of the corresponding point from the vantage point
     is less or equal to the median distance or not */
    int counter1=0;
    int counter2=0;
    for(int i=0; i<n-1;i++){
        if(distances[i] <= vantageptree->md)
        {
          IDSin[counter1]=IDS[i];
          for(int j=0;j<d;j++){
            *(XInner + counter1*d + j) = *(X + i*d +j);
          }
          counter1++;
        }else
        {
          IDSout[counter2]=IDS[i];
          for(int j=0;j<d;j++){
            *(XOuter + counter2*d + j) = *(X + i*d +j);
          }
          counter2++;
        }
      }
      // recursive in order to build all the sub-trees
      vantageptree->inner = cilk_spawn Vpt(XInner,IDSin,innerelements,d);
      vantageptree->outer = Vpt(XOuter,IDSout,outerelements,d);
      cilk_sync;
      /* after we calculate the vantage point tree we free the allocated memory
      from the variables we created to call the function*/
      free(XInner);
      free(XOuter);
      free(IDSin);
      free(IDSout);
    }
  return vantageptree; // returning the vantage point tree
}


vptree * buildvp(double *X , int n , int d){
  /* setting the IDS for the first call of the Vpt function
  each point has its initial place */
  int *IDS=(int *) malloc(n*sizeof(int));
  for(int i=0;i<n;i++)
  {
    IDS[i] = i;
  }
  return (Vpt(X,IDS,n,d)); // basically returns the vantage point tree
}

vptree * getInner(vptree * T){ // returns the inner sub-tree
  return (T->inner);
}

vptree * getOuter(vptree * T){ // returns the outer sub-tree
  return (T->outer);
}
double * getVP(vptree * T){ // returns the vantage point
  return (T->vp);
}
double getMD(vptree * T){ // returns the median distance of the corresponding sub-tree
  return (T->md);
}
int getIDX(vptree * T){ // returns the ID of the vantage point
  return (T->vpid);
}


double quickselect(double *arr, int n, int k) // where arr is the array n is the number of points and k is the k-th smallest element
{
  unsigned long i,ir,j,l,mid;
  double a,temp;

  l=0;
  ir=n-1;
  for(;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir]);
      }
      return arr[k];
    }
    else {
      mid=(l+ir) >> 1;  //binary shift by 1; devides by 2 keeps the lowest int ,the other parts are like quicksort
      SWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      if (j >= k) ir=j-1;
      if (j <= k) l=i;
    }
  }
}
// where arr is the array we want to find the median value
// n is the number of elements arr contains
double median(double *arr, int n){
  if(n%2==1)// if n modulo 2 equals 1 that means that n is odd and we return the quickselect
  {
    return quickselect(arr,n,n/2);
  }else// if n modulo 2 equals 0 we add 2 quickselects for n and n/2-1 and for n and n/2 we add them and divide by 2
  {
    return 0.5*(quickselect(arr,n,n/2-1)+quickselect(arr,n,n/2));
  }

}
