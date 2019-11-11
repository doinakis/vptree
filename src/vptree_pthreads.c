/*
*Implementation of Vantage point tree using pthreads
*Doinakis Michail
*e-mail: doinakis@ece.auth.gr
*/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include "../inc/vptree.h"

int activeThreads = 0, maxThreads = 4;
pthread_mutex_t mux;
/*In order to pass arguments into the threads we need to create structures. So we create
*one structure for the distances and one structure for the vantage point tree
*/

/* where X is the n by d matrix with the points
*IDS is a matrix containing the indexes of the original set of points
* n is the number of points
* d is the number of dimensions*/
typedef struct pthreads_dist{
  double *X;
  int n;
  int d;
  double *dist; //the distance that its gonna be calculated
  int pointer; // points which calculation is being done
}pthreads_dist;

typedef struct pthreads_vpt{
  double *X;
  int n;
  int d;
  int *IDS;
}pthreads_vpt;
//function that calculates the distance between a point and the vantage point
void* pthreads_calculatedistances(void *thread_argument) // the argument is a pthreads_dist struct
{
  pthreads_dist *corrthread = (pthreads_dist *) thread_argument;//initialize a struct with the argument of the thread
  *(corrthread->dist) = 0; // initialize the distance value to 0
  double temp = 0;
  for(int i=0; i<corrthread->d;i++){//calculate the distance of a single point with vantage point
    temp = temp+ pow((*(corrthread->X + corrthread->pointer * corrthread->d + i)-(*(corrthread->X + (corrthread->n-1)*corrthread->d + i))),2);
  }
  *(corrthread->dist)=sqrt(temp);//save the value the the adress of the struct variable
  return NULL;
}


/*This function takes a thread argument that contains basically a struct of pthreads_vpt type
and returns void pointer which we are gonna cast to a struct value again of pthreads_vpt*/
void * pthreads_Vpt(void * thread_argument){
  vptree *vantageptree = (vptree *) malloc(sizeof(vptree) *1);
  pthreads_vpt *vpt = (pthreads_vpt *) thread_argument;
  double *X = vpt->X; // initialize variables for the tree the points the dimensions and the IDS
  int n = vpt->n;
  int d = vpt->d;
  int *IDS = vpt->IDS;
  if(n == 0){ // if there are no points in the X array return NULL
    return NULL;
  }else if(n == 1){
    vantageptree->vp = (double *) malloc(d *sizeof(double)); // allocate the space for the vantage point
    vantageptree->vpid = IDS[n-1]; // set the value of it to be the last point of the IDS matrix
    vantageptree->md = 0;
    vantageptree->inner = NULL; // since its only one it has no inner nor outer sub-tree
    vantageptree->outer = NULL;
    for(int i=0; i<d; i++){
      vantageptree->vp[i] = *(X + (n-1)* d + i); //setting the vantage point
    }
  }else {
    vantageptree->vp = (double *) malloc(d *sizeof(double));  //allocating the space for the vantage point
    vantageptree->vpid = IDS[n-1]; // we set the vantage point to be the last point of the corresponding tree
    double *distances = (double *) malloc((n -1)*sizeof(double));  // a matrix to keep the distances we are gonna calculate
    double *help_distances = (double *) malloc((n -1)*sizeof(double)); // contains the same valuse as the distances matrix(used just for help)
    pthread_t thread_In,thread_Out; // the two threads for the recursive call
    void * return_In; // the value that its thread will return
    void * return_Out;
    int parallel = 0;
    int innerelements = 0; // counter for the inner sub-tree elements
    int outerelements = 0; // same counter for the outer sub-tree elements
    for(int i=0; i<d; i++){ // setting the last point in the X array to be the vantage point
      vantageptree->vp[i] = *(X + (n -1)*d + i);
    }
    // checks if the activeThreads are lower than the maxThreads and if the number of points is lower than their difference
    if(activeThreads<maxThreads && n<maxThreads-activeThreads){
      pthread_t *threads = (pthread_t *) malloc((n -1)*sizeof(pthread_t));
      pthreads_dist *dist_data = (pthreads_dist *) malloc((n -1)*sizeof(pthreads_dist));
      //parrallel calculation of distances
      for(int i=0; i<n-1; i++){ //each iteration calculates one distance between the vantage point and one element
        dist_data[i].X = X;
        dist_data[i].n = n;
        dist_data[i].d = d;
        dist_data[i].dist = &distances[i];
        dist_data[i].pointer = i;
        pthread_create(&threads[i],NULL,pthreads_calculatedistances,(void *) &dist_data[i]);
      }
      for(int i=0; i<n-1; i++){ // wait for all the threads to end with the calculation
        pthread_join(threads[i],NULL);
      }
      for(int i=0; i<n-1; i++){
        help_distances[i] = distances[i];
      }
      free(dist_data);
      free(threads);
    }else{ // if we cant allocate that many threads the calculation is done sequentially
      double temp=0;
      for(int i = 0; i < n-1; i++){
        for(int j = 0; j < d; j++){
          temp = temp + pow((*(X+i*d+j)-vantageptree->vp[j]),2);
        }
        distances[i] = sqrt(temp);
        help_distances[i] = sqrt(temp);
        temp=0;
      }
    }
    vantageptree->md = median(help_distances,n-1); // calculating the median value
    free(help_distances);
    /* calculation of the inner and outer elements and allocation of the right space for the IDS matrixes
     and the Inner and Outer set of elements */
    for(int i=0; i<n-1; i++){
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
    double *XInner = (double *) malloc(innerelements*d *sizeof(double));
    double *XOuter = (double *) malloc(outerelements*d *sizeof(double));
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
      /*create structures of pthreads_vpt for the recursive calls
      and assigning the values we calculated to the corresponding variables*/
      pthreads_vpt *Inner = (pthreads_vpt *) malloc(sizeof(pthreads_vpt));
      pthreads_vpt *Outer = (pthreads_vpt *) malloc(sizeof(pthreads_vpt));
      Inner->X = XInner;
      Outer->X = XOuter;
      Inner->n = innerelements;
      Outer->n = outerelements;
      Inner->d = d;
      Outer->d = d;
      Inner->IDS = IDSin;
      Outer->IDS = IDSout;
      parallel = 0;
      if(activeThreads < maxThreads){ //checks the number of available threads.
        //if the number is lower
        pthread_mutex_lock (&mux);// lock the activeThreads variable so only the running thread can write on it
        activeThreads += 2; // we increase the number of active threads by 2
        pthread_mutex_unlock (&mux);// unlock of the variable
        pthread_create(&thread_In,NULL,pthreads_Vpt,(void *) Inner);//create 2 threads
        pthread_create(&thread_Out,NULL,pthreads_Vpt,(void *) Outer);
        parallel = 1;//auxilary variable to check if the calculation is done in parallel
      }
      if(parallel){//if it is done in parallel
        pthread_join(thread_In,&return_In); //wait til both of the threads return their values
        pthread_join(thread_Out,&return_Out);
        pthread_mutex_lock (&mux);
        activeThreads -= 2; // decrease the number of activeThreads by 2
        pthread_mutex_unlock (&mux);
        vantageptree->inner = (vptree *) return_In; // assign the returned values to the inner and outer subtree
        vantageptree->outer = (vptree *) return_Out;
      }else{//if there are no more available threads we make the recursion sequentially
        vantageptree->inner = (vptree *) pthreads_Vpt((void *) Inner);
        vantageptree->outer = (vptree *) pthreads_Vpt((void *) Outer);
      }
      free(IDSin);
      free(IDSout);
      free(XInner);
      free(XOuter);
   }
  return (void *) vantageptree;
}


vptree * buildvp(double *X , int n , int d){
  /* setting the IDS for the first call of the Vpt function
  each point has its initial place */
  pthread_mutex_init (&mux, NULL);
  int *IDS=(int *) malloc(n*sizeof(int));
  for(int i=0; i<n; i++)
  {
    IDS[i] = i;
  }
  pthread_t threadid; // creating a thread id to start a thread
  pthreads_vpt *thread_tree;
  void *thread_return;
  thread_tree = (pthreads_vpt *) malloc(sizeof(pthreads_vpt)); // allocating the struct size
  thread_tree->X = X; // and setting the values of the struct
  thread_tree->n = n;
  thread_tree->d = d;
  thread_tree->IDS = IDS;
  pthread_create(&threadid,NULL,pthreads_Vpt,(void *) thread_tree);
  pthread_join(threadid,&thread_return);
  pthread_mutex_destroy(&mux);
  return (vptree *)thread_return; // basically returns the vantage point tree
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
