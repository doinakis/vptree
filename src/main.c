/*καππα
*Implementation of Vantage point tree
*Main Implementation for test Usage
*Doinakis Michail
*e-mail: doinakis@ece.auth.gr
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <sys/times.h>
#include "../inc/vptree.h"
#define N 1000000
#define D 100

struct timeval startwtime,endwtime;
double p_time;



int main()
{
  srand(time(NULL));
  double *X=(double *)malloc(N*D*sizeof(double));
  int i=0;
  int j=0;
  for(i=0;i<N;i++){
    for(j=0;j<D;j++){
      *(X+i*D+j)=(rand() % 10001) / 10000.0; //creating random double points
    }
  }
  vptree * T;
  T=(vptree *) malloc( sizeof(vptree)*1);
  gettimeofday(&startwtime,NULL);
  T=buildvp(X,N,D);
  gettimeofday(&endwtime,NULL);
  p_time = (double)((endwtime.tv_usec-startwtime.tv_usec)/1.0e6+endwtime.tv_sec-startwtime.tv_sec);
  printf("Time %f\n",p_time);

   free(T);
   free(X);
  return 0;
}
