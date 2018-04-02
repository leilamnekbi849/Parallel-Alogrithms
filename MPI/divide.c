#include <stdio.h>
#include <mpi.h>
#define np 8


char * toArray(int number) {

  int n = log10(number) + 1,i;
  char *numberArray = calloc(n, sizeof(char));
  
  for ( i = 0; i < n; ++i, number /= 10 )
    numberArray[i] = number % 10;
  
  return numberArray;
}

int main (int argv, char ** argc) {
  MPI_Init(&argv,&argc);
  int node,csize,i;
  MPI_Comm_rank(MPI_COMM_WORLD,&node);
  MPI_Comm_size(MPI_COMM_WORLD,&csize);
  char * sbegin = "/MPIinputs/input", *send =".txt";  
  char * mid = toArray(node);
  
  printf("%s %s %s\n", sbegin,mid,send);
  MPI_Finalize();
  return 0;
}
