#include <stdio.h>
#include <mpi.h>
#define np 8

int main (int argv, char ** argc) {
  MPI_Init(&argv,&argc);
  int node,csize;
  MPI_Comm_rank(MPI_COMM_WORLD,&node);
  MPI_Comm_size(MPI_COMM_WORLD,&csize);
  int val = 1;
  
  int to = (node+1) % np, from = (node-1+np)%np;       
     
  if(node==0){
    val = node;
    MPI_Send(&val,1,MPI_INT,to,0,MPI_COMM_WORLD);
    MPI_Recv(&val,1,MPI_INT,from,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    printf("%d is node and val of parent is %d\n",node,val);
    // printf("%d is node and val of to,from is %d %d\n",node,to,from);
    // break;
  } else {
    MPI_Recv(&val,1,MPI_INT,from,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    printf("%d is node and val of parent is %d\n",node,val);
    val = node;
    MPI_Send(&val,1,MPI_INT,to,0,MPI_COMM_WORLD);
    // printf("%d is node and val of to,from is %d %d\n",node,to,from);

  }
 
  printf("out of here %d\n",node);
   
  MPI_Finalize();
  return 0;
}

