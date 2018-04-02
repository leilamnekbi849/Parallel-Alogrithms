#include <stdio.h>
#include <mpi.h>
#define np 16  // no of processors 
#define nfiles 128  // no of files 
#define filesize 128  // no of inputs in each file 

void  get_file_name (char * sbegin, int number){
  // return the file name thats to be read by the node passed as an argument
  char *send =".txt";  
  char * mid = toArray(number);
  sbegin[0]='\0';
  strcat(sbegin, "MPIinputs/input");
  strcat(sbegin,mid);
  strcat(sbegin,send);
}


void read_integers(char * fname){
  // reads integers from file to data array
  printf("%s||\n",fname);
  
  FILE* F =  fopen(fname,"r");
  if(!F)
    printf("no file found for file name fname %s\n",fname);
 
  
  while(fscanf(F, "%d", &data[dataCount]) != EOF && dataCount<filesize)
    dataCount++;
      
  printf("Read succesfully %d integers\n", dataCount);
 


}

void populate_cluster (int node){
  // reads all the files the current node is responsible for 

  char  fname[1000];

  int block_size = nfiles/np; // the number of files each processor will have to read 

  for(i = block_size*node; i < (block_size+)*(node+1); i++){

    get_file_name(fname,i);
    read_integers(fname);
  
  }


}

int main (int argv, char ** argc) {

  MPI_Init(&argv,&argc);
  int node,csize,i,temp;
  MPI_Comm_rank(MPI_COMM_WORLD,&node);
  MPI_Comm_size(MPI_COMM_WORLD,&csize);
  
  populate_cluster(node);

  
  
  printf("sum for node#%d = %d\n",node,sum);

  int othersSum = 0;
  if(node)
    MPI_Send(&sum,1,MPI_INT,0,0,MPI_COMM_WORLD);
  else 
    for(i=1;i<np;i++){
      MPI_Recv(&othersSum,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      sum = sum + othersSum;
    }
 
  if(!node)
    printf("NET_SUM = %d\n",sum);
  
  MPI_Finalize();
  return 0;
}
