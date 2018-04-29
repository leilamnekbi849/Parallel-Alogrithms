#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#define max(x,y) ((x>y)? x:y)
#define min(x,y) ((x<y)? x:y)

#define np 16  // no of processors 
#define nfiles 128  // no of files 
#define filesize 128  // no of inputs in each file 
#define cK 32 // this is number of clusters K cK and it has to be a multiple of np 
#define  max_transfers nfiles*filesize / (cK) // assuming there aren't more than 30% transfers in between clusters is a wrong assumption
#define range 128 // this is the max int on any file


int cluster[cK/np][nfiles*filesize], cluster_pointer[cK/np];

//these are buffers used to send and recieve points 
// NOTE : the 10,000 constant is a current assumption that there is not more movement than that. This shouldnt be hard coded and should be set to avg no of points in a cluster/5 or something
// so a better thing might be ~ (nfiles*filesize / cK) / 4 ~= 128 currently 
// the ptr buffers are just to keep a track of how many points need to be sent and recieved 


int send_points[np + 1][cK/np + 1][max_transfers], recv_points[np + 1][cK/np + 1][max_transfers],send_points_ptr[np+1][cK/np + 1], recv_points_ptr[np+1][cK/np+1];
  

// cluster pointer keeps track of how many elements are present in a given  cluster 
// each node has k/np clusters and max size of a cluster is total no of inputs  
// possibility here is to use dynamic arrays instead for cluster for better memory instead of just alloting entire size of input 
// another suggestion is that since data is evenly distributed, you should assume that no cluster gets huge amount of data, i.e you can assume roughly 
// each cluster gets nfiles*filesize / __K data so you can multiply this by 4 to be on the safer side 

float centers[np + 1][cK/np + 1]; // centers[i][j] is the jth cluster center in the ith node 


int blocksize = nfiles/np; // the number of files each processor is responsible for 

int data[10001],dataCount=0;


char * toArray(int number) {
  // converts an integer to a character array
  //printf("file num is %d\n", number);
  if(!number) return "0";

  int n = log10(number) + 1,i = 0;
  char *numberArray = calloc(n, sizeof(char));
  
  for (i = 0; i < n; ++i, number /= 10 )
    numberArray[n-1-i] = (number % 10) + '0';
  
  return numberArray;
}



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
  //  printf("%s||\n",fname);
  
  FILE* F =  fopen(fname,"r");
  if(!F)
    printf("no file found for file name fname %s\n",fname);
 
  
  while(fscanf(F, "%d", &data[dataCount]) != EOF && dataCount<filesize*blocksize) // max number of integers read are blocksize * filesize
    dataCount++;
      
  //  printf("Read succesfully %d integers\n", dataCount);
 


}

void populate_data (int node){
  // reads all the files the current node is responsible for 
  
  char  fname[1000];

  int i = 0; // the number of files each processor will have to read 
  // printf("%d is the blocksize\n",blocksize);

  for(i = blocksize*node; i < blocksize*(node+1); i++){

    get_file_name(fname,i); // the name of the ith file is stored in fname 
    read_integers(fname);  // integers from fname are read into global data 
  
  }


}


void populate_clusters(int clustersPerNode,int node){
  // this function initalizes the clusters 
  int dataPerCluster = dataCount/clustersPerNode;  
  // printf("Clusters per node = %d data count = %d dataPerCluster = %d\n ",clustersPerNode,dataCount, dataPerCluster);
  
  int i = 0,j = 0;
  for(i = 0;i < clustersPerNode;i++) {
    for(j = i*dataPerCluster; j<(i+1)*dataPerCluster; j++){
      cluster[i][cluster_pointer[i]] = data[j]; // ith cluster adds new data 
      cluster_pointer[i] ++; // increment its pointer
      
      
    }
    // !!!!!!!!!!!!!!!!!!!!!!! CHECKPOINT 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! each cluster recieves appropirate amount of data 
    if(cluster_pointer[i] * cK != nfiles * filesize)
      printf("CHECKPOINT0 FAILED: %d th cluster has %d data at the beginning which is unexpected because total input (%d) != data per cluster(%d) * cK(%d)\n",i,cluster_pointer[i],nfiles*filesize,cK);
  }
  
}



void initialize_all_means(int clustersPerNode){

  int i=0,j=0;
  float step = range / cK, current = 1;
  for(i=0; i<np; i++)
    for(j=0; j<clustersPerNode; j++){
       centers[i][j] = current;
       current += step;
       //    printf("%.3f\n",centers[i][j]);
    }

}

void calculate_node_means(int clustersPerNode, int node){
  // calculates means for the given nodes clusters 
  
  int i=0,j=0;

  for(i=0; i<clustersPerNode; i++){
    int sum = 0;

    for(j=0; j<cluster_pointer[i]; j++)
      sum += cluster[i][j];
    
    float clusterMean = sum*1.0/cluster_pointer[i]; // using integers for approximation 
    
    printf("%.3f is the clusterMean for  %d#%d\n", clusterMean,node,i);
    
    centers[node][i] = clusterMean;
  }


}


void bcast_and_get_means(int clustersPerNode, int node){
  // function broadcasts means of given node's clusters and gets means from other nodes' clusters too

  int i=0,j=0,k=0;
  // float packets[cK];
  
  calculate_node_means(clustersPerNode, node);
  
  for(i=0; i<clustersPerNode; i++){
    
    float clusterMean = centers[node][i];
    
    // now this needs to be broadcast 
    // multiple options here : calculate all means and then broadcast or broadcast one by one i.e batch vs individual
    // note : this has to asychronous
    // broadcase needs to be of the form <mean, node#,  cluster#> which here is <clusterMean, node, i> --- could be an integer a[3]
    
    float packet[2] = {clusterMean,i};
     MPI_Status status;
    for(k=0;k<np;k++){
      if(k!=node){
	MPI_Sendrecv(&clusterMean,1,MPI_FLOAT,k,1,
		    &centers[k][i],1,MPI_FLOAT,k,1,MPI_COMM_WORLD,&status);
	
	/*MPI_Sendrecv(clusterMean,1,MPI_FLOAT,k,1,
		     centers[k][i],1,MPI_FLOAT,k,1,MPI_COMM_WORLD,&status);
	MPI_Sendrecv(i,1,MPI_INT,k,2,
		     ,3,MPI_FLOAT,k,1,MPI_COMM_WORLD,&status);
	MPI_Sendrecv(packet,3,MPI_FLOAT,k,1,
		     packets[k],3,MPI_FLOAT,k,1,MPI_COMM_WORLD,&status);
	*/
      }
    }
    
  }
    
  
  if(node==0)
    for(i=0;i<np;i++){
      printf("node %d -- centers -- \n",i);
      for(j=0;j<cK/np;j++){
	printf("%.3f ",centers[i][j]); // sahi hai kya babua  
      }
      printf("\n");

    }
}


void shift_buffers(int clustersPerNode,int node){
  // shift all the points in cluster [][]  marked as -1 as they have been removed 
  int i=0,j=0,k=0;
  
   for(i=0;i<clustersPerNode;i++){
    // we are talking about i th cluster in the given node
     int nulls = 0,filler = 0;
     for(j=0;j<cluster_pointer[i];j++)
       if(cluster[i][j]==-1)
	 nulls++;
       else 
	 cluster[i][filler++] = cluster[i][j];
     
     cluster_pointer[i]-=nulls;
     
    
     if(node==1)
       printf("%d %d - elements filled = NEW clustersize\n", filler,cluster_pointer[i]);

   }
  
}

void init_send_recv_buffers(int clustersPerNode){
  int i=0,j=0;
  for(i=0;i<np;i++)
    for(j=0;j<clustersPerNode;j++)
      send_points_ptr[i][j] = recv_points_ptr[i][j] = 0;
  

}


//int new_points[4*nfiles*filesize/cK];
void send_recv_points(int clustersPerNode, int node, int outgoing, int initial_points){

  // outgoing is the number of out going points 
  // now that you know which points are to be sent and which will stay, do the communication brother
  int in = 0, ic = 0, incoming = 0, dummy = 0;

  MPI_Status status,status1;
  for(in = 0; in<np; in++)
    if(in!=node)
      for(ic = 0; ic<clustersPerNode; ic++){
	int zero = 0,cnt = 0,tag1 = (node+in)*100 + (node-in)*(node-in) + 10*ic; // tags need to be matched -> exchange me involved elements are node*1000,ic*100,in*10,ic*1 
	
	MPI_Sendrecv(&send_points_ptr[in][ic],1,MPI_INT,in,tag1,
		     &cnt,1,MPI_INT,in,tag1,MPI_COMM_WORLD,&status);
	
	int tag2 = (node+in)*1000 + (node-in)*(node-in) + ic;
	if(cnt)
	  MPI_Sendrecv(send_points[in][ic],max(1,send_points_ptr[in][ic]),MPI_INT,in,tag2,
		       &cluster[ic][cluster_pointer[ic]],cnt,MPI_INT,in,tag2,MPI_COMM_WORLD,&status1);                 //     &new_nodes[incoming] -> this is an option to collect data
	else // recieve a dummy variable that is garbage value stored in send buffer
	  MPI_Sendrecv(send_points[in][ic],max(1,send_points_ptr[in][ic]),MPI_INT,in,tag2,
		       &dummy,1,MPI_INT,in,tag2,MPI_COMM_WORLD,&status1);
	
	cluster_pointer[ic] += cnt;
	incoming += cnt;
	
      }

  if(!node)
    printf(" incoming(%d) - outgoing(%d) = %d\n",incoming,outgoing,incoming-outgoing);  // need to put a check point here that sum of all the nodes for this value is 0
  

  int final_points = 0; // count
  int j = 0;
    for(ic=0;ic<clustersPerNode;ic++)
      final_points += cluster_pointer[ic];

    if(incoming-outgoing != final_points-initial_points) // makes sure that node recieves as many points as expected into every cluster
      printf("CHECKPOINT3 FAILED for node %d \n",node);

}


void re_clusterify(int clustersPerNode, int node){
  // initialize counts to 0 first  because we dont need to send anything before clusterifying so start with new buffers 
  init_send_recv_buffers(clustersPerNode);

  // calcluate clusters again
  // idea I am thinking is create an array like A[node][cluster][1000] and Aptr[node][cluster] . keep inserting into closest cluster in A  and updating Aptr 
  // later on send 2 things 
  // for every pair [node,cluster] send the size of points to be sent as stored in Aptr[node][cluster] 
  // recieve this in the respective node
  // if the number is >0 then recieve an array of points for every node,cluster pair
  // accept the points at recieving cluster and then add them to respective cluster array/data structure and update its count
  
  // take care of internal and external transfers : if there is an internal transfer then no need to used send
  
  int i,j,k,m;
  int initial = 0, external = 0, internal = 0, intact = 0, remaining = 0; 
// note that internal and intact are overlapping because some internal movements happen while the loop is being executed --> cluster_pointer is updated dynamically
// so they add to existing intact values <-> so its wrong to assume initial = external - (internal + intact)

  // for testing purposes, lets calculate the initial number of values at a given node
 

  // TEST LOOP
  for (i=0;i<clustersPerNode;i++) // ith cluster
    initial += cluster_pointer[i]; // initial number of values at cluster i before we re clusterify 
  

  // MAIN LOOP BEGINS ----------------------------------------------------------------------------------------------------------------------------------
  for (i=0;i<clustersPerNode;i++){
    // ith cluster
    for(j=0; j<cluster_pointer[i];j++){
	// j point in the ith cluster of current node

      float closest_distance = 100000000*1.0;
      int closest_node = -1, cluster_id = -1;

      // loop for finding nearest cluster center
	for (k=0;k<np;k++){
	  // the kth node 
	  for(m=0;m<clustersPerNode;m++){
	    // the m the cluster in the kth node
	    float distance_to_cluster = abs(centers[k][m]-cluster[i][j]);
	    if(distance_to_cluster < closest_distance){
	      closest_node = k;
	      cluster_id = m;
	      closest_distance = distance_to_cluster;
	    }
 
	    // end of for loop of m
	  }
	  // end of for loop for k 
	}
	
	// -- we found the closest cluster center so now time to update 
	// -- send the point to closest_node,cluster_id
	// -- remove the point from here -> for now just mark it as -1 and then do a shift operation at the end so this doesnt make the loop n^5 
	
	if(closest_node != node){
	  external ++;
	  int tmp_ptr = send_points_ptr[closest_node][cluster_id];
	  send_points[closest_node][cluster_id][tmp_ptr] = cluster[i][j];
	  send_points_ptr[closest_node][cluster_id] ++; // add point to buffer <closest node, closest cluster id> 
	  cluster[i][j] = -1; // remove point from this cluster
	}
	else if (cluster_id != i){ 
	  // move the point from i th cluster to the cluster_id th cluster
	  internal ++;
	  int tmp_ptr = cluster_pointer[cluster_id];
	  cluster[cluster_id][tmp_ptr] = cluster[i][j]; // add this point to cluster id cluster 
	  cluster_pointer[cluster_id]++; 
	  cluster[i][j] = -1; // remove this point from current cluster 
	}
	else 
	  intact ++;
	//if(node==1)
	//  printf("point %d closest to  %f in node %d\n",cluster[i][j], centers[closest_node][cluster_id], closest_node);
	// the end of all points in given cluster 
    }
    
    // end i loop i.e end of movement of points from all clusters in given node
  }
  shift_buffers(clustersPerNode,node); // remove all the points set to -1 from the nodes clusters 
 
  // END OF MAIN LOOP--------------------------------------------------------------------------------------------------------------------------------------------


 // TEST LOOP 
  for (i=0;i<clustersPerNode;i++) // ith cluster
    remaining += cluster_pointer[i]; // remaining number of values at cluster i before we re clusterify 
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHECKPOINT 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for node, Initial data = Data exiting + Data staying 
  if(initial != external + remaining)
    printf("CHECKPOINT1 FAILED:\nInitial values(%d) == Exiting values(%d) + Remaining values(%d) ? \n",initial,external,remaining); // check whether we have same number of elements that we started
    // printf("%d %d %d ext,int,static for node #%d\n", external,internal,intact,node);
  
  int test_outgoing = 0, flag = 1;
  for(i=0;i<np;i++){
    if(i!=node){
      for(j=0;j<clustersPerNode;j++){
	test_outgoing += send_points_ptr[i][j];
	for(k=0;k<send_points_ptr[i][j];k++)
	  if(send_points[i][j][k]<0 || send_points[i][j][k]>range)
	    flag = 0;
      }
    }    
   
  }
  //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHECKPOINT 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! data outgoing is fed properly to the outgoing buffer 
  if((test_outgoing != external) || (flag==0)) 
    printf("CHECKPOINT2 FAILED:\n2.1 -> outgoingbuffersize(%d) == expectedsize(%d)?  AND \n2.2 -> all values in expected range (1) == %d \n",test_outgoing,external,flag);
  
  send_recv_points(clustersPerNode,node,external,initial);


}



void check_stop_condition(){



}

int main (int argv, char ** argc) {

  MPI_Init(&argv,&argc);
  int node,csize,i,temp;
  MPI_Comm_rank(MPI_COMM_WORLD,&node);
  MPI_Comm_size(MPI_COMM_WORLD,&csize);
  
  populate_data(node);
  populate_clusters(cK/np,node); // clusters per node is the parameter
  // recieve_means(cK/np,node); 
  initialize_all_means(cK/np);
  re_clusterify(cK/np, node);


  //  while(some condition) 
  // bcast_and_get_means(cK/np,node); // same as above
  // re_clusterify(cK/np,node);

  int sum = 0;
  for(i=0;i<dataCount;i++)
    sum += data[i];

  
  //  printf("sum for node#%d = %d for %d integers\n",node,sum,dataCount);

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
