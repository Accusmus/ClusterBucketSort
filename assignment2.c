//first.cpp Adding numbers using two nodes C++ version
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <climits>
#include <iostream>
#include <fstream>
#include <sstream>

#include "mpi.h"

typedef unsigned long ULONG;

using namespace std;

int main(int argc,char* argv[])
{
  MPI_Status Stat;

  MPI::Init(argc, argv);

  int nitems;
  const float xmin = 0.0; 
  const float xmax = 100.;
  int nbuckets = 1000;

  int myid = MPI::COMM_WORLD.Get_rank();
  int numproc = MPI::COMM_WORLD.Get_size(); 

  if(argc == 2){ 
    nitems = atoi(argv[1]);
  }else{
    cout << "Please Provide the number of items as arguement" << endl;
    exit(1);
  }

  int itemsperproc = nitems/numproc;

  float *data=(float*)malloc(nitems*sizeof(float));
  float *sub_data=(float*)malloc((nitems/numproc)*sizeof(float));

  // master process
  if(myid == 0){
    //cout << "Number of random numbers: " << nitems << endl;
    //cout << "{";
    for(int i = 0; i < nitems; i++){
      data[i] = drand48()*(xmax-xmin-1)+xmin;
      //cout << data[i] << " ";
    }
    //cout << "}" << endl;
  }
  
  //send a small subset of the random numbers to each process
  MPI_Scatter(data, itemsperproc, MPI_FLOAT, sub_data, itemsperproc, MPI_FLOAT, 0, MPI_COMM_WORLD); 
  //cout << myid << ": {";
  for(int i = 0; i < itemsperproc; i++){
    //cout << sub_data[i] << " ";
  }
  //cout << "}" << endl;

  //No need to sort the data into a number of small buckets
  //the number of buckets is the number of processors
  int bucketsize[numproc];
  int offset[numproc];

  for(int i = 0; i < numproc; i++){
    bucketsize[i] = 0;
    offset[i] = 0;
  }

  float * smallbuckets = (float *)malloc(itemsperproc * sizeof(float));
  float * bigbucket = (float *)malloc(itemsperproc * numproc * sizeof(float));

  int count = 0;
  
  for(int i = 0; i < numproc; i++){ 
    float min = (((xmax - xmin)/numproc) + xmin) * i;
    float max = (((xmax - xmin)/numproc) + xmin) * (i + 1);

    for(int j = 0; j < itemsperproc; j++){
      if(sub_data[j] > min && sub_data[j] < max){
	smallbuckets[count] = sub_data[j];
	sub_data[j] = -1.0;
	bucketsize[i] ++;
	count++;
      }
    }
  }

  int recvcount[numproc];
  int recvdisp[numproc];

  for(int i = 0; i < numproc; i++){
    recvcount[i] = itemsperproc;
    recvdisp[i] = i * itemsperproc;
  }

  for(int i = 0; i < numproc; i++){
    if(i == 0) continue;
    offset[i] = offset[i-1] + bucketsize[i-1];
  }

  MPI_Alltoallv(smallbuckets, bucketsize, offset, MPI_FLOAT, bigbucket, recvcount, recvdisp, MPI_FLOAT, MPI_COMM_WORLD);

  cout << "{";
  for(int i = 0; i < itemsperproc*numproc; i++){
    cout << bigbucket[i] << " ";
  }
  cout << "}" << endl;
  

  MPI::Finalize();
}
