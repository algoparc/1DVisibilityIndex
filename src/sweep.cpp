/*
Copyright 2016-2018 Ben Karsin, Nodari Sitchinava

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/


#include <chrono>
#include<omp.h>
#include "SimpleRBI.h"
#include<fstream>
#include<iostream>
#include<string>
#include<climits>
#include"buildSynthetic.h"

#define ITERS 10

int main(int argc, char* argv[]) {

	if(argc != 2) {
		printf("Usage: ./sweep <size (2^x)>\n");
		exit(1);
	}

  int size = pow(2,atoi(argv[1]));

  int* xvals = (int*)malloc(size*sizeof(int));
  int* peaks = (int*)malloc(size*sizeof(int));
  int* visIdx = (int*)malloc(size*sizeof(int));

  int total=0;

  for(int iter=0; iter < ITERS; iter++) {
  srand(time(NULL));
  for(int i=0; i<size; i++) {
    xvals[i] = i; 
    visIdx[i]=0;
  }
  buildRandom(size, peaks);
//  buildFlat(size, peaks);
//  buildParabolic(size, peaks);


  std::chrono::steady_clock::time_point beginSeq;
  std::chrono::steady_clock::time_point endSeq; 

  beginSeq= std::chrono::steady_clock::now();

//  vis1D(xvals, peaks, visIdx, size, numThreads);
  sweepVis(xvals, peaks, visIdx, size);

  endSeq= std::chrono::steady_clock::now();

  total += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();
  }

  printf("%lf\n", ((double)total)/ITERS/1000.0);

/*
for(int i=0; i<size; i++) {
  printf("%d ", visIdx[i]);
}
printf("\n");
*/

  free(xvals);
  free(peaks);
  free(visIdx);
  return 0;
}
