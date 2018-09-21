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
#include "vis-succinct.h"
#include<fstream>
#include<iostream>
#include<string>
#include<climits>
#include<algorithm>
#include"buildSynthetic.h"

#define ITERS 10


int main(int argc, char* argv[]) {

	if(argc != 3) {
		printf("Usage: ./linPar <size (2^x)> <# threads>\n");
		exit(1);
	}

  int size = pow(2,atoi(argv[1]));
int* peaks = (int*)malloc(size*sizeof(int));

std::string fileName;
std::ifstream file;


  int runTimes[ITERS];

  int numThreads = atoi(argv[2]);

  int* xvals = (int*)malloc(size*sizeof(int));
  int* visIdx = (int*)malloc(size*sizeof(int));



  std::chrono::steady_clock::time_point beginSeq;
  std::chrono::steady_clock::time_point endSeq; 

for(int iter=0; iter<ITERS; iter++) {
  for(int i=0; i<size; i++) {
    xvals[i] = i; 
    visIdx[i]=0;
//printf("%d ", peaks[i]);
  }
  buildRandom(size, peaks);
//  buildFlat(size, peaks);
//  buildParabolic(size, peaks);

  beginSeq= std::chrono::steady_clock::now();

  vis1D(xvals, peaks, visIdx, size, numThreads);

  endSeq= std::chrono::steady_clock::now();
  runTimes[iter] = std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();

}

int mean;
float stdev=0;
int total=0;

std::sort(runTimes,runTimes+ITERS);
mean = runTimes[ITERS/2];
for(int i=0; i<ITERS-1; i++) {
  total += runTimes[i];
  stdev += (pow((float)(runTimes[i]-mean),2.0))/(float)ITERS;
}
stdev = sqrt(stdev);

printf("%lf\n", ((double)total)/ITERS/1000.0);

//printf("%lf\n", (double)mean/1000.0);
//printf("%.4f\n", stdev);


  free(xvals);
  free(peaks);
  free(visIdx);
  return 0;
}

