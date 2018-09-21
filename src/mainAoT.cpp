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
#include<fstream>
#include<iostream>
#include<string>
#include<climits>
#include<cstring>
#include<algorithm>

#include"visAoT.h"
#include"buildSynthetic.h"

//#define DEBUG 1

//#define SIZE 524288

#define ITERS 10

int main(int argc, char* argv[]) {

	if(argc != 3) {
		printf("Usage: ./visAoT <size (2^x)> <# threads>\n");
		exit(1);
	}

int size=pow(2,atoi(argv[1]));
int runTimes[ITERS];

int* peaks = (int*)malloc(size*sizeof(int));
int numThreads = atoi(argv[2]);

  int* xvals = (int*)malloc(size*sizeof(int));
  int* visIdx = (int*)malloc(size*sizeof(int));

  srand(time(NULL));

  std::chrono::steady_clock::time_point beginSeq;
  std::chrono::steady_clock::time_point endSeq; 

  for(int iter=1; iter<ITERS; iter++) {
    for(int i=0; i<size; i++) {
      xvals[i] = i; 
      visIdx[i]=0;
    }
    buildRandom(size, peaks);
//  buildFlat(size, peaks);
//  buildParabolic(size, peaks);

  beginSeq= std::chrono::steady_clock::now();

  vis1DAoT(xvals, peaks, visIdx, size, numThreads);

  endSeq= std::chrono::steady_clock::now();
  runTimes[iter] = std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();

}

int mean;
int total=0;
float stdev=0;

std::sort(runTimes,runTimes+ITERS);
mean = runTimes[ITERS/2];
for(int i=0; i<ITERS-1; i++) {
  stdev += (pow((float)(runTimes[i]-mean),2.0))/(float)ITERS;
  total += runTimes[i];
}
stdev = sqrt(stdev);
  
//printf("%d %.4f\n", mean, stdev);
printf("%lf\n", ((double)total/ITERS)/1000.0);


#ifdef DEBUG
int* testVals = (int*)malloc(size*sizeof(int));
for(int i=0; i<size; i++) testVals[i]=0;
float tempSlope;
for(int i=0; i<size; i++) {
  tempSlope = -99999.99;
  for(int j=i; j<size; j++) {
    if(calcSlope(xvals[i], peaks[i], xvals[j], peaks[j]) > tempSlope) {
      testVals[i]++;
      testVals[j]++;
      tempSlope = calcSlope(xvals[i],peaks[i],xvals[j],peaks[j]);
    }
  }
}
bool correct=true;
for(int i=0; i<size; i++) {
  if(visIdx[i] != testVals[i]){
    printf("vis[%d]:%d, test:%d\n", i, visIdx[i], testVals[i]);
    correct=false;
  }
}
if(!correct)printf("incorrect result!\n");
#endif
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
