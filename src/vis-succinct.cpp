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


//#include "AOT.h"
//#include "chazel/params.h"
#include "chazel/chazelBits.cpp"
//#include "convexhull.h"
#include "BuildConvexHull.h"
#include <chrono>
#include<stdlib.h>
//#include <jemalloc/jemalloc.h>
#include<omp.h>
#include<algorithm>
#include<parallel/algorithm>
#include<stdio.h>
#include<string.h>

//#define DEBUG 1

template<typename T, typename V>
struct timeAndVals {
  T time;
  V val;
};

struct timeAndLine {
  float time;
  Line line;
};

bool cmpStructUp(timeAndVals<float,float> a, timeAndVals<float,float> b) { return a.val < b.val; }
bool cmpStructDown(timeAndVals<float,float> a, timeAndVals<float,float> b) { return a.val > b.val; }
bool cmpLineStructDown(timeAndLine a, timeAndLine b) { return a.line.slope > b.line.slope; }

template<typename T>
void printTree(T* leaves, T* root, WORDTYPE** bitVec, int N) {
  int levels = (int)ceil(log2((float)(N)))+1; 
  WORDTYPE tempWord;
  int maskVal = 1 << (WORDSIZE-1);
  for(int i=0; i<N; i++) {
    printf("%f ", root[i]);
  }
  printf("\n");
  for(int i=0; i<levels; i++) {
    for(int j=0; j<N/WORDSIZE; j++) {
      tempWord = bitVec[i][j];
printf("- %d - ", __builtin_popcount(tempWord));
      for(int bit=0; bit<WORDSIZE; bit++) {
        printf("%d ", (tempWord & maskVal)/maskVal);
        tempWord = tempWord << 1;
      }
printf("  ");
    }
    printf("\n");
  }
  for(int i=0; i<N; i++) {
    printf("%f ", leaves[i]);
  }
  printf("\n");
}


struct cmpSlopes {
  bool operator()(Line a, Line b) {
  return(a.slope <= b.slope);
}};

struct cmpSlopesDown {
  bool operator()(Line a, Line b) {
  return(a.slope >= b.slope);
}};

struct cmpPi {
  bool operator()(float a, float b) { return a <= b; }
};

struct cmpPiDown {
  bool operator()(float a, float b) { return a >= b; }
};

struct cmpTime {
  bool operator()(float a, float b) { return a < b; }
};
struct cmpTimeDown {
  bool operator()(float a, float b) { return a > b; }
};

struct aboveness {
  bool operator()(float x, float y, Line line) { 
    return (y > ((line.icp)+(x*line.slope))); }
};

float calcSlope(int x1, int y1, int x2, int y2) {
  int deltaX = x2-x1;
  int deltaY = y2-y1;
 
  return (float)deltaY/(float)deltaX; 
}
float calcIcp(int x, int y, float slope) {
  return (float)y-((float)x*slope) ;
}

// Initialize critical rays and visIdx for groups of WORDSIZE
void initCritRays(Point* points, int* visIdx, Ray* Rcrit, Ray* Lcrit) {
  float tempTangent;
  int j;
  // Initilize critical rays
  for(int i=0; i<WORDSIZE; i++) {
    Rcrit[i].x = points[i].x;
    Lcrit[i].x = points[i].x;
    Rcrit[i].y = points[i].y;
    Lcrit[i].y = points[i].y;
    visIdx[i]=1; // All rays initially see themselves

    Rcrit[i].slope = -99999999.0;
    Lcrit[i].slope = 99999999.9;
/*
    Lcrit[i].slope = Lcrit[i].y+1;
    Lcrit[i].icp = calcIcp(Lcrit[i].x, Lcrit[i].y, Lcrit[i].slope);
    Rcrit[i].slope = -(Rcrit[i].y+1);
    Rcrit[i].icp = calcIcp(Rcrit[i].x, Rcrit[i].y, Rcrit[i].slope);
*/
//    Rcrit[i].icp = 9999.9;
  }
  Lcrit[0].slope = ((float)Lcrit[0].y)+1.0;
  Lcrit[0].icp = calcIcp(Lcrit[0].x, Lcrit[0].y, Lcrit[0].slope);
  Rcrit[WORDSIZE-1].slope = -(((float)Rcrit[WORDSIZE-1].y)+1.0);
  Rcrit[WORDSIZE-1].icp = calcIcp(Rcrit[WORDSIZE-1].x, Rcrit[WORDSIZE-1].y, Rcrit[WORDSIZE-1].slope);
  for(int i=0; i<WORDSIZE-1; i++) {
    for(int j=i+1; j<WORDSIZE; j++) {
      tempTangent = calcSlope(points[i].x, points[i].y, points[j].x, points[j].y);
      if(tempTangent > Rcrit[i].slope) {
        visIdx[i]++;
        Rcrit[i].slope = tempTangent;
        Rcrit[i].icp = calcIcp(points[i].x, points[i].y, tempTangent);       
      }
    }
  }
  for(int i=WORDSIZE-1; i>0; i--) {
    for(int j=i-1; j>=0; j--) {
      tempTangent = calcSlope(points[i].x, points[i].y, points[j].x, points[j].y);
      if(tempTangent < Lcrit[i].slope) {
        visIdx[i]++;
        Lcrit[i].slope = tempTangent;
        Lcrit[i].icp = calcIcp(points[i].x, points[i].y, tempTangent);       
      }
    }
  }
}

/*
void BFSIntersections(float* initTimes, Line* dualLines, Ray* Rcrit, Ray* Lcrit, int setSize) {
  int rCount[setSize*2];
  int aboveCount[setSize*2];
  int belowCount[setSize*2];

  for(int i=0; i<setSize; i++) {
    rCount[i] = 0;
    aboveCount[i] = 0;
    belowCount[i] = 0;
    for(int j=0; j<setSize; j++) {
      if(Rcrit[i].dualX() < Lcrit[j].dualX()) {
        rCount[i]++;
        if(aboveness(Rcrit[i].dualX(), Rcrit[i].dualY(), dualLines[setSize+j])) {
          belowCount[i]++;
        }
        if(aboveness(Lcrit[j].dualX(), Lcrit[j].dualY(), dualLines[i])) {
          aboveCount[i]++;
        } 
      }
    }
  }
  for(int i=0; i<setSize; i++) {
    rCount[setSize+i] = 0;
    aboveCount[setSize+i] = 0;
    belowCount[setSize+i] = 0;
    for(int j=0; j<setSize; j++) {
      if(Lcrit[i].dualX() > Rcrit[j].dualX()) {
        rCount[setSize+i]++;
        if(aboveness(Lcrit[i].dualX(), Lcrit[i].dualY(), dualLines[j])) {
          belowCount[setSize+i]++;
        }
        if(aboveness(Rcrit[j].dualX(), Rcrit[j].dualY(), dualLines[setSize+i])) {
          aboveCount[setSize+i]++;
        } 
      }
    }
  }
  for(int i=0; i<setSize*2; i++)
    printf("%d ", rCount[i]); 
  printf("\n");
  for(int i=0; i<setSize*2; i++)
    printf("%d ", belowCount[i]); 
  printf("\n");
  for(int i=0; i<setSize*2; i++)
    printf("%d ", aboveCount[i]); 
  printf("\n");
}
*/
void findRayIntersections(float* origTimes, float* times, float* slopes, Line* dualLines, Ray* rRays, Ray* bRays, WORDTYPE** bitVec, int** prefixSums, WORDTYPE** slopeBitVec, int** slopePrefixSums, int offset, Peak* counts, int setSize, int parThreads, int setLevels, timeAndVals<float,float>* combined) {
//printf("parThreads:%d\n", parThreads);
  int i,j;
//  int setLevels = (int)ceil(log2((float)(setSize)))+1;
// Pointers to locations in each array for offset of current work
  WORDTYPE** bitPtrs = (WORDTYPE**)malloc(setLevels*sizeof(WORDTYPE*));
  int** prefixPtrs = (int**)malloc(setLevels*sizeof(int*));
//  WORDTYPE** slopeBitPtrs = (WORDTYPE**)malloc(setLevels*sizeof(WORDTYPE*));
//  int** slopePrefixPtrs = (int**)malloc(setLevels*sizeof(int*));

  for(i=0; i<setLevels; i++) {
    bitPtrs[i] = &bitVec[i][offset/WORDSIZE];
    prefixPtrs[i] = &prefixSums[i][offset/WORDSIZE];
//    slopeBitPtrs[i] = &slopeBitVec[i][offset/WORDSIZE];
//    slopePrefixPtrs[i] = &slopePrefixSums[i][offset/WORDSIZE];
  }

// Pointers to locations of blue trees (setSize past each red tree)
  WORDTYPE** bbitPtrs = (WORDTYPE**)malloc(setLevels*sizeof(WORDTYPE*));
  int** bprefixPtrs = (int**)malloc(setLevels*sizeof(int*));
//  WORDTYPE** bslopeBitPtrs = (WORDTYPE**)malloc(setLevels*sizeof(WORDTYPE*));
//  int** bslopePrefixPtrs = (int**)malloc(setLevels*sizeof(int*));

  for(i=0; i<setLevels; i++) {
    bbitPtrs[i] = &bitVec[i][(offset+setSize)/WORDSIZE];
    bprefixPtrs[i] = &prefixSums[i][(offset+setSize)/WORDSIZE];
//    bslopeBitPtrs[i] = &slopeBitVec[i][(offset+setSize)/WORDSIZE];
//    bslopePrefixPtrs[i] = &slopePrefixSums[i][(offset+setSize)/WORDSIZE];
  }
/*
for(int i=0; i<setSize; i++) {
  printf("(%.2f, %.2f, %.2f)", dualLines[i].slope, dualLines[i].icp, times[i]); 
}
printf("\n");
*/

//  int loopThreads = std::min(parThreads,setLevels);

//parThreads=1;
// Initialize nodeMem for BOTH red and blue trees
// Build RED TREE
  buildChazelTree<float, cmpTime>(setSize, bitPtrs, 0,times, prefixPtrs, parThreads);

//printTree<float>(origTimes, times, bitPtrs, setSize);


// Build BLUE TREE 
  buildChazelTree<float, cmpTimeDown>(setSize, bbitPtrs, 0,times+setSize, bprefixPtrs, parThreads);

/*
   for(int i=0; i<setSize; i++) {
  printf("(%.2f, %.2f, %.2f)", dualLines[i].slope, dualLines[i].icp, times[i]); 
}
printf("\n\n");
*/
//printf("blue points\n");

// Find RED rays below BLUE points
#pragma omp parallel private(i) shared(counts,setSize,setLevels,rRays,bRays) num_threads(parThreads) 
{
  #pragma omp for
  for(i=0; i<setSize; i++) {
    counts[setSize+i].below += seqQueryUnder<float,Line,cmpTime,aboveness>(setSize, setLevels, bRays[i].dualX(), bRays[i].dualY(), dualLines, times, bitPtrs, prefixPtrs, &counts[setSize+i].relCount, &counts[setSize+i].pi);
  }
/*
timeAndLine* lineSet = (timeAndLine*)malloc(setSize*sizeof(timeAndLine));
for(int i=0; i<setSize; i++) {
  lineSet[i].line = dualLines[i+setSize];
  lineSet[i].time = times[i+setSize];
}

std::sort(lineSet, lineSet+setSize, cmpLineStructDown);

for(int i=0; i<setSize; i++) {
  dualLines[i+setSize] = lineSet[i].line;
  times[i+setSize] = lineSet[i].time;
}
*/
// Find BLUE rays below RED points
  #pragma omp for
  for(i=0; i<setSize; i++) {
    counts[i].below += seqQueryUnder<float,Line,cmpTimeDown,aboveness>(setSize, setLevels, rRays[i].dualX(), rRays[i].dualY(), dualLines+setSize, times+setSize, bbitPtrs, bprefixPtrs, &counts[i].relCount, &counts[i].pi);
  }
}

#ifdef DEBUG
// Check that relCount is correctly computed for all setSize*2 points
bool correct=true;
int bruteVal=0;
for(int i=0; i<setSize; i++) {
bruteVal=0;
  for(int j=0; j<setSize; j++) {
    if(rRays[i].dualX() <= times[j+setSize]) bruteVal++;
  }
  if(bruteVal != counts[i].relCount) {
    correct=false;
  }
bruteVal=0;
  for(int j=0; j<setSize; j++) {
    if(bRays[i].dualX() >= times[j]) bruteVal++;
  }
  if(bruteVal != counts[i+setSize].relCount) {
    correct=false;
  }
}
if(!correct) printf("relCount incorrect!\n");

// Check that below is correctly computed
aboveness aboveCmp;
correct=true;
/*
for(int i=0; i<setSize; i++) {
  bruteVal=0;
  for(int j=0; j<setSize; j++) {
    if((rRays[i].dualX() <= times[j+setSize]) && aboveCmp(rRays[i].dualX(),rRays[i].dualY(),dualLines[j+setSize])) {
      printf("below:%d\n", j);
      bruteVal++;
    }
  }
  if(bruteVal != counts[i].below) {
    printf("bruteVal:%d, counts[%d]:%d\n", bruteVal, i, counts[i].below);
    printf("X:%.2f, Y:%.2f\n", rRays[i].dualX(), rRays[i].dualY());
    for(int i=0; i<setSize; i++) {
      printf("(%.2f + %.2f, %.2f) ", dualLines[setSize+i].slope, dualLines[setSize+i].icp, times[setSize+i]); 
    }
    correct=false;
  }
  bruteVal=0;
  for(int j=0; j<setSize; j++) {
    if((bRays[i].dualX() >= times[j]) && aboveCmp(bRays[i].dualX(), bRays[i].dualY(), dualLines[j]))
      bruteVal++;
  }
  if(bruteVal != counts[i+setSize].below) {
    printf("bruteVal:%d, counts[%d]:%d\n", bruteVal, i+setSize, counts[i+setSize].below);
    correct=false;
  }
}
*/
for(int i=0; i<setSize; i++) {
  bruteVal=0;
  for(int j=0; j<setSize; j++) {
    if((bRays[i].dualX() >= origTimes[j]) && aboveCmp(bRays[i].dualX(), bRays[i].dualY(), dualLines[j])) {
//      printf("below:%d\n", j);
      bruteVal++;
    }
  }
  if(bruteVal != counts[i+setSize].below) {
    printf("i:%d, bruteVal:%d, counts:%d\n", i, bruteVal, counts[i+setSize].below);
    printf("X:%.2f, Y:%.2f\n", bRays[i].dualX(), bRays[i].dualY());
    for(int idx=0; idx<setSize; idx++) {
      printf("(%.2f + %.2f, %.2f) ", dualLines[idx].slope, dualLines[idx].icp, origTimes[idx]); 
    }
    correct=false;
  }
}
if(!correct) printf("below incorrect!\n");
#endif


#pragma omp parallel for private(i) shared(combined, origTimes, counts) num_threads(parThreads) 
for(i=0; i<setSize*2; i++) {
  combined[i].time = origTimes[i];
//  combined[i].val = counts[i].pi;
  combined[i].val = dualLines[i].slope;
}



#pragma omp parallel for private(i) shared(combined, times,slopes) num_threads(parThreads) 
for(i=0; i<setSize*2; i++) {
  times[i] = combined[i].time;
  slopes[i] = combined[i].val;
}

#ifdef DEBUG
// Make sure sets properly sorted
#endif


// Build tree using structure so that vals are sorted.
//buildChazelTree<float, cmpTime>(setSize, slopeBitPtrs, 0,times, slopePrefixPtrs, parThreads);
buildChazelTree<float, cmpTime>(setSize, bitPtrs, 0,times, prefixPtrs, parThreads);

// Build tree using structure so that vals are sorted.
//buildChazelTree<float, cmpTimeDown>(setSize, bslopeBitPtrs, 0,times+setSize, bslopePrefixPtrs, parThreads);
buildChazelTree<float, cmpTimeDown>(setSize, bbitPtrs, 0,times+setSize, bprefixPtrs, parThreads);


#pragma omp parallel private(i) shared(counts,setSize,setLevels,rRays,bRays) num_threads(parThreads) 
{
  #pragma omp for
  for(i=0; i<setSize; i++) {
//    counts[i+setSize].above += seqQueryAbove<float,float,cmpTime,cmpPiDown>(setSize, setLevels, bRays[i].dualX(), bRays[i].dualSlope(), slopes, times, slopeBitPtrs, slopePrefixPtrs);
    counts[i+setSize].above += seqQueryAbove<float,float,cmpTime,cmpPiDown>(setSize, setLevels, bRays[i].dualX(), bRays[i].dualSlope(), slopes, times, bitPtrs, prefixPtrs);
  }

  #pragma omp for
  for(i=0; i<setSize; i++) {
//    counts[i].above += seqQueryAbove<float,float,cmpTimeDown,cmpPi>(setSize, setLevels, rRays[i].dualX(), rRays[i].dualSlope(), slopes+setSize, times+setSize, bslopeBitPtrs, bslopePrefixPtrs);
    counts[i].above += seqQueryAbove<float,float,cmpTimeDown,cmpPi>(setSize, setLevels, rRays[i].dualX(), rRays[i].dualSlope(), slopes+setSize, times+setSize, bbitPtrs, bprefixPtrs);
  } 
}

#ifdef DEBUG
// Check that above is correct for all points
#endif

  free(bitPtrs);
//  free(slopeBitPtrs);
  free(prefixPtrs);
  //  free(slopePrefixPtrs);
  free(bbitPtrs);
//  free(bslopeBitPtrs);
  free(bprefixPtrs);
//  free(bslopePrefixPtrs);

//  free(combined);
/*
  printf("Red SlopeTree\n");
  for(int i=0; i<setLevels; i++) {
    for(int j=0; j<setSize; j++) {
      printf("%f ", rPiRoot[i][j].val);
    }
  printf("\n");
  }

  printf("Blue SlopeTree\n");
  for(int i=0; i<setLevels; i++) {
    for(int j=0; j<setSize; j++) {
      printf("%f ", bPiRoot[i][j].val);
    }
  printf("\n");
  }

printf("Red Rays\n");
for(int i=0; i<setSize; i++) {
  printf("x:%f, y:%f, slope:%f, icp:%f\n", rRays[i].dualX(), rRays[i].dualY(), rRays[i].dualSlope(), rRays[i].dualIcp());
}

printf("Blue Rays\n");
for(int i=0; i<setSize; i++) {
  printf("x:%f, y:%f, slope:%f, icp:%f\n", bRays[i].dualX(), bRays[i].dualY(), bRays[i].dualSlope(), bRays[i].dualIcp());
}
for(int i=0; i<setSize*2; i++) 
  printf("%d ", counts[i].relCount);
printf("\n");
for(int i=0; i<setSize*2; i++) 
  printf("%d ", counts[i].below);
printf("\n");
for(int i=0; i<setSize*2; i++) 
  printf("%d ", counts[i].above);
printf("\n");
*/

}


void vis1D(int* xvals, int* peaks, int* visIdx, int SIZE, int numThreads) {
  omp_set_nested(1); // 1 - enables nested parallelism; 0 - disables nested parallelism.
  int retval;
  std::chrono::steady_clock::time_point beginSeq;
  std::chrono::steady_clock::time_point endSeq;
  int preTime=0;

//  beginSeq= std::chrono::steady_clock::now();
// Prepare structure for chazel sturcture
  int levels = (int)ceil(log2((float)(SIZE)))+1;

  timeAndVals<float,float>* combined = (timeAndVals<float,float>*)malloc(SIZE * sizeof(timeAndVals<float,float>));

  float* origTimes = (float*)malloc(SIZE*sizeof(float));
  float* times = (float*)malloc(SIZE*sizeof(float));
  float* slopeTimes = (float*)malloc(SIZE*sizeof(float));

  WORDTYPE** bitVec = (WORDTYPE**)malloc(levels*sizeof(WORDTYPE*));
  int** prefixSums = (int**)malloc(levels*sizeof(int*));
  WORDTYPE** slopeBitVec = (WORDTYPE**)malloc(levels*sizeof(WORDTYPE*));
  int** slopePrefixSums = (int**)malloc(levels*sizeof(int*));
  for(int i=0; i<levels; i++) {
    bitVec[i] = (WORDTYPE*)malloc((SIZE/WORDSIZE)*sizeof(WORDTYPE));
    prefixSums[i] = (int*)malloc((SIZE/WORDSIZE)*sizeof(int));
    slopeBitVec[i] = (WORDTYPE*)malloc((SIZE/WORDSIZE)*sizeof(WORDTYPE));
    slopePrefixSums[i] = (int*)malloc((SIZE/WORDSIZE)*sizeof(int));

    for(int j=0; j<SIZE/WORDSIZE; j++) {
      bitVec[i][j]=0;
      prefixSums[i][j]=0;
      slopeBitVec[i][j]=0;
      slopePrefixSums[i][j]=0;
    }

  }

  
  Ray* Rcrit = (Ray*)malloc(SIZE*sizeof(Ray));
  Ray* Lcrit = (Ray*)malloc(SIZE*sizeof(Ray));
  Peak* counts = (Peak*)malloc(SIZE*sizeof(Peak));
  for(int i=0; i<SIZE; i++) {
    counts[i].relCount=0;
    counts[i].below=0;
    counts[i].above=0;
    counts[i].pi = 0;
  }

//  Hull** convexHulls = (Hull**)malloc(numThreads*sizeof(Hull*));
  Hull** convexHulls[2];
  bool hullBit=false;

  Point* points = (Point*)malloc(SIZE*sizeof(Point));
  for(int i=0; i<SIZE; i++) {
    points[i].x=xvals[i];
    points[i].y=peaks[i];
  }

int i=0;
int threadsPerGroup=1;
#pragma omp parallel for private(i,threadsPerGroup) shared(xvals, peaks, visIdx, Rcrit, Lcrit) num_threads(numThreads) 
  for(i=0; i<SIZE; i+=WORDSIZE) {
    initCritRays(&points[i], &visIdx[i], &Rcrit[i], &Lcrit[i]);
  }
#ifdef DEBUG
// Check that initCritRays correctly calculates visibility for base case
#endif
  
  Line* dualLines = (Line*)malloc(SIZE*sizeof(Line));
int j=0;
  
  int initTime=0;
  int intersectionTime=0;
  int convexTime=0;
  int workGroups;
  int parGroups;

//  endSeq= std::chrono::steady_clock::now();
//  preTime += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();
//  convexHulls[1] = (Hull**)malloc(sizeof(Hull*));

  int setLevels;
  for(int step=WORDSIZE; step <= SIZE/2; step = step*2){ // Each merge step of overall algorithm
  setLevels = (int)ceil(log2((float)(step)))+1;
  workGroups = SIZE/(step*2);
  threadsPerGroup = std::max(numThreads/workGroups, 1);
  parGroups = std::min(numThreads, workGroups);

//  convexHulls[(int)hullBit] = (Hull**)malloc(workGroups*sizeof(Hull*));
#pragma omp parallel for private(i) shared(origTimes, Rcrit, dualLines, Lcrit, counts) num_threads(parGroups)
    for(i=0; i<=(SIZE-(step*2)); i+=(step*2)) { // Calc intersection of each pair of groups

#pragma omp parallel private(j) num_threads(threadsPerGroup)
{
//      beginSeq= std::chrono::steady_clock::now();
      #pragma omp for nowait
      for(j=0; j<step; j++) { // Use x-coordinates as times in AOT
        origTimes[i+j] = Rcrit[i+j].dualX();
        times[i+j] = origTimes[i+j];
        slopeTimes[i+j] = origTimes[i+j];
        dualLines[i+j].slope = Rcrit[i+j].dualSlope();
        dualLines[i+j].icp = Rcrit[i+j].dualIcp();
      }
      #pragma omp for nowait
      for(j=step; j<(step*2); j++) { // Use x-coordinates as times in AOT
        origTimes[i+j] = Lcrit[i+j].dualX();
        times[i+j] = origTimes[i+j];
        slopeTimes[i+j] = origTimes[i+j];
        dualLines[i+j].slope = Lcrit[i+j].dualSlope();
        dualLines[i+j].icp = Lcrit[i+j].dualIcp();
      }
}
      int offset=i;
//      endSeq= std::chrono::steady_clock::now();
//      initTime += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();

//      beginSeq= std::chrono::steady_clock::now();
      findRayIntersections(&origTimes[i], &times[i], &slopeTimes[i], &dualLines[i], &Rcrit[i], &Lcrit[i+step], bitVec, prefixSums, slopeBitVec, slopePrefixSums, offset, &counts[i], step,threadsPerGroup, setLevels, combined+i);
//      endSeq= std::chrono::steady_clock::now();
//      intersectionTime += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();
  

// Update Critical Rays
//      beginSeq= std::chrono::steady_clock::now();
//      convexHulls[i] = updateCriticalRaysPar(&Rcrit[i], &Lcrit[i+step], convexHulls[i], convexHulls[i+step], step, threadsPerGroup);

/*
      if(workGroups == numThreads) {
        convexHulls[(int)hullBit][i/(step*2)] = updateCriticalRaysReturnHull(&Rcrit[i], &Lcrit[i+step], step, 1);
      } else if(workGroups < numThreads) {
        convexHulls[(int)hullBit][i/(step*2)] = updateCriticalRaysPar(&Rcrit[i], &Lcrit[i+step], convexHulls[(int)!hullBit][(i/(step*2))*2], convexHulls[(int)!hullBit][(i/(step*2))*2+1], step, threadsPerGroup);
//        free(convexHulls[(int)!hullBit]);
      } else {
        updateCriticalRays(&Rcrit[i], &Lcrit[i+step], step, threadsPerGroup);
      }
*/
    updateCriticalRays(&Rcrit[i], &Lcrit[i+step], step, threadsPerGroup);
////      endSeq= std::chrono::steady_clock::now();
//     convexTime += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();
    }

//        free(convexHulls[(int)!hullBit]);
//        hullBit=!hullBit;
#ifdef DEBUG
// Check that this round properly calculated intersection
#endif
  }

//printf("pretime:%d, init:%d, intersection:%d, convex:%d\n", preTime, initTime, intersectionTime, convexTime);
/*
for(int i=0; i<SIZE; i++)
  printf("%d ", counts[i].relCount);
printf("\n");
for(int i=0; i<SIZE; i++)
  printf("%d ", counts[i].below);
printf("\n");
for(int i=0; i<SIZE; i++)
  printf("%d ", counts[i].above);
printf("\n");
*/

  for(int i=0; i<SIZE; i++) {
    visIdx[i] += (counts[i].below+counts[i].above)-counts[i].relCount;
  }

  free(origTimes);
  free(slopeTimes);
  for(int i=0; i<levels; i++) {
    free(bitVec[i]);
    free(prefixSums[i]);
    free(slopeBitVec[i]);
    free(slopePrefixSums[i]);
  }
  free(bitVec);
  free(prefixSums);
  free(slopeBitVec);
  free(slopePrefixSums);
  free(combined);
  free(Rcrit);
  free(Lcrit);
  free(counts);
  free(points);
  free(dualLines);
} 

  
