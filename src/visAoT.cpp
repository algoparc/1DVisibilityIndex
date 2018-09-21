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


#include "AOT.h"
#include "convexhull.h"
#include <chrono>
#include<stdlib.h>
//#include <jemalloc/jemalloc.h>
#include<omp.h>



/*
void print_counters(int EventSet) {
  long long values[10];
  int retval;
  retval=PAPI_read(EventSet,values);
  if(retval!=PAPI_OK)
    printf("read_counters failed!\n");

//printf("L1 misses:%lld\n", values[0]);
  printf("Total cycles: %lld\n", values[0]);
  printf("CYCLES LDM:%lld (%f)\n", values[1], (float)values[1]/(float)values[0]);
  printf("CYCLES NO EXECUTE:%lld (%f)\n", values[2], (float)values[2]/(float)values[0]);
  printf("STALLS LDM:%lld\n", values[3]);
  printf("Cache invalids:%lld\n", values[4]);
  printf("%lld %lld %lld %lld %lld\n", values[0], values[1], values[2], values[3], values[4]);


}
*/

inline bool cmpSlopes(Line a, Line b) {
  return(a.slope <= b.slope);
}

inline bool cmpSlopesDown(Line a, Line b) {
  return(a.slope >= b.slope);
}

inline bool cmpPi(float a, float b) {
  return (a <= b);
}

inline bool cmpPiDown(float a, float b) {
  return (a >= b);
}

inline bool cmpTime(float a, float b) {
  return(a < b);
}
inline bool cmpTimeDown(float a, float b) {
  return(a > b);
}

inline bool aboveness(float x, float y, Line line) {
  return (y > ((line.icp)+(x*line.slope)));
}

float calcSlope(int x1, int y1, int x2, int y2) {
  int deltaX = x2-x1;
  int deltaY = y2-y1;
 
  return (float)deltaY/(float)deltaX; 
}
float calcIcp(int x, int y, float slope) {
  return (float)y-((float)x*slope) ;
}

/*
void initCritRays(int* xvals, int* peaks, int* visIdx, Ray* Rcrit, Ray* Lcrit) {
  for(int i=0; i<4; i++) {
    Rcrit[i].x=1;
    Lcrit[i].x=1;
    Rcrit[i].y=1;
    Lcrit[i].y=1;
    Lcrit[i].slope=1;
    Rcrit[i].slope=1;
    Lcrit[i].icp=1;
    Rcrit[i].icp=1;
  }
}
*/
// Initialize critical rays and visIdx for groups of 4
void initCritRaysAoT(int* xvals, int* peaks, int* visIdx, Ray* Rcrit, Ray* Lcrit) {
  float tempTangent;
  int j;
  // Initilize critical rays
  for(int i=0; i<4; i++) {
    Rcrit[i].x = xvals[i];
    Lcrit[i].x = xvals[i];
    Rcrit[i].y = peaks[i];
    Lcrit[i].y = peaks[i];
    visIdx[i]=1; // All rays initially see themselves
  }
  Lcrit[0].slope = 9999.9;
  Lcrit[0].icp = 9999.9;
  Rcrit[4].slope = -9999.9;
  Rcrit[4].icp = -9999.9;
  for(int i=0; i<3; i++) {
    tempTangent = calcSlope(xvals[i], peaks[i], xvals[i+1], peaks[i+1]);
    visIdx[i]++;
    for(j=i+2; j<4; j++) {
      if(tempTangent < calcSlope(xvals[i], peaks[i], xvals[j], peaks[j])) {
        visIdx[i]++;
        tempTangent = calcSlope(xvals[i], peaks[i], xvals[j], peaks[j]);
      }
    }
    Rcrit[i].slope = tempTangent;
    Rcrit[i].icp = calcIcp(xvals[i], peaks[i], tempTangent);
  }

  for(int i=3; i>0; i--) {
    tempTangent = calcSlope(xvals[i], peaks[i], xvals[i-1], peaks[i-1]);
    visIdx[i]++;
    for(j=i-2; j>=0; j--) {
      if(tempTangent > calcSlope(xvals[i], peaks[i], xvals[j], peaks[j])) {
        visIdx[i]++;
        tempTangent = calcSlope(xvals[i], peaks[i], xvals[j], peaks[j]);
      }
    }
    Lcrit[i].slope = tempTangent;
    Lcrit[i].icp = calcIcp(xvals[i], peaks[i], tempTangent);
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
void findRayIntersectionsAoT(float* initTimes, Line* dualLines, Ray* rRays, Ray* bRays, AOTNode<Line>** nodeMem, AOTNode<float>** slopeNodeMem, int offset, Peak* counts, int setSize, int parThreads) {
//printf("parThreads:%d\n", parThreads);
  int i,j;
  int setLevels = (int)ceil(log2((float)(setSize)))+1;
  AOTNode<Line>* rootPtr;


  int loopThreads = std::min(parThreads,setLevels);
  AOTNode<Line>** memPtr = (AOTNode<Line>**)malloc(setSize*sizeof(AOTNode<Line>*));
// Initialize nodeMem for BOTH red and blue trees
//#pragma omp parallel for private(i,j) shared(memPtr,nodeMem,offset,dualLines,initTimes) num_threads(loopThreads)
  for(i=0; i<setLevels; i++) {
    memPtr[i] = &nodeMem[i][offset];
#pragma omp parallel for private(i,j) shared(memPtr,dualLines,initTimes) num_threads(parThreads)
    for(j=0; j<setSize*2;j++) {
      memPtr[i][j].initializeNode(dualLines[j], initTimes[j]);
    }
  }
// Build RED TREE
  AOTNode<Line>** rRoot;
  if(parThreads == 1) {
    rRoot = buildAOT<Line>(dualLines, initTimes, setSize, setLevels, cmpSlopes, cmpTime, 999, 0, 1,memPtr);
  }
  else {
    rRoot = buildAOT<Line>(dualLines, initTimes, setSize, setLevels, cmpSlopes, cmpTime, 4, 1, parThreads ,memPtr);
  }


  AOTNode<Line>** bMemPtr = (AOTNode<Line>**)malloc(setLevels*sizeof(AOTNode<Line>*));
  for(int i=0; i<setLevels; i++)
    bMemPtr[i] = &memPtr[i][setSize]; // Move nodeMem to location for blue tree

  AOTNode<Line>** bRoot;
  if(parThreads == 1) {
    bRoot = buildAOT<Line>(&dualLines[setSize], &initTimes[setSize], setSize, setLevels, cmpSlopesDown, cmpTimeDown, 999, 0, 1, bMemPtr);
  }
  else {
    bRoot = buildAOT<Line>(&dualLines[setSize], &initTimes[setSize], setSize, setLevels, cmpSlopesDown, cmpTimeDown, 4, 1, parThreads, bMemPtr);
  }

// Find RED rays below BLUE points
#pragma omp parallel private(i,rootPtr) shared(rRoot,counts,setSize,setLevels,rRays,bRoot,bRays) num_threads(parThreads) 
{
//#pragma omp parallel for private(i,rootPtr) shared(rRoot,counts,setSize,setLevels,bRays) num_threads(parThreads)
  #pragma omp for
  for(i=0; i<setSize; i++) {
    rootPtr = searchTimeEq<Line>(rRoot[setLevels-1], setSize, bRays[i].dualX(), cmpTime);
    if(rootPtr != NULL) {// In case it has time < all
//printf("qTime:%f, root:%f\n", bRays[i].dualX(), rootPtr->t);
      counts[setSize+i].relCount += rootPtr->subTreeSize;
      counts[setSize+i].below += underCount<Line>(rootPtr, setSize, setLevels, bRays[i].dualX(), bRays[i].dualY(), &counts[setSize+i].pi, aboveness);
//printf("relCount:%d, below:%d\n", counts[setSize+i].relCount, counts[setSize+i].below);
    }
  }

// Find BLUE rays below RED points
//#pragma omp parallel for private(i,rootPtr) shared(bRoot,counts,setSize,setLevels,rRays) num_threads(parThreads)
  #pragma omp for
  for(i=0; i<setSize; i++) {
    rootPtr = searchTimeEq<Line>(bRoot[setLevels-1], setSize, rRays[i].dualX(), cmpTimeDown);
    if(rootPtr != NULL) { // In case it has time < all
//printf("qTime:%f, root:%f\n", bRays[i].dualX(), rootPtr->t);
      counts[i].relCount += rootPtr->subTreeSize;
      counts[i].below += underCount<Line>(rootPtr, setSize, setLevels, rRays[i].dualX(), rRays[i].dualY(), &counts[i].pi, aboveness);
//printf("relCount:%d, below:%d\n", counts[i].relCount, counts[i].below);
    }
  }
}
/*
int bruteVal;
int correct=true;
for(int i=0; i<setSize; i++) {
  bruteVal=0;
  for(int j=0; j<setSize; j++) {
    if((bRays[i].dualX() >= initTimes[j]) && aboveness(bRays[i].dualX(), bRays[i].dualY(), dualLines[j])) {
      printf("below:%d\n", j);
      bruteVal++;
    }
  }
  if(bruteVal != counts[i+setSize].below) {
    printf("i:%d, bruteVal:%d, counts:%d\n", i, bruteVal, counts[i+setSize].below);
    printf("X:%.2f, Y:%.2f\n", bRays[i].dualX(), bRays[i].dualY());
    for(int idx=0; idx<setSize; idx++) {
      printf("(%.2f + %.2f, %.2f) ", dualLines[idx].slope, dualLines[idx].icp, initTimes[idx]);
    }
    correct=false;
  }
}
if(!correct) printf("below incorrect!\n");
*/

// MAKE NEW TREES.  TIME=X, BST ON PI
// search ray into pi-tree to get above.
  AOTNode<float>** slopeMemPtr = (AOTNode<float>**)malloc(setSize*sizeof(AOTNode<float>*));
//#pragma omp parallel for private(i,j) shared(slopeMemPtr,slopeNodeMem,offset,counts,initTimes) num_threads(loopThreads)
  for(i=0; i<setLevels; i++) {
    slopeMemPtr[i] = &slopeNodeMem[i][offset];
#pragma omp parallel for private(i,j) shared(memPtr,dualLines,initTimes) num_threads(parThreads)
    for(j=0; j<setSize*2;j++) {
      slopeMemPtr[i][j].initializeNode(counts[j].pi, initTimes[j]);
    }
  }
  float* piVals = (float*)malloc(setSize*2*sizeof(float));
  for(i=0; i<setSize*2; i++) {
    piVals[i] = counts[i].pi;
  }
/*
printf("initTimes\n");
for(int i=0; i<setSize*2; i++) {
  printf("%f ", initTimes[i]);
}
printf("\n");

printf("\n");
for(int i=0; i<2*setSize; i++ ){
  printf("%f ", piVals[i]);
}
printf("\n");
 */   
// Build RED TREE

  AOTNode<float>** rPiRoot;
  if(parThreads == 1) {
    rPiRoot = buildAOT<float>(piVals, initTimes, setSize, setLevels, cmpPiDown, cmpTime, 999, 0, 1, slopeMemPtr);
  }
  else {
    rPiRoot = buildAOT<float>(piVals, initTimes, setSize, setLevels, cmpPiDown, cmpTime, 5, 1, parThreads, slopeMemPtr);
  }
  AOTNode<float>* piRootPtr = NULL;


  AOTNode<float>** bSlopeMemPtr = (AOTNode<float>**)malloc(setLevels*sizeof(AOTNode<float>*));
  for(i=0; i<setLevels; i++)
    bSlopeMemPtr[i] = &slopeMemPtr[i][setSize]; // Move nodeMem to location for blue tree

  AOTNode<float>** bPiRoot;
  if(parThreads == 1) {
    bPiRoot = buildAOT<float>(&piVals[setSize], &initTimes[setSize], setSize, setLevels, cmpPi, cmpTimeDown, 999, 0, 1, bSlopeMemPtr);
  }
  else {
    bPiRoot = buildAOT<float>(&piVals[setSize], &initTimes[setSize], setSize, setLevels, cmpPi, cmpTimeDown, 5, 1, parThreads, bSlopeMemPtr);
  }
  piRootPtr = NULL;
/*
printf("\ntimes:\n");
for(int i=0; i<setSize; i++) {
  printf("%f ", rPiRoot[setLevels-1][i].t);
}
printf("\nvals\n");
for(int i=0; i<setSize; i++) {
  printf("%f ", rPiRoot[setLevels-1][i].val);
}
*/
//#pragma omp parallel for private(i,piRootPtr) shared(rPiRoot,counts,setSize,setLevels,bRays) num_threads(parThreads)
#pragma omp parallel private(i,piRootPtr) shared(rPiRoot,counts,setSize,setLevels,bRays,bPiRoot,rRays) num_threads(parThreads) 
{
  #pragma omp for
  for(i=0; i<setSize; i++) {
    piRootPtr = searchTimeEq<float>(rPiRoot[setLevels-1], setSize, bRays[i].dualX(), cmpTime);
    if(piRootPtr != NULL) {
//      printf("qTime:%f.0, root:%f.0\n", bRays[i].dualX(), piRootPtr->t);
      counts[setSize+i].above += aboveCount<float>(piRootPtr, setSize, setLevels, bRays[i].dualSlope(), cmpPiDown);
    }
  }

//#pragma omp parallel for private(i,piRootPtr) shared(bPiRoot,counts,setSize,setLevels,rRays) num_threads(parThreads)
  #pragma omp for
  for(i=0; i<setSize; i++) {
    piRootPtr = searchTimeEq<float>(bPiRoot[setLevels-1], setSize, rRays[i].dualX(), cmpTimeDown);
    if(piRootPtr != NULL) {
//      printf("qTime:%f.0, root:%f.0\n", rRays[i].dualX(), piRootPtr->t);
      counts[i].above += aboveCount<float>(piRootPtr, setSize, setLevels, rRays[i].dualSlope(), cmpPi);
    }
  }
}

  free(memPtr);
  free(bMemPtr);
  free(slopeMemPtr);
  free(piVals);
  free(bSlopeMemPtr);

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


void vis1DAoT(int* xvals, int* peaks, int* visIdx, int SIZE, int numThreads) {
  int retval;
//  std::chrono::steady_clock::time_point beginSeq;
//  std::chrono::steady_clock::time_point endSeq;

//  beginSeq= std::chrono::steady_clock::now();
// Prepare structure for chazel sturcture
  int levels = (int)ceil(log2((float)(SIZE)))+1;

  float* initTimes = (float*)malloc(SIZE*sizeof(float));
  AOTNode<Line>** nodeMem = (AOTNode<Line>**)malloc(levels*sizeof(AOTNode<Line>*));
  AOTNode<float>** slopeNodeMem = (AOTNode<float>**)malloc(levels*sizeof(AOTNode<float>*));
  for(int i=0; i<levels; i++) {
    nodeMem[i] = (AOTNode<Line>*)malloc(SIZE*sizeof(AOTNode<Line>));
    slopeNodeMem[i] = (AOTNode<float>*)malloc(SIZE*sizeof(AOTNode<float>));
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
int i=0;
int threadsPerGroup=1;
//printf("numThreads:%d\n", numThreads);
#pragma omp parallel for private(i,threadsPerGroup) shared(xvals, peaks, visIdx, Rcrit, Lcrit) num_threads(numThreads) 
  for(i=0; i<SIZE; i=i+4) {
    initCritRaysAoT(&xvals[i], &peaks[i], &visIdx[i], &Rcrit[i], &Lcrit[i]);
  }
  
  int step=4;

  Line* dualLines = (Line*)malloc(SIZE*sizeof(Line));
//  AOTNode<Line>** memPtr = (AOTNode<Line>**)malloc(levels*sizeof(AOTNode<Line>*));
//  AOTNode<float>** slopeMemPtr = (AOTNode<float>**)malloc(levels*sizeof(AOTNode<float>*));
int j=0;
  
  int initTime=0;
  int intersectionTime=0;
  int convexTime=0;
  int workGroups;
  int parGroups;

//  endSeq= std::chrono::steady_clock::now();
//  preTime += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();

  for(step=4; step <= SIZE/2; step = step*2){ // Each merge step of overall algorithm
  workGroups = SIZE/(step*2);
  threadsPerGroup = std::max(numThreads/workGroups, 1);
//  while(threadsPerGroup < numThreads/workGroups)
//    threadsPerGroup *= 2;
  parGroups = std::min(numThreads, workGroups);
//printf("parGroups:%d\n", parGroups);
//  printf("workGroups:%d, threadsPerGroup:%d, numThreads:%d\n", workGroups, threadsPerGroup, numThreads);
#pragma omp parallel for private(i) shared(initTimes, Rcrit, dualLines, Lcrit, nodeMem, slopeNodeMem, counts) num_threads(parGroups)
    for(i=0; i<=(SIZE-(step*2)); i+=(step*2)) { // Calc intersection of each pair of groups

#pragma omp parallel private(j) num_threads(threadsPerGroup)
{
//      beginSeq= std::chrono::steady_clock::now();
//#pragma omp parallel for private(j) num_threads(threadsPerGroup)
      #pragma omp for nowait
      for(j=0; j<step; j++) { // Use x-coordinates as times in AOT
        initTimes[i+j] = Rcrit[i+j].dualX();
        dualLines[i+j].slope = Rcrit[i+j].dualSlope();
        dualLines[i+j].icp = Rcrit[i+j].dualIcp();
      }
//#pragma omp parallel for private(j) num_threads(threadsPerGroup)
      #pragma omp for nowait
      for(j=step; j<(step*2); j++) { // Use x-coordinates as times in AOT
        initTimes[i+j] = Lcrit[i+j].dualX();
        dualLines[i+j].slope = Lcrit[i+j].dualSlope();
        dualLines[i+j].icp = Lcrit[i+j].dualIcp();
      }
}
      int offset=i;

//      endSeq= std::chrono::steady_clock::now();
//      initTime += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();

//      beginSeq= std::chrono::steady_clock::now();
//      findRayIntersections(&initTimes[i], &dualLines[i], &Rcrit[i], &Lcrit[i+step], memPtr, slopeMemPtr, offset, &counts[i], step);
      findRayIntersectionsAoT(&initTimes[i], &dualLines[i], &Rcrit[i], &Lcrit[i+step], nodeMem, slopeNodeMem, offset, &counts[i], step,threadsPerGroup);
//      endSeq= std::chrono::steady_clock::now();
//      intersectionTime += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();
  

// Update Critical Rays
//      beginSeq= std::chrono::steady_clock::now();
      updateCriticalRays(&Rcrit[i], &Lcrit[i+step], step, threadsPerGroup);
//      updateCriticalRays(&Rcrit[i], &Lcrit[i+step], step, 1);
//      endSeq= std::chrono::steady_clock::now();
//      convexTime += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();
    }
  }

//printf("preLoop:%d, init:%d, intersection:%d, convex:%d\n", preTime, initTime, intersectionTime, convexTime);
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

//  print_counters(EventSet);

  for(int i=0; i<SIZE; i++) {
    visIdx[i] += (counts[i].below+counts[i].above)-counts[i].relCount;
  }

  free(initTimes);
  for(int i=0; i<levels; i++) {
    free(nodeMem[i]);
    free(slopeNodeMem[i]);
  }
  free(nodeMem);
  free(slopeNodeMem);
} 

  
