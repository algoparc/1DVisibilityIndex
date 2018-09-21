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


#include "SimpleRBI.h"
#include "convexhull.h"
#include<stdio.h>
#include<chrono>
#include<stdlib.h>
//#include<adds.h>
#include<algorithm>
#include<cassert>
#include <ext/pb_ds/tag_and_trait.hpp>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/hash_policy.hpp>
#include <utility>

using namespace std;
using namespace __gnu_pbds;

float calcSlope(int x1, int y1, int x2, int y2) {
  int deltaX = x2-x1;
  int deltaY = y2-y1;
 
  return (float)deltaY/(float)deltaX; 
}
float calcIcp(int x, int y, float slope) {
  return (float)y-((float)x*slope) ;
}

void calcAbove(DualRay** red, DualRay** blue, int leftsize, int rightsize);

inline bool aboveness(DualRay* a, DualRay* b) {
  double aY = (double)a->getEndY();
  double bY = (((double)b->getA())*((double)a->getEndX())) + (double)b->getB();
  return (aY > bY);
//  double dist = aY - bY;
//printf("aY:%lf, bY:%lf, dist:%lf\n", aY, bY,dist);
//  return dist >0;
}

inline bool cmpEndXUp(DualRay* a, DualRay* b) {
  return a->getEndX() < b->getEndX();
}

inline bool cmpEndXDown(DualRay* a, DualRay* b) {
  return a->getEndX() > b->getEndX();
}

inline float cmpEndXDistDown(DualRay* a, DualRay* b) {
  return (float)b->getEndX() - (float)a->getEndX();
}

inline float cmpEndXDistUp(DualRay* a, DualRay* b) {
  return (float)a->getEndX() - (float)b->getEndX();
}

inline bool cmpSlopeUp(DualRay* a, DualRay* b) {
//  return(a->getA() < b->getA());
  return (a->a < b->a);
}
inline bool cmpSlopeDown(DualRay* a, DualRay* b) {
//  return(a->getA() > b->getA());
  return(a->a > b->a);
}

inline bool cmpImmBelowUp(DualRay* a, DualRay* b) {
  if(a->immBelow == NULL)
    return false;
  else if(b->immBelow == NULL)
    return true;
  else {
//    printf("%d, %d\n", a->immBelow->getA(), b->immBelow->getA());
    return (a->immBelow->getA() < b->immBelow->getA()); 
  }
}

inline bool cmpImmBelowDown(DualRay* a, DualRay* b) {
  if(a->immBelow == NULL)
    return false;
  else if(b->immBelow == NULL)
    return true;
  else
    return (a->immBelow->getA() > b->immBelow->getA()); 
}

inline bool cmpTimeUp(double a, double b) {
  return a < b;
}

inline bool cmpTimeDown(double a, double b) {
  return b < a;
}

inline bool cmpTimeEqUp(double a, double b) {
  return a <= b;
}

inline bool cmpTimeEqDown(double a, double b) {
  return b <= a;
}

void initCritRays(int* xvals, int* peaks, int* visIdx, Ray* Rcrit, Ray* Lcrit) {
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


/******************************************************************/
/*   PLANESWEEP SEQUENTIAL METHOD                                 */
/******************************************************************/

//EITHER CONVERT DATA IN SWEEP TO DUALRAYS OR CONVERT THIS TO USE MY FORMAT INSTEAD OF DUALRAYS!
void calcSimpleRBI(DualRay** red, int leftsize, DualRay** blue, int rightsize, PNode* redData, PNode* blueData) {
  int i,j;

/**************************************************
*  Perform planesweep with RED rays and BLUE points
**************************************************/

  std::sort(red, red+leftsize, cmpEndXUp);
  std::sort(blue, blue+rightsize, cmpEndXUp);

for(i=0; i<rightsize; i++) {
//  printf("b%d: slope:%d, intercept:%d, endX:%lf, endY:%lf\n", i, blue[i]->getA(), blue[i]->getB(), blue[i]->getEndX(), blue[i]->getEndY());
  blue[i]->Under = 0;
  blue[i]->Above = 0;
  blue[i]->Overlap=0;
  blue[i]->immBelow=NULL;
}

//printf("NEWR ED TREE\n");
  PNode* redTree=NULL;
  i=0;
  j=0;
//  PNode* redData = (PNode*)malloc((leftsize)*sizeof(PNode));
  while(i<leftsize && j<rightsize) {
    if(red[i]->getEndX() < blue[j]->getEndX()) {
      if(redTree == NULL) {
        redTree = &redData[0];
        redTree->initialize(red[i]);
      }
//        redTree = new PNode(red[i]);
      else

        insertPNode(red[i], redTree, &redData[i], cmpSlopeUp);
      i++;
//printf("red Tree!\n");
//printTree(redTree);
    }
    else {
      blue[j]->Overlap = i;
      if(i>0) {
        blue[j]->Under = underCount(blue[j], redTree);
      }
      j++;
    }
  }
  while(j<rightsize) {
    blue[j]->Overlap = i;
    if(i>0) {
      blue[j]->Under = underCount(blue[j], redTree);
    }
    j++;
  }
//  free(data);

/**************************************************
*  Perform planesweep with BLUE rays and RED points
**************************************************/
// Blue tree is built on decreasing slope and based on decreasing x-value

for(i=0; i<rightsize; i++) {
//  printf("b%d: slope:%d, intercept:%d, endX:%lf, endY:%lf\n", i, blue[i]->getA(), blue[i]->getB(), blue[i]->getEndX(), blue[i]->getEndY());
  red[i]->Under = 0;
  red[i]->Above = 0;
  red[i]->Overlap=0;
  red[i]->immBelow=NULL;
}

//  PNode* blueData = (PNode*)malloc(rightsize*sizeof(PNode));
  int dataIdx=0;
  PNode* blueTree=NULL;
  i=rightsize-1;
  j=leftsize-1;
//printf("NEW BLUE TREE\n");
  while(i>=0 && j>=0) {
//printf("r%d, i:%d, rightsize-1:%d\n", red[j]->label, i, rightsize-1);
    if(blue[i]->getEndX() > red[j]->getEndX()) {
      if(blueTree == NULL) {
        blueTree = &blueData[0];
        blueTree->initialize(blue[i]);
      }
//        blueTree = new PNode(blue[i]);
      else
        insertPNode(blue[i], blueTree, &blueData[dataIdx], cmpSlopeDown);
//printf("blue tree!\n");
//printTree(blueTree);
      i--;
      dataIdx++;
    }
    else {
      red[j]->Overlap = (rightsize-1)-i;
      if(red[j]->Overlap > 0) {
        red[j]->Under = underCount(red[j], blueTree);
      }
      j--;
    }
  }
  while(j>=0) {
    red[j]->Overlap = (rightsize-1)-i;
    if(red[j]->Overlap > 0) {
      red[j]->Under = underCount(red[j], blueTree);
    }
    j--;
  }
//  free(redData);
//  free(blueData);

//  deletePTree(redTree);
//  deletePTree(blueTree);

  calcAbove(red, blue, leftsize, rightsize);

/*
  for(i=0; i<leftsize; i++) {
    red[i]->getCell()->addToVisibilityCountRight(red[i]->Above+red[i]->Under-red[i]->Overlap);
  }
  for(i=0; i<rightsize; i++) {
    blue[i]->getCell()->addToVisibilityCountLeft(blue[i]->Above+blue[i]->Under-blue[i]->Overlap);
  }
*/
}




void calcAbove(DualRay** red, DualRay** blue, int leftsize, int rightsize) {
  int i,j;
  std::sort(red, red+leftsize, cmpSlopeDown);
  std::sort(blue, blue+rightsize, cmpImmBelowDown);

  PNode* blueData = (PNode*)malloc(rightsize*sizeof(PNode));
  
  PNode* blueSlopeTree=NULL;
  i=0;
  j=0;
  while(i<leftsize && j < rightsize && blue[j]->immBelow != NULL) {
    if(blueSlopeTree == NULL) {
      if(blue[j]->immBelow->getA() >= red[i]->getA()) {
//        blueSlopeTree = new PNode(blue[j]);
        blueSlopeTree = &blueData[0];
        blueSlopeTree->initialize(blue[j]);
        j++;
      }
      else
        i++;
    }
    else if(red[i]->getA() > blue[j]->immBelow->getA()) { // Search Red ray
//      printTree(blueSlopeTree);
      red[i]->Above = aboveCount(red[i], blueSlopeTree, cmpEndXUp);
      i++;
    } 
    else {
      insertPNodeDups(blue[j], blueSlopeTree, &blueData[j], cmpEndXDistUp);
      j++;
    }
  }
  while(i<leftsize && blueSlopeTree != NULL) {
//      printTree(blueSlopeTree);
    red[i]->Above = aboveCount(red[i], blueSlopeTree, cmpEndXUp);
    i++;
  }
//  deletePTree(blueSlopeTree);
  free(blueData);

  std::sort(blue, blue+rightsize, cmpSlopeUp);
  std::sort(red, red+leftsize, cmpImmBelowUp);

  PNode* redData = (PNode*)malloc(leftsize*sizeof(PNode));
  PNode* redSlopeTree=NULL;
  i=0;
  j=0;
  while(i<rightsize && j < leftsize && red[j]->immBelow != NULL) {
    if(redSlopeTree == NULL) {
      if(red[j]->immBelow->getA() <= blue[i]->getA()) {
//        redSlopeTree = new PNode(red[j]);
        redSlopeTree = &redData[0];
        redSlopeTree->initialize(red[j]);
        j++;
      }
      else
        i++;
    }
    else if(blue[i]->getA() < red[j]->immBelow->getA()) { // Search Red ray
//      printTree(blueSlopeTree);
      blue[i]->Above = aboveCount(blue[i], redSlopeTree, cmpEndXDown);
      i++;
    } 
    else {
      insertPNodeDups(red[j], redSlopeTree, &redData[j], cmpEndXDistDown);
      j++;
    }
  }
  while(i<rightsize && redSlopeTree != NULL) {
//      printTree(blueSlopeTree);
    blue[i]->Above = aboveCount(blue[i], redSlopeTree, cmpEndXDown);
    i++;
  }
//  deletePTree(redSlopeTree);
  free(redData);

}

void sweepVis(int* xvals, int* peaks, int* visIdx, int SIZE) {


//  std::chrono::steady_clock::time_point beginSeq;
//  std::chrono::steady_clock::time_point endSeq;

//  beginSeq= std::chrono::steady_clock::now();
// Prepare structure for chazel sturcture
  int levels = (int)ceil(log2((float)(SIZE)))+1;

  float* initTimes = (float*)malloc(SIZE*sizeof(float));

  Ray* Rcrit = (Ray*)malloc(SIZE*sizeof(Ray));
  Ray* Lcrit = (Ray*)malloc(SIZE*sizeof(Ray));
  Peak* counts = (Peak*)malloc(SIZE*sizeof(Peak));
  DualRay** dualRays = (DualRay**)malloc(SIZE*sizeof(DualRay*));
  DualRay* dualData = (DualRay*)malloc(SIZE*sizeof(DualRay));
  PNode* treeData = (PNode*)malloc(SIZE*sizeof(PNode));
  for(int i=0; i<SIZE; i++) {
    counts[i].relCount=0;
    counts[i].below=0;
    counts[i].above=0;
    counts[i].pi = 0;
    dualRays[i] = &dualData[i];
  }
int i=0;
//printf("numThreads:%d\n", numThreads);
  for(i=0; i<SIZE; i=i+4) {
    initCritRays(&xvals[i], &peaks[i], &visIdx[i], &Rcrit[i], &Lcrit[i]);
  }
  
  int step=4;

//  Line* dualLines = (Line*)malloc(SIZE*sizeof(Line));
//  AOTNode<Line>** memPtr = (AOTNode<Line>**)malloc(levels*sizeof(AOTNode<Line>*));
//  AOTNode<float>** slopeMemPtr = (AOTNode<float>**)malloc(levels*sizeof(AOTNode<float>*));
int j=0;
  
  int initTime=0;
  int intersectionTime=0;
  int convexTime=0;
  int workGroups;
  int parGroups;
  DualRay* redPtr;
  DualRay* bluePtr;

//  endSeq= std::chrono::steady_clock::now();
//  preTime += std::chrono::duration_cast<std::chrono::milliseconds>(endSeq - beginSeq).count();


  for(step=4; step <= SIZE/2; step = step*2){ // Each merge step of overall algorithm
    for(i=0; i<=(SIZE-(step*2)); i+=(step*2)) { // Calc intersection of each pair of groups

      for(j=0; j<step; j++) { // Use x-coordinates as times in AOT
        dualRays[i+j]->endX = Rcrit[i+j].dualX(); 
        dualRays[i+j]->setA(Rcrit[i+j].dualSlope());
//        dualLines[i+j].icp = Rcrit[i+j].dualIcp();
      }
      for(j=step; j<(step*2); j++) { // Use x-coordinates as times in AOT
        dualRays[i+j]->endX = Lcrit[i+j].dualX(); 
        dualRays[i+j]->setA(Lcrit[i+j].dualSlope());
//        dualLines[i+j].icp = Rcrit[i+j].dualIcp();
      }

      int offset=i;
/*
      redPtr = &dualRays[i];
      bluePtr = &dualRays[i+step];
for(int i=0; i<step; i++) 
  printf("red[%d]=%f\n", i, redPtr[i].getEndX());
for(int i=0; i<step; i++) 
  printf("blue[%d]=%f\n", i, bluePtr[i].getEndX());
*/
      calcSimpleRBI(&dualRays[i], step, &dualRays[i+step], step, &treeData[i], &treeData[i+step]);
  

// Update Critical Rays
      updateCriticalRays(&Rcrit[i], &Lcrit[i+step], step, 1);
    }
  }


//  print_counters(EventSet);

  for(int i=0; i<SIZE; i++) {
    visIdx[i] += (counts[i].below+counts[i].above)-counts[i].relCount;
  }

  free(initTimes);
  free(Lcrit);
  free(Rcrit);
  free(counts);
  free(dualData);
  free(dualRays);
  free(treeData);
} 

  
