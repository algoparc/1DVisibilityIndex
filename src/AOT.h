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

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>

template<typename T>
class AOTNode {
public:
  AOTNode<T>* left;
  AOTNode<T>* right;
  float t;
  int subTreeSize;
  T val;
  int height;
//  int minSlope;

//template<typename T>
  AOTNode<T>(T inVal, float inT) {
    left=NULL;
    right=NULL;
    t=inT;
    val=inVal;
    subTreeSize=1;
    height=0;
  }
    
  AOTNode<T>() {
    left=NULL;
    right=NULL;
    t=0.0;
    val=NULL;
    subTreeSize=1;
  }

  ~AOTNode<T>() {
  }

//  template<typename T>
  AOTNode(T inVal, float inT, AOTNode<T>* leftPtr, AOTNode<T>* rightPtr) {
    left=leftPtr;
    right=rightPtr;
    t=inT;
    val=inVal;
    subTreeSize=0; // Internal nodes don't count
    height=0;
    if(leftPtr != NULL) {
      subTreeSize+=leftPtr->subTreeSize;
      height = leftPtr->height+1;
    }
    if(rightPtr != NULL) {
      subTreeSize+=rightPtr->subTreeSize;
      height = rightPtr->height+1;
    }
  }
//  template<typename T>
  void initializeNode(T inVal, float inT, AOTNode<T>* leftPtr, AOTNode<T>* rightPtr) {
    left=leftPtr;
    right=rightPtr;
    t=inT;
    val=inVal;
    subTreeSize=0; // Internal nodes don't count
    height=0;
    if(leftPtr != NULL) {
      subTreeSize+=leftPtr->subTreeSize;
      height = leftPtr->height+1;
    }
    if(rightPtr != NULL) {
      subTreeSize+=rightPtr->subTreeSize;
      height = rightPtr->height+1;
    }
  }
//  template<typename T>
  void initializeNode(T inVal, float inT) {
    left=NULL;
    right=NULL;
    t=inT;
    val=inVal;
    subTreeSize=1;
    height=0;
  }
};
/*
template<typename T>
AOTNode<T>** buildAOT(T* vals, double* times, int size, int levels, bool (*f)(int,int), bool (*timeCmp)(double,double));


//int underCount(AOTNode** root, int size, int levels, DualRay* qVal, bool (*cmp)(DualRay*,DualRay*), bool (*timeCmp)(double,double));

template<typename T>
void printTree(AOTNode<T>** root, int size);

template<typename T>
AOTNode<T>* searchTime(AOTNode<T>* group, int size, double time, bool (*timeCmp)(double,double));

//int AboveCount(AOTNode** AOTRoot, int size, int levels, DualRay* qVal, bool (*cmp)(DualRay*,DualRay*), bool (*timeCmp)(double,double));

template<typename T>
void deleteAOTNode(AOTNode<T>* node);

template<typename T>
void seqMergeNodes(AOTNode<T>* a, AOTNode<T>* b, int aSize, int bSize, T** minVal, AOTNode<T>* outMem, bool (*key)(T*,T*), bool (*time)(double,double));
*/

#include "AOT-utils.h"
