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


#include <climits>
#include <math.h>
#include<algorithm>
#include<queue>
#include<pthread.h>
#include<unistd.h>

//#define NUM_THREADS 8
//#define SEQ_LEVELS 5
//#define PAR 2

//#define malloc mymalloc

void* mymalloc(int size) {
  printf("allocating %d\n", size);
  return malloc(size);
}

int done;
template<typename T>
struct pthreadMerge {
  bool (*key)(T,T);
  bool (*time)(float,float);
  int size;
  int id;
  int levels;
  T* vals;
  AOTNode<T>** nodes;
};

template<typename T>
struct mergePathData {
  int id;
  AOTNode<T>* A;
  AOTNode<T>* B;
  int aSize;
  int bSize;
  AOTNode<T>* result;
  bool (*key)(T,T);
  bool (*time)(float,float);
  T* vals;
  int numThreads;
};


template<typename T>
void seqMergeNodes(AOTNode<T>* a, AOTNode<T>* b, int aSize, int bSize, T* minVal, AOTNode<T>* outMem, bool (*key)(T,T), bool (*time)(float,float)) {
  int leftIdx=0;
  int rightIdx=0;
  int i=0;
  bool leftItem;
  AOTNode<T>* mergeNode;
  T minLeft = minVal[0];
  T minRight = minVal[aSize];

i=0;
// Step through each node group
  while(leftIdx < aSize && rightIdx < bSize) {
    if((*time)(a[leftIdx].t, b[rightIdx].t)) { // Merge from left set
      mergeNode = &a[leftIdx];
      if(rightIdx==0) { // Special case where right pointer goes to left set
        outMem[i].initializeNode(minLeft, mergeNode->t, NULL, mergeNode);
      }

      else {
        outMem[i].initializeNode(minRight, mergeNode->t, mergeNode, &b[rightIdx-1]);
      }
      leftIdx++;
    } 
    else { // Merge from right set
      mergeNode = &b[rightIdx];
      if(leftIdx==0) {
        outMem[i].initializeNode(minRight, mergeNode->t, NULL, mergeNode);
      }
      else {
        outMem[i].initializeNode(minRight, mergeNode->t, &a[leftIdx-1], mergeNode);
      }
      rightIdx++;
    }
    i++;
  } 

  // Finish up remain list
  while(leftIdx < aSize) {
    mergeNode = &a[leftIdx];
//printf("%d, %d\n", i, leftIdx);
    outMem[i].initializeNode(minRight, mergeNode->t, mergeNode, &b[rightIdx-1]);
    leftIdx++;
    i++;
  }
  while(rightIdx < bSize) {
    mergeNode = &b[rightIdx];
    outMem[i].initializeNode(minRight, mergeNode->t, &a[leftIdx-1], mergeNode);
    rightIdx++;
    i++;
  }
}

// Expects minVal to be a 2-elm array with the first value from left and right set, respectively.
template<typename T>
void seqMergeBackwards(AOTNode<T>* a, AOTNode<T>* b, int aSize, int bSize, int workSize, T* minVal, AOTNode<T>* outMem, bool (*key)(T,T), bool (*time)(float,float)) {
  int leftIdx=aSize-1;
  int rightIdx=bSize-1;
  int i=0;
  bool leftItem;
  AOTNode<T>* mergeNode;
  T minLeft = minVal[0];
  T minRight = minVal[1];

  for(i=workSize-1; i>=0; i--) {
//printf("i:%d, a[%d]:%d, b[%d]:%d\n", i, leftIdx, (int)a[leftIdx].t, rightIdx, (int)b[rightIdx].t);
    if((rightIdx < 0 || !time(a[leftIdx].t,b[rightIdx].t)) && leftIdx >=0) {
      mergeNode = &a[leftIdx];
//printf("A writing %d to location:%d (%04x)\n", (int)mergeNode->t, i, &outMem[i]);
// Step through each node group
      if(rightIdx<0) { // Special case where right pointer goes to left set
        outMem[i].initializeNode(minLeft, mergeNode->t, NULL, mergeNode);
      }
      else {
        outMem[i].initializeNode(minRight, mergeNode->t, mergeNode, &b[rightIdx]);
      }
      leftIdx--;
    } 
    else { // Merge from right set
      mergeNode = &b[rightIdx];
//printf("B writing %d to location:%d (%04x)\n", (int)mergeNode->t, i, &outMem[i]);
      if(leftIdx<0) {
        outMem[i].initializeNode(minRight, mergeNode->t, NULL, mergeNode);
      }
      else {
        outMem[i].initializeNode(minRight, mergeNode->t, &a[leftIdx], mergeNode);
      }
      rightIdx--;
    }
//    i--;
  } 
/*
  // Finish up remain list
  while(leftIdx < aSize) {
    mergeNode = &a[leftIdx];
//printf("%d, %d\n", i, leftIdx);
    outMem[i].initializeNode(minRight, mergeNode->t, mergeNode, &b[rightIdx-1]);
    leftIdx++;
    i++;
  }
  while(rightIdx < bSize) {
    mergeNode = &b[rightIdx];
    outMem[i].initializeNode(minRight, mergeNode->t, &a[leftIdx-1], mergeNode);
    rightIdx++;
    i++;
  }
*/
}


template<typename T>
int* findPivot(AOTNode<T>* A, AOTNode<T>* B, int aSize, int bSize, int* idx, int workSize, bool (*time)(float,float)) {
//  int* idx = (int*)malloc(2*sizeof(int));  
//  idx[0] = (int)ceil((float)aSize/2.0);
//  idx[1] = (int)ceil((float)bSize/2.0);
//  for(int step=(int)ceil((float)std::min(aSize,bSize)/4.0);step >= 1; step = step/2) {
int id = idx[0];
  for(int step=workSize/2; step >= 1; step = step/2) {

    if(time(A[idx[0]].t,B[idx[1]].t)) {
      idx[0] += step;
      idx[1] -= step;
    }
    else {
      idx[0] -= step;
      idx[1] += step;
    }
  }

  if(idx[0] == 0 || idx[1] == bSize) {
    if(time(A[idx[0]].t,B[idx[1]-1].t)) {
      idx[0]++;
      idx[1]--;
    }
  }
  else if(idx[1] == 0 || idx[0] == aSize) {
    if(time(B[idx[1]].t,A[idx[0]-1].t)) {
      idx[0]--;
      idx[1]++;
    }
  } 
  else if(time(B[idx[1]].t,A[idx[0]-1].t) && time(B[idx[1]-1].t,A[idx[0]].t)) {
    if(idx[1] < bSize) {
      idx[0]--;
      idx[1]++;
    }
  }
  else if(!time(B[idx[1]].t,A[idx[0]+1].t) && !time(B[idx[1]-1].t,A[idx[0]].t)) {
    if(idx[0] < aSize) {
      idx[0]++;
      idx[1]--;
    }
  }
/*
  for(int i=0; i<2; i++) {
    if(time(A[idx[0]].t,B[idx[1]].t) && idx[1] >0) {
      idx[0]++;
      idx[1]--;
      if(!time(A[idx[0]].t,B[idx[1]].t)) {
        idx[0]--;
        idx[1]++;
      }
    }
    if(!time(A[idx[0]].t,B[idx[1]].t) && idx[0] > 0) {
      idx[0]--;
      idx[1]++;
      if(time(A[idx[0]].t, B[idx[1]].t)) {
        idx[0]++;
        idx[1]--;
      }
    }
  }
*/
  return idx;
}

template<typename T>
void* threadMergePath(void* arguments) {
  mergePathData<T>* args = (mergePathData<T>*)arguments;
//printf("Allocating pivots %d\n", 2*sizeof(int));
  int* pivots = (int*)malloc(2*sizeof(int));
  int workSize = (args->aSize+args->bSize)/args->numThreads;
  T minVal[2]; // Smallest values from left and right subtrees
  minVal[0] = args->vals[0];
  minVal[1] = args->vals[args->aSize];
  pivots[0] = args->id*workSize/2;
  pivots[1] = args->id*workSize/2;
  int offset = (args->id-1)*workSize;
 
//  if(args->id*2 <= NUM_THREADS) {
  if(args->id < args->numThreads) {
    if(args->id > args->numThreads/2) 
      pivots = findPivot(args->A, args->B, args->aSize, args->bSize, pivots, std::min(workSize,args->aSize-pivots[0]), args->time);
    else
      pivots = findPivot(args->A, args->B, args->aSize, args->bSize, pivots, workSize, args->time);
//    printf("thread:%d, workSize:%d, pivots[0]:%d, pivots[1]:%d\n", args->id, workSize, pivots[0],pivots[1]);
//printf("writing to position starting at:%d, size:%d\n", aOffset+bOffset, workSize);
    seqMergeBackwards(args->A, args->B, pivots[0], pivots[1], workSize, minVal, &args->result[offset], args->key, args->time);
//  }

//  else if(args->id < NUM_THREADS){
//    args->result = &args->result[(args->id-1) * workSize];

//    pivots = findPivot(&args->A[(args->aSize-aSizeNew)], &args->B[args->bSize-bSizeNew], aSizeNew, bSizeNew, pivots, args->time);
//    pivots = findPivot(args->A, args->B, aSizeNew, bSizeNew, pivots, workSize, args->time);
  }
  else {
//    printf("thread:%d, workSize:%d, writing to:%d\n", args->id, workSize, args->aSize+args->bSize-workSize);
    seqMergeBackwards(args->A, args->B, args->aSize, args->bSize, workSize, minVal, &args->result[args->aSize+args->bSize-workSize], args->key, args->time);
  }
//printf("Freeing pivots %d\n", sizeof(pivots));
  //free(pivots);
}

template<typename T>
void mergepathMergeNodes(AOTNode<T>* a, AOTNode<T>* b, int aSize, int bSize, T* minVal, AOTNode<T>* outMem, bool (*key)(T,T), bool (*time)(float,float), int numThreads) {

//printf("aSize:%d, bSize:%d\n", aSize, bSize);
//  pthread_t threads[numThreads];
//printf("Allocating threads %d\n", numThreads*sizeof(pthread_t));
  pthread_t* threads = (pthread_t*)malloc(numThreads*sizeof(pthread_t));
//printf("Allocating threadArgs %d\n", numThreads*sizeof(mergePathData<T>));
  mergePathData<T>* threadArgs = (mergePathData<T>*)malloc(numThreads*sizeof(mergePathData<T>));
  // Set up data to send to parallel function
  for(int i=0; i<(numThreads); i++) {
    threadArgs[i].id = i+1;
    threadArgs[i].A = a;
    threadArgs[i].B = b;
    threadArgs[i].result = outMem;
    threadArgs[i].aSize = aSize;
    threadArgs[i].bSize = bSize;
    threadArgs[i].key = key;
    threadArgs[i].time = time;
    threadArgs[i].vals = minVal;
    threadArgs[i].numThreads = numThreads;

   // threadMergePath((void*)&threadArgs[i]);
  }
  for(int i=1; i<numThreads; i++) {
    pthread_create(&threads[i], NULL, threadMergePath<T>, (void *)&threadArgs[i]);
  }

  threadMergePath<T>((void*)&threadArgs[0]);

  for(int i=1; i<numThreads; i++)
    pthread_join(threads[i],NULL);
//printf("Freeing:%d\n",sizeof(threadArgs));
  //free(threadArgs);
//printf("Freeing:%d\n",sizeof(threads));
  //free(threads); 
}

template<typename T>
void* parMerge(void* arguments) {

  pthreadMerge<T>* args = (pthreadMerge<T>*)arguments;
  int offset = args->size*args->id;
  int step=1;
//printf("thread:%d, parallel size:%d, offset:%d\n", args->id, args->size, offset);
  for(int i=0; i<args->levels; i++) {
    for(int j=offset; j<(args->size+offset); j+=(2*step)) {
//printf("Thread:%d merging %d-%d with %d-%d\n", args->id, j,j+step-1, j+step, j+step+step-1);
//printf("thread:%d at level:%d, merging %d and %d\n", args->id, i, j, j+step);
      seqMergeNodes<T>(&args->nodes[i][j], &args->nodes[i][j+step], step, step, &args->vals[j], &args->nodes[i+1][j], args->key, args->time);
    }
    step = step *2;
  }
  if(args->id == 1) {
//printf("setting done to true!\n");
    done = done + 1;
  }
//  return NULL;
//  if((int)arguments->id > 0)
//    pthread_exit(NULL); 

return NULL;
}


template<typename T>
AOTNode<T>** buildAOT(T* vals, float* times, int size, int levels, bool (*f)(T,T), bool (*timeCmp)(float,float), int seqLevels, int parallel, int numThreads, AOTNode<T>** nodes) {
  int i,j;
  int rc;
//printf("levels:%d, size:%d\n", levels, size);
//for(i=0; i<size; i++) {
//printf("%f.0 ", vals[i].slope);
//}
//printf("\n");
//for(i=0; i<size; i++) {
//printf("%f.0 ", times[i]);
//}
//printf("\n");
//printf("Allocating nodePtrs %d\n", levels*sizeof(AOTNode<T>*));
//  AOTNode<T>** nodes = (AOTNode<T>**)malloc(levels*sizeof(AOTNode<T>*));
//  DualRay** minVals = (DualRay**)malloc(size*sizeof(DualRay*));
 
//printf("levels:%d, sizeofAOT:%d\n", levels, sizeof(AOTNode<T>));
//  for(i=0; i<levels; i++) {
//printf("Allocating nodes %d\n", size*sizeof(AOTNode<T>));
//    nodes[i] = (AOTNode<T>*)malloc(size*sizeof(AOTNode<T>));
//  }
/*
printf("\nBuilding tree\n");
for(int i=0; i<size; i++) {
printf("%f ", vals[i]);
}
printf("\n");
for(int i=0; i<size; i++) {
printf("%f ", times[i]);
}
*/
  int parLevels=0;
  int remLevels;
  int step=1;

//  for(i=0; i<levels; i++) {
//    for(j=0; j<size;j++) {
//      nodes[i][j].initializeNode(vals[j], times[j]);
//    }
//  }

//  pthread_t threads[NUM_THREADS];
//printf("Allocating threads %d\n", numThreads*sizeof(pthread_t));
  pthread_t* threads = (pthread_t*)malloc(numThreads*sizeof(pthread_t));
if(parallel==1) {
  if(seqLevels >= levels) {
//printf("building tree size:%d, levels:%d, parLevels:0\n", size, levels);
    remLevels = levels;
    parLevels = 0;
    step = 1;
  }
  else {
    parLevels = (levels - seqLevels);
//printf("building tree size:%d, levels:%d, parLevels:%d, threads:%d\n", size, levels, parLevels, numThreads);
    done=0;
//    parLevels = 1;
    remLevels = seqLevels;
//    remLevels = levels - 1;
//printf("Allocating threadArgs %d\n", numThreads*sizeof(pthreadMerge<T>));
    pthreadMerge<T>* threadArgs = (pthreadMerge<T>*)malloc(numThreads*sizeof(pthreadMerge<T>));
    for(i=0; i<numThreads; i++) {
      threadArgs[i].key = f;
      threadArgs[i].time = timeCmp;
      threadArgs[i].size = size/numThreads;
      threadArgs[i].id = i;
      threadArgs[i].vals = vals;
      threadArgs[i].nodes = nodes;
      threadArgs[i].levels = parLevels;
    }
    for(i=1; i<numThreads; i++) {
      rc = pthread_create(&threads[i], NULL, parMerge<T>, (void *)&threadArgs[i]);
//      pthread_create(&threads[i], NULL, parMerge<T>, NULL);
   }
    parMerge<T>((void *)&threadArgs[0]);

    for(i=1; i<numThreads; i++)
      pthread_join(threads[i],NULL);

//printf("Freeing threadArgs %d\n", sizeof(threadArgs));
//    free(threadArgs);
//printf("Freeing threads %d\n", sizeof(threads));
//    free(threads);

    step= 1 << (parLevels);
//step=1;
//    free(threadArgs);
  }
}
int temp=1;
//printf("done with parallel!\n");
// Finish top levels sequentially

if(parallel == 0) {
  for(i=parLevels; i<levels-1; i++) {
//  for(i=0; i<levels-1; i++) {
    for(j=0; j<size; j+=(2*step)) {
//      for(int k=0; k<size*size*size; k++)
//        temp+=i;
      seqMergeNodes<T>(&nodes[i][j], &nodes[i][j+step], step, step, &vals[j], &nodes[i+1][j], f, timeCmp);
    }
    step = step*2;
  }
}
else {
  for(i=parLevels; i<levels-1; i++) {
//  for(i=0; i<levels-1; i++) {
    for(j=0; j<size; j+=(2*step)) {
//      for(int k=0; k<size*size*size; k++)
//        temp+=i;
      mergepathMergeNodes<T>(&nodes[i][j], &nodes[i][j+step], step, step, &vals[j], &nodes[i+1][j], f, timeCmp, numThreads);
    }
    step = step*2;
  }
}


// Finish top levels with mergepath merging
/*
  for(i=parLevels; i<levels-1; i++) {
  for(i=0; i<levels-1; i++) {

if(parallel == 1) 
      mergepathMergeNodes<T>(&nodes[i][j], &nodes[i][j+step], step, step, &vals[j], &nodes[i+1][j], f, timeCmp, numThreads);
else
      seqMergeNodes<T>(&nodes[i][j], &nodes[i][j+step], step, step, &vals[j], &nodes[i+1][j], f, timeCmp);
    }
    step = step*2;
  }
*/
//    for(j=0; j<size; j+=(2*step)) {
//      for(int k=0; k<size*size*size; k++)
//        temp+=i;
//printTree(nodes,8);
/*
printf("A=[");
for(int k=0; k<step; k++) {
  printf("%d, ", (int)nodes[i][j+k].t);
}
printf("]\n");
printf("B=[");
for(int k=0; k<step; k++) {
  printf("%d, ", (int)nodes[i][j+step+k].t);
}
printf("]\n");
*/
/*
printf("Result=[");
for(int k=0; k<step*2; k++) {
//  printf("%d (%04x), ", (int)nodes[i+1][j+k].t, &nodes[i+1][j+k]);
  printf("%d, ", (int)nodes[i+1][j+k].t);
}
printf("]\n");
*/

  return nodes;
//  return NULL;

}



template<typename T>
AOTNode<T>* searchTime(AOTNode<T>* group, int size, double time, bool (*timeCmp)(float,float)) {
  int position=size/2;
  int step=size/4;

  if((*timeCmp)(time, group[0].t) || time == group[0].t) {
    return NULL;
  }

  while(step > 0) {
    if((*timeCmp)(group[position].t, time))
      position += step;
    else
      position -= step;
    step = step/2;
  }
  if(!(*timeCmp)(group[position].t, time))
    position--;

  return &group[position];
}

template<typename T>
AOTNode<T>* searchTimeEq(AOTNode<T>* group, int size,double time, bool (*timeCmp)(float,float)) {
  int position=size/2;
  int step=size/4;

  if((*timeCmp)(time, group[0].t)) {
    return NULL;
  }

  while(step > 0) {
    if((*timeCmp)(group[position].t, time))
      position += step;
    else
      position -= step;
    step = step/2;
  }
  if(!(*timeCmp)(group[position].t, time))
    position--;

  return &group[position];
}



template<typename T>
int underCount(AOTNode<T>* rootNode, int size, int levels, float qX, float qY, float* pi, bool (*cmp)(float,float,T)) {

  int count=0;
  int depth=0;
//  if(rootNode != NULL) {
  if(rootNode != NULL) {
    while(depth < levels-1) {
//printf("qX:%f, qY:%f, slope:%f, icp:%f, result:%d\n", qX, qY, rootNode->val.slope, rootNode->val.icp, (*cmp)(qX, qY, rootNode->val));
      if((*cmp)(qX, qY, rootNode->val)) { // right branch
        if(rootNode->left != NULL) {
//printf("Going right, leftSize:%d\n", rootNode->left->subTreeSize);
          count += rootNode->left->subTreeSize;
        }
        rootNode = rootNode->right;
      }
      else { // left branch
//printf("Going left, leftNull?:%d\n", (rootNode->left == NULL));
        if(rootNode->left != NULL)
          rootNode = rootNode->left;
        else
          break;
      }
      depth++;
    }

//printf("Last Node - qX:%f, qY:%f, slope:%f, icp:%f, result:%d\n", qX, qY, rootNode->val.slope, rootNode->val.icp, (*cmp)(qX, qY, rootNode->val));
    if((*cmp)(qX, qY, rootNode->val)) { // Check that above at least 1
      *pi = rootNode->val.slope; 
      count++;
    }
  }
  return count;
}

template<typename T>
int aboveCount(AOTNode<T>* rootNode, int size, int levels, float qSlope, bool (*cmp)(float,float)) {

  int count=0;
  int depth=0;
  if(rootNode != NULL) {
    while(depth < levels-1) {
//printf("qSlope:%f, nodeSlope:%f, result:%d\n", qSlope, rootNode->val, (*cmp)(qSlope, rootNode->val));
      if((*cmp)(qSlope, rootNode->val)) { // right branch
        if(rootNode->left != NULL)
          count += rootNode->left->subTreeSize;
        rootNode = rootNode->right;
      }
      else { // left branch
        if(rootNode->left != NULL)
          rootNode = rootNode->left;
        else
          break;
      }
      depth++;
    }
//printf("Leaf - qSlope:%f, nodeSlope:%f, result:%d\n", qSlope, rootNode->val, (*cmp)(qSlope, rootNode->val));
    if((*cmp)(qSlope,rootNode->val)) { // Check that above at least 1
      count++;
    }
  }
  return count;
}

/*
int AboveCount(AOTNode** AOTRoot, int size, int levels, DualRay* qVal, bool (*cmp)(DualRay*,DualRay*), bool (*timeCmp)(double,double)) {

  AOTNode* rootNode;
//  double epsilon = (-0.00001);

//  if((*timeCmp)(-0.0001, 0.0001))
//    epsilon += 0.00002;

  rootNode = searchTime(AOTRoot[levels-1], size,((double)qVal->getA()), timeCmp);

//printf("b%d, slope:%lf, rootTime:%lf\n", qVal->label, (double)qVal->getA(), rootNode->t);
  int count=0;
  if(rootNode != NULL && rootNode->t != INT_MIN && rootNode->t != INT_MAX) {
//  printf("b%d, slope:%d, rootNodeTime:%lf\n", qVal->label, qVal->getA(), rootNode->t);
    while(rootNode->right != NULL) {
      if((*cmp)(qVal, rootNode->val)) { // right branch
//printf("true! right branch\n");
        if(rootNode->left != NULL)
          count += rootNode->left->subTreeSize;
        rootNode = rootNode->right;
      }
      else { // left branch
//printf("false! left branch\n");
        if(rootNode->left != NULL)
          rootNode = rootNode->left;
        else
          break;
      }
    }
//printf("At leaf, comparing %lf and %lf, cmp:%d\n", qVal->getEndX(), rootNode->val->getEndX(), (*cmp)(qVal, rootNode->val));
    if((*cmp)(qVal, rootNode->val)) { // Check that above at least 1
      count++;
    }
  }
  return count;
}
*/


template<typename T>
void printTree(AOTNode<T>** root, int size) {
  int levels = (int)log2((float)size)+1;
printf("in printTree.  levels:%d\n", levels);
  for(int i=levels-1; i>=0; i--) {
    for(int j=0; j<size; j++) {
      printf("(t:%d, k:%d ", (int)root[i][j].t, (int)root[i][j].val);
      if(root[i][j].left != NULL)
        printf("L");
      if(root[i][j].right != NULL)
        printf("R");
      printf(") - ");
    }
    printf("\n");
  }

}

