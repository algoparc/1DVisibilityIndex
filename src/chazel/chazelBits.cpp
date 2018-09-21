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


#include "chazelBits.h"
#include <omp.h>

/******************************************************************/
/* IN PROGRESS - persistent BST packing bits into WORDTYPE type   */
/******************************************************************/

template<typename T,typename cmp>
void mergePath(T* A, T* B, int size, int threads, int threadId, int* left, int* right) {
  int step, initPos; 
  cmp cmpTime;
  if(threadId == 0) {
    *left=0;
    *right=0;
/*
for(int i=0; i<size; i++) {
  printf("%d ", A[i]);
}
printf("\n");
for(int i=0; i<size; i++) {
  printf("%d ", B[i]);
}
printf("\n");
*/
  }  
  else {
    initPos=(threadId)*(size/threads);
    *left = initPos;
    *right = initPos;
    if(initPos > size/2) 
      step = (size-initPos)/2;
    else
      step = initPos/2;
    
    while(step >= 1) {
//      if(A[*left] < B[*right]) {
      if(cmpTime(A[*left], B[*right])) {
        *left += step;
        *right -= step;
      }
      else {
        *left -= step;
        *right += step;
      }
      step /= 2;
    }
//    if(*left > 0 && (A[*left-1] >= B[*right])) {
    if(*left > 0 && cmpTime(B[*right],A[*left-1])) {
      (*left)--;
      (*right)++;
    }
//    if(*right > 0 && (B[*right-1] > A[*left])) {
    if(*right > 0 && cmpTime(A[*left],B[*right-1])) {
      (*left)++;
      (*right)--;
    }   
  }
//    printf("thread:%d, left:%d, right:%d\n", threadId, *left, *right);
}

// If need to use multiple threads per merge
template<typename T, typename cmp>
void mergeLevelLarge(int N, T* input, T* output, WORDTYPE* bitVec, int step, int tasks, int numThreads) {
  int threadsPerTask = numThreads/tasks;
  omp_set_nested(1); // 1 - enables nested parallelism; 0 - disables nested parallelism.
//printf("step:%d\n", step);
  // Parallel by the number of merge tasks first
  #pragma omp parallel num_threads(tasks) shared(input, output, bitVec)
//for(int i=0; i<tasks; i++)
{
  int taskSize = step*2;
  int taskOffset = omp_get_thread_num() * taskSize;
//  int taskOffset = i * taskSize;
  int rightOffset = taskOffset+step;
// printf("task:%d, taskOffset:%d, rightOffset:%d\n", omp_get_thread_num(), taskOffset, rightOffset); 
  #pragma omp parallel num_threads(threadsPerTask) shared(input, output, bitVec)
  {
    int leftIdx, rightIdx;
    leftIdx=0;
    rightIdx=0;
//    mergePath<T,cmp>(input+taskOffset, input+rightOffset, step, threadsPerTask, omp_get_thread_num(), &leftIdx, &rightIdx);

    cmp cmpTime;
    int writeOffset = taskOffset + (omp_get_thread_num()*(taskSize/threadsPerTask));
    int writeIdx;
    int wordOffset = writeOffset/WORDSIZE;
    int bitOffset=0;
//printf("thread:%d, writeOffset:%d\n", omp_get_thread_num(), writeOffset);
    // Each thread merges X elements total
    for(writeIdx=0; writeIdx<(taskSize/threadsPerTask); writeIdx++) {
      bitVec[wordOffset] = bitVec[wordOffset] << 1; // shift by 1 place and put a 0 in LSB
//      if(input[taskOffset+leftIdx] < input[rightOffset+rightIdx]) {
      if(cmpTime(input[taskOffset+leftIdx], input[rightOffset+rightIdx])) {
        output[writeOffset+writeIdx]=input[taskOffset+leftIdx];
        leftIdx++;
        bitOffset++;
        if(bitOffset == WORDSIZE) {
          bitOffset=0;
          wordOffset++;
        }
        if(leftIdx == step) {
        //printf ("breaking!\n"); 
          break;
        }
      }
      else {
//printf("thread:%d, %d >= %d, right!\n", omp_get_thread_num(), input[taskOffset+leftIdx], input[rightOffset+rightIdx]);
        bitVec[wordOffset]++;
        output[writeOffset+writeIdx]=input[rightOffset+rightIdx];
        rightIdx++;
        bitOffset++;
        if(bitOffset == WORDSIZE) {
          bitOffset=0;
          wordOffset++;
        }
        if(rightIdx == step) {
//          printf("breaking!\n");
          break;
        }
      }
    }
    writeIdx++;
    if(leftIdx < step) {
      while(writeIdx < (taskSize/threadsPerTask)) {
        bitVec[wordOffset] = bitVec[wordOffset] << 1; // shift by 1 place and put a 0 in LSB
        output[writeOffset+writeIdx] = input[taskOffset+leftIdx];
        leftIdx++;
        bitOffset++;
        if(bitOffset == WORDSIZE) {
          bitOffset=0;
          wordOffset++;
        }
        writeIdx++;
      }
    }
    if(rightIdx < step) {
      while(writeIdx < (taskSize/threadsPerTask)) {
        bitVec[wordOffset] = bitVec[wordOffset] << 1; // shift by 1 place and put a 0 in LSB
        bitVec[wordOffset]++;
        output[writeOffset+writeIdx] = input[rightOffset+rightIdx];
        rightIdx++;
        bitOffset++;
        if(bitOffset == WORDSIZE) {
          bitOffset=0;
          wordOffset++;
        }
        writeIdx++;
      }
    }
  }
}
}

// If have enough tasks to sequentially merge...
template<typename T, typename cmp>
void mergeLevelSmall(int N, T* input, T* output, WORDTYPE* bitVec, int step, int tasks, int numThreads) {

//  int leftIdx,rightIdx;
//  int totalIdx=0;
//  int wordOffset=0;
//  int bitOffset=0;
int parThreads =numThreads;
if(tasks < numThreads) parThreads=tasks;
//#pragma omp parallel num_threads(parThreads) shared(input, output, bitVec)
//{
  int threadOffset=0;
//  int threadOffset = omp_get_thread_num() * (N/parThreads);
  int leftIdx,rightIdx;
  int totalIdx=threadOffset;
  int wordOffset=threadOffset/WORDSIZE;
  int bitOffset=0;
  cmp cmpTime;

  // Merge pairs of length step through entire level
//  for(int pairOffset=0; pairOffset < N; pairOffset += 2*step) {
  for(int pairOffset=0; pairOffset < (N/parThreads); pairOffset += 2*step) {
//printf("thread:%d, threadOffset:%d, pairOffset:%d\n", omp_get_thread_num(), threadOffset, pairOffset);
    leftIdx=0;
    rightIdx=step;
    // For each merge pair, continue until either side reaches end
    while(leftIdx < step && rightIdx < (2*step)) {
        bitVec[wordOffset] = bitVec[wordOffset] << 1; // shift by 1 place and put a 0 in LSB
//      if(input[pairOffset+leftIdx+threadOffset] < input[pairOffset+rightIdx+threadOffset]) {
      if(cmpTime(input[pairOffset+leftIdx+threadOffset], input[pairOffset+rightIdx+threadOffset])) {
        output[totalIdx]=input[pairOffset+leftIdx+threadOffset];
        leftIdx++;
      }
      else {
        bitVec[wordOffset]++; // make LSB a 1
        output[totalIdx]=input[pairOffset+rightIdx+threadOffset];
        rightIdx++;
      } 
      totalIdx++;
      bitOffset++;
      if(bitOffset == WORDSIZE) {
        bitOffset=0;
        wordOffset++;
      }
    }
    // If leftover elms from left set
    while(leftIdx < step) {
//printf("%d ", 0);
      bitVec[wordOffset] = bitVec[wordOffset] << 1;
      output[totalIdx]=input[pairOffset+leftIdx+threadOffset];
      leftIdx++;
      totalIdx++;
      bitOffset++;
      if(bitOffset == WORDSIZE) {
        bitOffset=0;
        wordOffset++;
      }
    }
    // If leftover elms from right set
    while(rightIdx < (2*step)) {
//printf("%d ", 1);
      bitVec[wordOffset] = bitVec[wordOffset] << 1;
      bitVec[wordOffset]++;
      output[totalIdx]=input[pairOffset+rightIdx+threadOffset];
      rightIdx++;
      totalIdx++;
      bitOffset++;
      if(bitOffset == WORDSIZE) {
        bitOffset=0;
        wordOffset++;
      }
    }
//    pairOffset = pairOffset + (2*step);
  }
//}
}

// This is where we calc prefix sum offsets for every logn elts.
void calcPrefixBlocks(WORDTYPE* bitVec, int* prefixSum, int N, int numThreads) {
//  int eltsPerThread = (N/WORDSIZE)/numThreads;
  int eltsPerThread = N/WORDSIZE;
//#pragma omp parallel num_threads(numThreads) shared(prefixSum,bitVec)
//  {
    int runningPrefix=0;
//    int threadOffset = eltsPerThread*omp_get_thread_num();
    int threadOffset=0;
    for(int i=0; i<eltsPerThread; i++) {
//printf("thread:%d, i:%d, idx:%d\n", omp_get_thread_num(), i, threadOffset+i);
//      prefixSum[threadOffset+i] = runningPrefix;
//      runningPrefix += __builtin_popcountll(bitVec[threadOffset+i]);
      prefixSum[i] = runningPrefix;
      runningPrefix += __builtin_popcountll(bitVec[i]);
      
    }
//    #pragma omp barrier
/*
    runningPrefix = 0;
    for(int i=0; i<omp_get_thread_num(); i++) {
      runningPrefix += prefixSum[((i+1)*eltsPerThread)-1];
      runningPrefix += __builtin_popcountll(bitVec[((i+1)*eltsPerThread)-1]);
    }
//    #pragma omp barrier
    if(omp_get_thread_num() > 0) {
      for(int i=0; i<eltsPerThread; i++) {
        prefixSum[threadOffset+i] += runningPrefix;
      }
    }
*/
//  }
/*
  int* testPrefixSum = (int*)malloc(N/WORDSIZE*sizeof(int));
  int runningPrefix=0;
  for(int i=0; i<(N/WORDSIZE); i++) {
    testPrefixSum[i] = runningPrefix;
    runningPrefix += __builtin_popcountll(bitVec[i]);
  }
  bool correct=true;
  for(int i=0; i<N/WORDSIZE; i++) {
    if(prefixSum[i] != testPrefixSum[i]) {
      correct=false;
printf("prefix:%d, correct:%d\n", prefixSum[i], testPrefixSum[i]);
    }
  }
  if(!correct)printf("WRONG PREFIX SUM!\n");
*/
}


template<typename T, typename cmp>
//void seqBuildTree(int N, WORDTYPE** bitVec, int bitOffset, T** initTimes, int** prefixSums, int numThreads) {
void buildChazelTree(int N, WORDTYPE** bitVec, int bitOffset, T* initTimes, int** prefixSums, int numThreads) {
  T* tempTimes[2];
  bool tempBit = false;
  tempTimes[0] = initTimes;
  tempTimes[1] = (T*)malloc(N*sizeof(T));

  int step=1;
  int level=0;
  int tasks;

  while(step < N) {
    tasks = N/(step*2);
//    if(tasks >= numThreads) {
//    if(tasks >= 0) {
      mergeLevelSmall<T,cmp>(N, tempTimes[(int)tempBit], tempTimes[(int)!tempBit], &bitVec[level][bitOffset], step, tasks, numThreads);
//    }
//    else {
//      mergeLevelLarge<T,cmp>(N, tempTimes[(int)tempBit], tempTimes[(int)!tempBit], &bitVec[level][bitOffset], step, tasks, numThreads);
//    }

    calcPrefixBlocks(bitVec[level], prefixSums[level], N, numThreads);
  
    step = step*2;
    level++;
    tempBit = !tempBit;
  }
  if(tempBit) {
    for(int i=0; i<N; i++)
      tempTimes[0][i] = tempTimes[1][i];
  }

//  *initTimes = tempTimes[(int)tempBit];
//  return tempTimes[(int)tempBit];
  free(tempTimes[1]);
//  *initTimes = tempTimes[(int)!tempBit];
} 

/*
template<typename T>
void parBuildTree(int N, WORDTYPE** bitVec, int bitOffset, T** initTimes, int** prefixSums) {
  T* tempTimes[2];
  bool tempBit = false;
  tempTimes[0] = *initTimes;
  tempTimes[1] = (int*)malloc(N*sizeof(T));

  int step=1;
  int level=0;

  while(step < N) {
    mergeLevelSmall<T>(N, tempTimes[(int)tempBit], tempTimes[(int)!tempBit], &bitVec[level][bitOffset], step);

    calcPrefixBlocks(bitVec[level], prefixSums[level], N, 1);
  
    step = step*2;
    level++;
    tempBit = !tempBit;
  }

  *initTimes = tempTimes[(int)tempBit];
//  free(tempTimes[(int)!tempBit];
}

*/
template<typename T>
int bSearch(int N, T query, T* list) {
  int idx = N/2;
  int step = N/4;
  while(step >= 1) {
    if(query < list[idx])
      idx -= step;
    else
      idx += step;
    step = step/2;
  }
  if(query < list[idx]) {
    return idx-1;
  }
  return idx;
}

template<typename T, typename cmpType>
int bSearchGen(int N, T query, T* list) {

  int idx = N/2;
  int step = N/4;
  cmpType cmp;
  while(step >= 1) {
//    if(query < list[idx])
    if(cmp(query, list[idx]))
      idx -= step;
    else
      idx += step;
    step = step/2;
  }
//  if(query < list[idx]) {

  if(cmp(query, list[idx])) {
    idx -= 1;
  }
  if(!cmp(query, list[idx])) {
    idx += 1;
  }

  return idx;
}

// Initial naive prefix sum for correctness, calculate word and bit offets internally
int prefixSum(int startIdx, int timeIdx, WORDTYPE* bitVec, int* prefixSums) {
//  timeIdx++;
  int startWord = startIdx/(WORDSIZE);
//  int startBit = startIdx & (WORDSIZE-1);
  int startBit = startIdx %(WORDSIZE);
  int idxWord = timeIdx/(WORDSIZE);
//  int idxBit = timeIdx & (WORDSIZE-1);
  int idxBit = timeIdx %(WORDSIZE);
  WORDTYPE mask;

  int total=0;
//printf("startIdx:%d, startWord:%d, startBit:%d, idx:%d, idxWord:%d, idxBit:%d\n", startIdx, startWord, startBit, timeIdx, idxWord, idxBit);

  // Special case where all contained in 1 word
  if(startWord != idxWord && idxBit == 0) {
    idxWord--;
    idxBit=WORDSIZE-1; 
  }
  if(startWord == idxWord) {
    WORDTYPE startMask;

    if(startBit == 0) startMask = MAX_INT;
    else
      startMask = (1U << (WORDSIZE-(startBit)))-1;

//    int endMask = ((1U << (WORDSIZE-(idxBit+1)))-1) ^ ((1U << WORDSIZE)-1);
    WORDTYPE endMask = ((1U << (WORDSIZE-(idxBit+1)))-1) ^ MAX_INT;
//    total = __builtin_popcountll(bitVec[startWord] & (startMask & endMask));
  }
  else {
// Count relevant bits in start word

    if(startBit == 0) mask = MAX_INT;
    else
      mask = ((1U << (WORDSIZE-(startBit)))-1);

    total += __builtin_popcountll((bitVec[startWord] & mask));

/*
    for(int i=startWord+1; i<idxWord; i++) { // loop through each relevant word
      total += __builtin_popcountll(bitVec[i]);
    }
*/
    total+= prefixSums[idxWord] - prefixSums[startWord+1];

// Count bits at beginning of last word
//    mask = ((1U << (WORDSIZE-(idxBit+1)))-1) ^ ((1U << WORDSIZE)-1);
    mask = ((1U << (WORDSIZE-(idxBit+1)))-1) ^ MAX_INT;
    total += __builtin_popcountll((bitVec[idxWord] & mask));
  }
//printf("startWord:%d, startBit:%d, timeIdx:%d, idxWord:%d, idxBit:%d, total:%d\n", startWord, startBit, timeIdx, idxWord, idxBit,total);
/*
WORDTYPE temp=bitVec[startWord];
int brute=0;
for(int i=0; i<WORDSIZE-startBit; i++) {
  brute += (temp & 1);
  temp = temp >> 1;
}
brute += prefixSums[idxWord] - prefixSums[startWord+1];
temp = bitVec[idxWord];
  temp = temp >> idxBit;
for(int i=idxBit; i<WORDSIZE; i++) {
  brute += (temp & 1);
  temp = temp >> 1;
} 
  printf("%d %d\n", total, brute);
*/
  return total;
}

template<typename T, typename V>
int seqQuery(int N, int iLevels, T qTime, V qVal, V* vals, T* sortedTimes, WORDTYPE** bitVec, int** prefixSums) {
  int level=iLevels-1;
  int timeIdx = bSearch<T>(N, qTime, sortedTimes);
  int rightTimeElts = prefixSum(0, timeIdx, bitVec[iLevels-1], prefixSums[iLevels-1]);
  int leftTimeElts = (timeIdx+1) - rightTimeElts;
  int totalCount=0;

//printf("initial timeIdx:%d, rightElts:%d, leftElts:%d\n", timeIdx, rightTimeElts, leftTimeElts);
  int valIdx = N/2;
  int step = N/4;
  int sumStart=0;


// perform binary search while keeping track of prefix sums of bits
  while(step >= 1) {
    level--;
    if(qVal >= vals[valIdx]) {
      totalCount += leftTimeElts;
      timeIdx = rightTimeElts-1;
      if(timeIdx<0) return totalCount;
      sumStart += (step*2);
      rightTimeElts = prefixSum(sumStart, sumStart+timeIdx, bitVec[level], prefixSums[level]);
      leftTimeElts = (timeIdx+1) - rightTimeElts;
      valIdx += step;
    } 
    else {
      timeIdx = leftTimeElts-1;
      if(timeIdx<0) return totalCount;
      valIdx -= step;
      rightTimeElts = prefixSum(sumStart, sumStart+timeIdx, bitVec[level], prefixSums[level]);
      leftTimeElts = (timeIdx+1) - rightTimeElts;
    }
    step = step/2;
  }
  if(leftTimeElts>0) {
    if(qVal >= vals[valIdx-1])
      totalCount++;
  }
  if(rightTimeElts>0) {
    if(qVal >= vals[valIdx])
      totalCount++;
  }

//  totalCount += leftTimeElts;
  return totalCount;
   
}

template<typename T, typename V, typename timeCmp, typename valCmp>
int seqQueryAbove(int N, int iLevels, T qTime, V qVal, V* vals, T* sortedTimes, WORDTYPE** bitVec, int** prefixSums) {
  int level = iLevels-1;
  int timeIdx = 0;
  timeIdx = bSearchGen<T, timeCmp>(N, qTime, sortedTimes);
if(timeIdx > N) printf("D - %d\n", timeIdx);
  int rightTimeElts = prefixSum(0, timeIdx, bitVec[iLevels-1], prefixSums[iLevels-1]);
  int leftTimeElts = (timeIdx+1) - rightTimeElts;
  int totalCount=0;

  int valIdx = N/2;
  int step = N/4;
  int sumStart=0;

  valCmp cmp;

// perform binary search while keeping track of prefix sums of bits
  while(step >= 1) {
    level--;
    if(cmp(qVal, vals[valIdx])) {
      totalCount += leftTimeElts;
      timeIdx = rightTimeElts-1;
      if(timeIdx<0) return totalCount;
      sumStart += (step*2);
      if((sumStart+timeIdx) >=N) printf("N:%d, E - %d, %d\n", N, sumStart, timeIdx);
      rightTimeElts = prefixSum(sumStart, sumStart+timeIdx, bitVec[level], prefixSums[level]);
      leftTimeElts = (timeIdx+1) - rightTimeElts;
      valIdx += step;
    } 
    else {
      timeIdx = leftTimeElts-1;
      if(timeIdx<0) return totalCount;
      valIdx -= step;
      if((sumStart+timeIdx) >=N) printf("N:%d, F - %d, %d!\n", N, sumStart, timeIdx);
      rightTimeElts = prefixSum(sumStart, sumStart+timeIdx, bitVec[level], prefixSums[level]);
      leftTimeElts = (timeIdx+1) - rightTimeElts;
    }
    step = step/2;
  }
  if(leftTimeElts>0) {
    if(cmp(qVal, vals[valIdx-1])) {
      totalCount++;
    }
  }
  if(rightTimeElts>0) {
    if(cmp(qVal, vals[valIdx])) {
      totalCount++;
    }
  }
//  totalCount += leftTimeElts;
  return totalCount;
}

template<typename T, typename V, typename timeCmp, typename valCmp>
int seqQueryUnder(int N, int iLevels, T qTime, T qVal, V* vals, T* sortedTimes, WORDTYPE** bitVec, int** prefixSums, int* relResult, float* pi) {
  int level=iLevels-1;
  int timeIdx=0;
  timeIdx = bSearchGen<T, timeCmp>(N, qTime, sortedTimes);
if(timeIdx > N) printf("A - %d\n", timeIdx);
  *relResult = timeIdx;
  int rightTimeElts = prefixSum(0, timeIdx-1, bitVec[iLevels-1], prefixSums[iLevels-1]);
  int leftTimeElts = (timeIdx) - rightTimeElts;
  int totalCount=0;

//printf("initial timeIdx:%d, rightElts:%d, leftElts:%d\n", timeIdx, rightTimeElts, leftTimeElts);
  int valIdx = N/2;
  int step = N/4;
  int sumStart=0;

  valCmp cmp;

// perform binary search while keeping track of prefix sums of bits
  while(step >= 1) {
    level--;
//    if(qVal >= vals[valIdx]) {
//printf("valSlope:%f, valIcp:%f, cmp:%d\n", vals[valIdx].slope, vals[valIdx].icp, cmp(qTime, qVal, vals[valIdx]));
    if(cmp(qTime, qVal, vals[valIdx])) {
//printf("above %d! qT:%.2f, qV:%.2f, vS:%.2f, vI:%.2f\n", valIdx, qTime, qVal, vals[valIdx].slope, vals[valIdx].icp);
      totalCount += leftTimeElts;
      timeIdx = rightTimeElts-1;
      if(timeIdx<0) return totalCount;
      sumStart += (step*2);
      rightTimeElts = prefixSum(sumStart, sumStart+timeIdx, bitVec[level], prefixSums[level]);
      leftTimeElts = (timeIdx+1) - rightTimeElts;
      valIdx += step;
    } 
    else {
      timeIdx = leftTimeElts-1;
      if(timeIdx<0) return totalCount;
      valIdx -= step;
      rightTimeElts = prefixSum(sumStart, sumStart+timeIdx, bitVec[level], prefixSums[level]);
      leftTimeElts = (timeIdx+1) - rightTimeElts;
    }
    step = step/2;
  }
  *pi = vals[valIdx].slope;

/*
  if(leftTimeElts>0) {
//    if(qVal >= vals[valIdx-1])
    if(cmp(qTime, qVal, vals[valIdx-1])) {
      totalCount++;
    }
  }
  if(rightTimeElts>0) {
//    if(qVal >= vals[valIdx])
    if(cmp(qTime, qVal, vals[valIdx])) {
      totalCount++;
    }
  }
*/
//  totalCount += leftTimeElts;
  return totalCount;
   
}
