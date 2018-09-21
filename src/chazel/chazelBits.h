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
#include<math.h>
#include "params.h"


// If blocks are small, only store bitVec values
template<typename T>
void mergeLevelSmall(int N, T* input, T* output, bool* bitVec, int step);

// If blocks are large, need to also store prefixBlocks array values
void mergeLevelLarge(int N, int* input, int* output, bool* bitVec, int step, int* prefixBlocks, int blockSize);

// This is where we calc prefix sum offsets for every logn elts.
void calcPrefixBlocks(bool* bitVec, int* prefixSum); 

template<typename T>
void seqBuildTree(int N, WORDTYPE** bitVec, int offset, T** initTimes, int** prefixSums);

int bSearch(int N, int query, int* list); 

// Initial naive prefix sum for correctness
int prefixSum(int timeIdx, bool* bitVec);

int seqQuery(int N, int iLevels, int qTime, int qVal, int* vals, int* sortedTimes, bool** bitVec, int** prefixSums); 
