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


#include "DualRay.h"
#include <cstddef>
#include <functional>

class PNode {

  public:
    PNode* parent;
    PNode* leftChild;
    PNode* rightChild;
    DualRay* theRay;
int depth;
    int subTreeSize;
    int weight;
    int balance;

    PNode(DualRay* newRay) {
      theRay = newRay;
      leftChild = NULL;
      rightChild = NULL;
      parent= NULL;
      depth = 0;
      subTreeSize=1;
      weight=1;
      balance=0;
    }
    void initialize(DualRay* newRay) {
      theRay = newRay;
      leftChild = NULL;
      rightChild = NULL;
      parent = NULL;
      depth = 0;
      subTreeSize=1;
      weight=1;
      balance=0;
    }

};

void printTree(PNode* root);

void deletePTree(PNode* root);

PNode* searchPTree(DualRay* target, PNode* root, bool update);

int underCount(DualRay* target, PNode* root);

int aboveCount(DualRay* target, PNode* root, bool (*f)(DualRay*,DualRay*));

void insertPNode(DualRay* newRay, PNode* root, PNode* newNode, bool (*f)(DualRay*,DualRay*));

void insertPNodeDups(DualRay* newRay, PNode* root, PNode* newNode, float (*f)(DualRay*,DualRay*));
