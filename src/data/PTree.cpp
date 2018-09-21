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


#include "PTree.h"
#include <stdio.h>
#include <queue>

void printTree(PNode *root) {
   std::queue<PNode*> theQueue;
   theQueue.push(root);
   int currDepth=0;
   PNode* nodePtr=NULL;

   while(theQueue.size() > 0) {
     nodePtr = theQueue.front();
     theQueue.pop();
     if(currDepth != nodePtr->depth) {
       printf("\n");
       currDepth = nodePtr->depth;
     }
//     printf("%d(%d)(%d) ", nodePtr->theRay->label, nodePtr->subTreeSize, nodePtr->weight);
    // printf(" %04x(%d) ", nodePtr, nodePtr->balance);
     if(nodePtr->leftChild != NULL) {
       printf("L");
       theQueue.push(nodePtr->leftChild);
     }
     if(nodePtr->rightChild != NULL) {
       printf("R");
       theQueue.push(nodePtr->rightChild);
     }
  }
printf("\n");
}

double compareRays(DualRay* a, DualRay* b) {
  double aY = (double)a->getEndY();
  double bY = (((double)b->getA())*((double)a->getEndX())) + (double)b->getB();
  double dist = aY - bY;
//printf("aY:%lf, bY:%lf\n", aY, bY);
  if(dist==0.0) {
    dist = dist + (a->getA() < b->getA()) - (a->getA() > b->getA());
  }
  return dist;
  
}

double comparePointAndRay(DualRay* a, DualRay* b) {
  double aY = (double)a->getEndY();
  double bY = (((double)b->getA())*((double)a->getEndX())) + (double)b->getB();
  double dist = aY - bY;
  if(dist==0.0) {
    dist = -1.0;
  }

  return dist;
}

void deletePTree(PNode *root){
    if(root!=NULL)
    {
        if(root->leftChild != NULL)
          deletePTree(root->leftChild);
        if(root->rightChild != NULL)
          deletePTree(root->rightChild);

        if(root->leftChild!=NULL)
            root->leftChild=NULL;
        if(root->rightChild!=NULL)
            root->rightChild=NULL;
        delete root;
        root=NULL;
    }
}

PNode* searchPTree(DualRay* target, PNode* root, bool updateSize) {
  PNode* nodePtr = root;
  double dist;
  while (nodePtr != NULL) {
    if(updateSize)
      nodePtr->subTreeSize = nodePtr->subTreeSize+1;
    dist = compareRays(target, nodePtr->theRay);
    if(dist == 0)
      break;
    else if(dist < 0) {
      if(nodePtr->leftChild == NULL)
        break;
      else
        nodePtr = nodePtr->leftChild;
    }
    else {
      if(nodePtr->rightChild == NULL)
        break;
      else
        nodePtr = nodePtr->rightChild;
    }
  }
  return nodePtr;
}

int underCount(DualRay* target, PNode* root) {
  int count = 0;
  double dist=0.0;
  PNode* nodePtr = root;
  while(nodePtr != NULL) {
    dist = comparePointAndRay(target, nodePtr->theRay);
    if(dist == 0){
      if(nodePtr->leftChild != NULL) {
        count = count + (nodePtr->leftChild->subTreeSize);
      }
      break;
    }
    else if(dist < 0) {
      if(nodePtr->leftChild == NULL) 
        break;
      else {
        nodePtr = nodePtr->leftChild;
      }
    }
    else {
      target->immBelow = nodePtr->theRay;
      if(nodePtr->leftChild != NULL) {
        count += nodePtr->leftChild->subTreeSize;
}
      if(nodePtr->rightChild == NULL) {
        count++;
        break;
      }
      else {
        nodePtr = nodePtr->rightChild;
        count++;
      }
    }
  }
  return count;
}

int aboveCount(DualRay* target, PNode* root, bool (*f)(DualRay*,DualRay*)) {
bool debug=false;
  int count=0; 
  PNode* nodePtr = root;
  bool cmpResult;
  while(nodePtr != NULL) {
    cmpResult = (*f)(target,nodePtr->theRay);
if(debug)
printf("Compare with:r%d\n", nodePtr->theRay->label);
    if(!cmpResult) {
      if(nodePtr->leftChild == NULL) 
        break;
      else {
        nodePtr = nodePtr->leftChild;
      }
    }
    else {
      if(nodePtr->leftChild != NULL) {
        count += nodePtr->leftChild->subTreeSize;
if(debug)
printf("1 incrementing by %d, now %d\n", nodePtr->leftChild->subTreeSize, count);
      }
      if(nodePtr->rightChild == NULL) {
        count+= nodePtr->weight;
if(debug)
printf("2 incrementing by %d, now %d\n", nodePtr->weight, count);
        break;
      }
      else {
        count+= nodePtr->weight;
        nodePtr = nodePtr->rightChild;
if(debug)
printf("3 incrementing by %d, now %d\n", nodePtr->weight, count);
      }
    }
  }
  return count; 
}
/*
// Insert rays into BST by their slope.
void insertPNode(DualRay* newRay, PNode* root, PNode* newNode, bool (*f)(DualRay*,DualRay*)){
//  PNode* newNode = new PNode(newRay);
  newNode->initialize(newRay);
  PNode* nodePtr = root;
  int depth=1;
  while (nodePtr != NULL) {
    nodePtr->subTreeSize = nodePtr->subTreeSize+1;
//    if(newRay->getA() < nodePtr->theRay->getA()){
    if((*f)(newRay, nodePtr->theRay)) {
      if(nodePtr->leftChild == NULL) {
        newNode->depth=depth;
        nodePtr->leftChild = newNode;
        newNode->parent = nodePtr;
        nodePtr=NULL;
      }
      else
        nodePtr = nodePtr->leftChild;
    }
    else {
      if(nodePtr->rightChild == NULL) {
        newNode->depth=depth;
        nodePtr->rightChild = newNode;
        newNode->parent = nodePtr;
        nodePtr=NULL;
      }
      else
        nodePtr = nodePtr->rightChild;
    }
    depth++;
  }
*/
void increaseDepth(PNode* node) {
  if(node != NULL) {
    node->depth++;
    increaseDepth(node->leftChild);
    increaseDepth(node->rightChild);
  }
}

void decreaseDepth(PNode* node) {
  if(node != NULL) {
    node->depth--;
    decreaseDepth(node->leftChild);
    decreaseDepth(node->rightChild);
  }
}
    


void rotate_left(PNode* node) {
//printf("rotate_left node:%04x\n", node);
//  node = node->parent;
//printf("node:%04x, rightChild:%04x\n", node, node->rightChild);
//  PNode* parent = node->parent;
  PNode* pivot = node->rightChild;
//if(node->parent != NULL) {
//  node = pivot->parent;
if(node->parent != NULL) {
  if(node->parent->leftChild == node)
    node->parent->leftChild = pivot;
  else
    node->parent->rightChild = pivot;
}

increaseDepth(node->leftChild);
decreaseDepth(pivot->rightChild);
  
//  if(child != NULL) {
// Fix parent node links
  if(pivot->leftChild != NULL) {
    node->rightChild->parent = node;
 }
  node->rightChild = pivot->leftChild;
  pivot->leftChild = node;
  pivot->parent = node->parent;
  node->parent = pivot;
  //printf("increment %04x and %04x", node, pivot);
  pivot->balance += 1;
  node->balance += 1;
  node->depth++;
  pivot->depth--;
//  } 
// Fix subtree counts
//  int temp = node->subTreeSize;
//  node->subTreeSize -= child->subTreeSize;
//  if(child->leftChild != NULL) {
//    printf("leftChild not null\n");
//    node->subTreeSize += child->leftChild->subTreeSize;
//    child->leftChild->parent = node;
//    node->rightChild = child->leftChild;
//  }
//  else {
//    node->rightChild = NULL;
//   printf("Setting right child to NULL!\n");
  
//  child->subTreeSize = temp;

// Set new pointers
//  child->leftChild = node;
}

void rotate_right(PNode* node) {
//printf("rotate_right node:%04x\n", node);
  PNode* pivot = node->leftChild;
//if(pivot->parent != NULL) {
//  node = pivot->parent;
  
//  if(child != NULL) {
// Fix parent node links

if(node->parent != NULL) {
  if(node->parent->leftChild == node)
    node->parent->leftChild = pivot;
  else
    node->parent->rightChild = pivot;
}

increaseDepth(node->rightChild);
decreaseDepth(pivot->leftChild);

  if(pivot->rightChild != NULL) {
    node->leftChild->parent = node;
  }
  node->leftChild = pivot->rightChild;
  pivot->rightChild = node;
  pivot->parent = node->parent;
  node->parent = pivot;
  //printf("decrement %04x and %04x", node, pivot);
  pivot->balance -= 1;
  node->balance -= 1;
  node->depth++;
  pivot->depth--;
//}
//  node = node->parent;
  PNode* child = node->leftChild;
//  if(child != NULL) {
// Fix parent node links
  child->parent = node->parent;
  node->parent = child;

// Fix subtree counts
  int temp = node->subTreeSize;
  node->subTreeSize -= child->subTreeSize;
  if(child->rightChild != NULL) {
    node->subTreeSize += child->rightChild->subTreeSize; 
    child->rightChild->parent = node;
    node->leftChild = child->rightChild;
  }
  else
    node->leftChild = NULL;
  child->subTreeSize = temp;
}
// Set new pointers
//  child->rightChild = node;
//printf("Before exiting:%04x\n", node->parent);
//}


void balancePTree(PNode* newLeaf) {
  PNode* N = newLeaf;
  PNode* P = newLeaf->parent;
// N is the child of P whose height increases by 1.
 while(P->parent != NULL) {
   // balance_factor(P) has not yet been updated!
   if (N == P->leftChild) { // the left subtree increases
     if (P->balance == 1) { // The left column in the picture
       // ==> the temporary balance_factor(P) == 2 ==> rebalancing is required.
       if ((N->balance) == -1) { // Left Right Case
          rotate_left(N); // Reduce to Left Left Case
       }
       // Left Left Case
       rotate_right(P);
       break; // Leave the loop
     }
     if (P->balance == -1) {
       P->balance = 0; // N’s height increase is absorbed at P.
       break; // Leave the loop
     }
     P->balance = 1; // Height increases at P
   } else { // N == right_child(P), the child whose height increases by 1.
     if (P->balance == -1) { // The right column in the picture
       // ==> the temporary balance_factor(P) == -2 ==> rebalancing is required.
       if (N->balance == 1) { // Right Left Case
          rotate_right(N); // Reduce to Right Right Case
       }
       // Right Right Case
       rotate_left(P);
       break; // Leave the loop
     }
     if (P->balance == 1) {
       P->balance = 0; // N’s height increase is absorbed at P.
       break; // Leave the loop
     }
     P->balance = -1; // Height increases at P
   }
   N = P;
   P = P->parent;
  }
// } while (P != NULL && P->parent != NULL); // Possibly up to the root
}


// Insert rays into BST by their slope.
void insertPNode(DualRay* newRay, PNode* root, PNode* newNode, bool (*f)(DualRay*,DualRay*)){
//  PNode* newNode = new PNode(newRay);
//printf("Insert.\n");
  newNode->initialize(newRay);
  PNode* nodePtr = root;
  int depth=1;
  while (nodePtr != NULL) {
    nodePtr->subTreeSize = nodePtr->subTreeSize+1;
//    if(newRay->getA() < nodePtr->theRay->getA()){
    if((*f)(newRay, nodePtr->theRay)) {
      if(nodePtr->leftChild == NULL) {
        newNode->depth=depth;
        nodePtr->leftChild = newNode;
        newNode->parent = nodePtr;
          //balancePTree(newNode);

        nodePtr=NULL;
      }
      else
        nodePtr = nodePtr->leftChild;
    }
    else {
      if(nodePtr->rightChild == NULL) {
        newNode->depth=depth;
        nodePtr->rightChild = newNode;
        newNode->parent = nodePtr;
          //balancePTree(newNode);
        nodePtr=NULL;
      }
      else
        nodePtr = nodePtr->rightChild;
    }
    depth++;
  }
}

void insertPNodeDups(DualRay* newRay, PNode* root, PNode* newNode, float (*f)(DualRay*,DualRay*)){
//  PNode* newNode = new PNode(newRay);
  newNode->initialize(newRay);
  PNode* nodePtr = root;
  int depth=1;
  float dist=0.0;
  while (nodePtr != NULL) {
    nodePtr->subTreeSize = nodePtr->subTreeSize+1;
//    if(newRay->getA() < nodePtr->theRay->getA()){
    dist = (*f)(newRay, nodePtr->theRay);
    if(dist == 0.0) { // Condense into 1 node
      nodePtr->weight++;
      nodePtr=NULL;
    }
    else if(dist > 0.0) {
      if(nodePtr->leftChild == NULL) {
        newNode->depth=depth;
        nodePtr->leftChild = newNode;
        nodePtr=NULL;
      }
      else
        nodePtr = nodePtr->leftChild;
    }
    else {
      if(nodePtr->rightChild == NULL) {
        newNode->depth=depth;
        nodePtr->rightChild = newNode;
        nodePtr=NULL;
      }
      else
        nodePtr = nodePtr->rightChild;
    }
    depth++;
  }
}
