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
#include<stdio.h>

//=====Constructors=============================================================

/*
DualRay::DualRay(CriticalRay* r, SliceCell* c) {
    leftOriented = r->isLeftOriented();
    cellPTR = c;
    a = r->getStartingPoint()->getX();
    b = -r->getStartingPoint()->getY();
    //even if there is no starting point there exist a vector
    startX = r->getStartingPoint()->getX();
    startY = r->getStartingPoint()->getY();
    otherX = r->getOtherPoint()->getX();
    otherY = r->getOtherPoint()->getY();
    num = r->getOtherPoint()->getY() - r->getStartingPoint()->getY();
    denom = r->getOtherPoint()->getX() - r->getStartingPoint()->getX();
    endX = r->getSlope();
    endY = -r->getIntercept();
    immBelow = NULL;
}
*/

/*
DualRay::DualRay(int x1, int y1, int x2, int y2, int a, int b, bool isL)
: startX(x1), startY(y1), otherX(x2), otherY(y2), a(a), b(b), leftOriented(isL) {
    num = y2 - y1;
    denom = x2 - x1;
    immBelow = NULL;
}

DualRay::DualRay(const DualRay& orig) {
    startX = orig.startX;
    startY = orig.startY;
    otherX = orig.otherX;
    otherY = orig.otherY;
    num = orig.num;
    denom = orig.denom;
    a = orig.a;
    b = orig.b;
    leftOriented = orig.leftOriented;
}
*/

DualRay::DualRay(float setA) {
    a = setA;
    immBelow = NULL;
}

DualRay& DualRay::operator=(const DualRay& orig) {
//    startX = orig.startX;
//    startY = orig.startY;
//    otherX = orig.otherX;
//    otherY = orig.otherY;
//    num = orig.num;
//    denom = orig.denom;
    a = orig.a;
    b = orig.b;
    leftOriented = orig.leftOriented;
    label = orig.label;
    return *this;
}


//bool DualRay::operator<(const DualRay& a) {
//  printf("%d vs %d, returning %d\n", this.startX, a.startX, this.startX < a.startX);
///  return this.startX < a.startX;
//}


DualRay::~DualRay() {
    //do not delete cell itself, only the pointer needs to be deleted
}

//=====Getters and Setters======================================================

void DualRay::setA(float setA) {
  a = setA;
  immBelow=NULL;
}

float DualRay::getA() {
    return a;
}

float DualRay::getB() {
    return b;
}

bool DualRay::isLeftOriented() {
    return leftOriented;
}
/*

int DualRay::getStartX() {
    return startX;
}

int DualRay::getStartY() {
    return startY;
}

int DualRay::getOtherX() {
    return otherX;
}

int DualRay::getOtherY() {
    return otherY;
}

int DualRay::getNum() {
    return num;
}

int DualRay::getDenom() {
    return denom;
}
*/

float DualRay::getEndX() {
    return endX;
}

float DualRay::getEndY() {
    return endY;
}

//=====Methods==================================================================

/*
bool DualRay::startsBeforeOther(DualRay* ray) {
    //we know that left oriented rays start at -infinity and therefore are always before the other

    if (leftOriented) {
        //this is before the other or equal to the other
        //either because this has no starting point or because it is left oriented
        return true;
    } else if (ray->isLeftOriented()) { //(&& !leftOriented)
        //the other is by definition before me
        //this has a starting point and is right oriented, while the other is left oriented
        return false;
    } else {
        //both are right oriented and have a starting point so this means that the start points determine it all
        return (num * ray->getDenom()) < (ray->getNum() * denom);
    }
}

bool DualRay::startingPointBeforeOther(DualRay* ray) {      
    if ((leftOriented && !ray->isLeftOriented()) || (!leftOriented && ray->isLeftOriented())) {
        //this is left oriented, sign flips (denom is < 0) or the other is left oriented, sign flips (ray.denom < 0)
        return (num * ray->getDenom()) > (ray->getNum() * denom);
    } else { 
        //none of them is left oriented, or both are left oriented
        return (num * ray->getDenom()) < (ray->getNum() * denom);
    }
}
*/

bool DualRay::rAbove(DualRay* b) {
    double aY = endY;
    double bY = (((double)b->getA())*((double)endX) + (double)b->getB());
    double dist = aY - bY;
    bool retVal;
    if(dist==0) {
      retVal = (a < b->getA());
    }
    retVal= (dist < 0.0);
//printf("a:r%d, slope:%d, intercept:%d, endX:%lf... b:r%d, slope:%d, intercept:%d, endX:%lf - dist:%lf, retVal:%d\n", label, a, getB(), endX, b->label, b->getA(), b->getB(), b->getEndX(), dist, retVal);
    return retVal;
}
