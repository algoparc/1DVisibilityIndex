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
#include <algorithm>
#include <unistd.h>
#include <chrono>
//#include "convexhull.h"
#include "BuildConvexHull.h"
using namespace std;

/*    
void SortPoint(Point* p, int n)
{
    Point *min,temp;
    int i,j;
    min=&p[0];
    for(i=1; i<n; i++)
        if(p[i].y<min->y)
            min=&p[i];
    temp=*min;
    *min=p[0];
    p[0]=temp;
    for(i=2; i<n; i++)
    {
        j=i-1;
        temp=p[i];
        while(j>0&&((temp.x-p[0].x)*(p[j].y-p[0].y)-(p[j].x-p[0].x)*(temp.y-p[0].y))>0)
        {
            p[j+1]=p[j];
            j--;
        }
        p[j+1]=temp;
    }
    p[n]=p[0];
}
*/
/*
int CreateConvexHull(int size, Point* p) {
    int j,Num,k,t;
    Num=size;
    for(k=3; k<size || k<=Num; k++)
    {
        j=2;
        for(j=2; j<=k-1; j++)
        {
            if(((p[k-j+1].x-p[k-j].x)*(p[k-j].y-p[0].y)-(p[k-j].x-p[0].x)*(p[k-j+1].y-p[k-j].y))*((p[k-j+1].x-p[k-j].x)*(p[k].y-p[k-j+1].y)-(p[k].x-p[k-j+1].x)*(p[k-j+1].y-p[k-j].y))>0)
            {
                for(t=k-1; t<Num; t++)
                    p[t]=p[t+1];
                Num--;
                k=k-1;
                j=j-1;
            }
        }
    }
    return Num;
}

*/
void findTangentToHull(Point* points, int hullSize, Ray* searchRay, bool up) {
  Point tangetPoint;
  float tangetSlope = calcSlope(searchRay->x, searchRay->y, points[0].x, points[0].y);
  tangetPoint = points[0];
  for(int i=0; i<hullSize; i++) {
    if(up) {
      if(tangetSlope < calcSlope(searchRay->x, searchRay->y, points[i].x, points[i].y)) {
        tangetSlope = calcSlope(searchRay->x, searchRay->y, points[i].x, points[i].y);
        tangetPoint = points[i];
      }
    }
    else {
      if(tangetSlope > calcSlope(searchRay->x, searchRay->y, points[i].x, points[i].y)) {
        tangetSlope = calcSlope(searchRay->x, searchRay->y, points[i].x, points[i].y);
        tangetPoint = points[i];
      }
    }
  }
  searchRay->slope = tangetSlope;
  searchRay->icp = calcIcp(searchRay->x, searchRay->y, tangetSlope);
}

void findTangentToHullPar(Hull* hull, Ray* searchRay, bool up) {
//  float tangetSlope = calcSlope(searchRay->x, searchRay->y, hull->crits[0]->x, hull->crits[0]->y);
  float tangetSlope = searchRay->slope;
  for(int i=0; i<hull->size; i++) {
    if(up) {
      if(tangetSlope < calcSlope(searchRay->x, searchRay->y, hull->crits[i]->x, hull->crits[i]->y)) {
        tangetSlope = calcSlope(searchRay->x, searchRay->y, hull->crits[i]->x, hull->crits[i]->y);
//        tangetRay = hull->crits[i];
      }
    }
    else {
      if(tangetSlope > calcSlope(searchRay->x, searchRay->y, hull->crits[i]->x, hull->crits[i]->y)) {
        tangetSlope = calcSlope(searchRay->x, searchRay->y, hull->crits[i]->x, hull->crits[i]->y);
//        tangetRay = hull->crits[i];
      }
    }
  }
  searchRay->slope = tangetSlope;
  searchRay->icp = calcIcp(searchRay->x, searchRay->y, tangetSlope);
}



void updateCriticalRays(Ray* Rcrit, Ray* Lcrit, int step, int threads) {

  int hullSize=0;
  Point* points = (Point*)malloc(step*sizeof(Point));
// Update Lcrit

  for(int i=0; i<step; i++) {
    points[i].x = Rcrit[i].x;
    points[i].y = Rcrit[i].y;
  }
 

  //SortPoint(points, step);
//  hullSize = CreateConvexHull(step, points);
  hullSize = buildConvexHull(points, step);

#pragma omp parallel for num_threads(threads)
  for(int i=0; i<step; i++) {
    findTangentToHull(points, hullSize, &Lcrit[i], 1);
  }

// Update Rcrit
  for(int i=0; i<step; i++) {
    points[i].x = Lcrit[i].x;
    points[i].y = Lcrit[i].y;
  }
//  SortPoint(points, step);
//  hullSize = CreateConvexHull(step, points);
  hullSize = buildConvexHull(points, step);

#pragma omp parallel for num_threads(threads)
  for(int i=0; i<step; i++) {
    findTangentToHull(points, hullSize, &Rcrit[i], 0);
  }
  free(points);
}

Hull* mergeHulls(Hull* lHull, Hull* rHull, int threads) {
  // for left side, binary search while checking slope
  // temp just scan to find point
  // Count number of points to copy to new hull
//  int lCount=1;
  int lCount=0;
#pragma omp parallel for num_threads(threads) shared(lCount)
  for(int i=1; i<lHull->size-1; i++) {
    if((lHull->crits[i-1]->slope <= lHull->crits[i]->slope) && (lHull->crits[i]->slope > lHull->crits[i+1]->slope)) {
      lCount=i+1;
    }
  }

  // for right side, binary search while checking slope
  // temp just scan
  int rCount=0;
#pragma omp parallel for num_threads(threads) shared(lCount)
  for(int i=1; i<rHull->size-1; i++) {
    if((rHull->crits[i-1]->slope >= rHull->crits[i]->slope) && (rHull->crits[i]->slope < rHull->crits[i+1]->slope)) {
      rCount=rHull->size-i;
    }
  }

//  for(int i=rHull->size-1; i>0 && (rHull->crits[i]->slope <= rHull->crits[i-1]->slope); i--) {
//    rCount++;
//  }

  // copy points from left and right into new hull
  Hull* newHull = (Hull*)malloc(sizeof(Hull));
  newHull->crits = (Ray**)malloc((lCount+rCount)*sizeof(Ray*));

  // Copy previous points from lHull and rHull

#pragma omp parallel for num_threads(threads)
  for(int i=0; i<lCount; i++) {
    newHull->crits[i] = lHull->crits[i];
  }
#pragma omp parallel for num_threads(threads)
  for(int i=0; i<rCount; i++) {
    newHull->crits[i+lCount] = rHull->crits[(rHull->size-rCount)+i];
  }
  newHull->size = lCount+rCount;
  free(lHull->crits);
  free(rHull->crits);
  return newHull;
}


Hull* updateCriticalRaysReturnHull(Ray* Rcrit, Ray* Lcrit, int step, int threads) {


// Update Lcrit
  Hull* lHull = buildAndSaveHull(Rcrit,step);

#pragma omp parallel for num_threads(threads)
  for(int i=0; i<step; i++) {
    findTangentToHullPar(lHull, &Lcrit[i], 0);
  }

// Update Rcrit
  Hull* rHull = buildAndSaveHull(Lcrit, step);

#pragma omp parallel for num_threads(threads)
  for(int i=0; i<step; i++) {
    findTangentToHullPar(rHull, &Rcrit[i], 1);
  }

  Hull* newHull = mergeHulls(lHull, rHull, threads);
  return newHull;
}

Hull* updateCriticalRaysPar(Ray* Rcrit, Ray* Lcrit, Hull* lHull, Hull* rHull, int step, int threads)
{
#pragma omp parallel for num_threads(threads)
  for(int i=0; i<step; i++) {
    findTangentToHullPar(rHull, &Rcrit[i], 1);
  }
#pragma omp parallel for num_threads(threads)
  for(int i=0; i<step; i++) {
    findTangentToHullPar(lHull, &Lcrit[i], 0);
  }

  Hull* newHull = mergeHulls(lHull, rHull, threads);
  return newHull;
}
