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


#include "BuildConvexHull.h"
#include<stdlib.h>
#include<stdio.h>

//=====Constructors=============================================================

/* Not needed to implement as there are no variables etc.
BuildConvexHull::BuildConvexHull() {
}

BuildConvexHull::BuildConvexHull(const BuildConvexHull& orig) {
}

BuildConvexHull& BuildConvexHull::operator=(const BuildConvexHull& orig) {
}

BuildConvexHull::~BuildConvexHull() {
}*/

//=====Getters and Setters======================================================


//=====Methods==================================================================

Hull* buildAndSaveHull(Ray* crits, int numPoints) {
    //The temporary list that contains the convex hull points (including the artificial points)
    Hull* tmpHull = (Hull*)malloc(sizeof(Hull));
    tmpHull->crits = (Ray**)malloc((numPoints)*sizeof(Ray*));
    tmpHull->size = 0; //Denotes the current size of the convex hull

    tmpHull->crits[0] = &crits[0];
    tmpHull->size += 1;

    for (int i = 0; i < numPoints; i++) {
        //remove all points that violate the convex hull property
        while (tmpHull->size >= 2 && //and the sequence of the last two points with the selected point does not make a counter-clockwise turn
               (1LL*(tmpHull->crits[tmpHull->size - 1]->y - tmpHull->crits[tmpHull->size - 2]->y) * (crits[i].x - tmpHull->crits[tmpHull->size - 1]->x)
                - 1LL*(tmpHull->crits[tmpHull->size - 1]->x - tmpHull->crits[tmpHull->size - 2]->x) * (crits[i].y - tmpHull->crits[tmpHull->size - 1]->y))
               <= 0LL) {
            //Remove the last point from tmpHull
//            tmpHull[tmpHull->size] = nullptr;
            tmpHull->size -= 1;
        }
        //append the last selected point to
        tmpHull->crits[tmpHull->size] = &crits[i];
        tmpHull->size += 1;
    }
    //append the last point
//    tmpHull->points[tmpHull->size] = *e;
//    tmpHull->crits[tmpHull->size] = &crits[numPoints - 1];
//    tmpHull->size += 1;
    //optimise array size and require dynamically allocated array, because of sharing
    Hull* newHull = (Hull*)malloc(sizeof(Hull));
    newHull->size = tmpHull->size;
    newHull->crits = (Ray**)malloc(tmpHull->size*sizeof(Ray*));
    for(int j=0; j<tmpHull->size; j++) {
      newHull->crits[j] = tmpHull->crits[j];
    }
    free(tmpHull);

    //actually use the ConvexHull class and return
//    l = hull;
//    return tmpHull->size;
    return newHull;
//    return new ConvexHull(hull, tmpHull->size);
}

int buildConvexHull(Point* l, int numPoints) {
    //The temporary list that contains the convex hull points (including the artificial points)
    Point* tmpHull = new Point[numPoints + 2];
    int hullSize = 0; //Denotes the current size of the convex hull

    Point* s = new Point(l[0].getX(), -1);
    Point* e = new Point(l[numPoints - 1].getX(), -1);

    tmpHull[0] = *s; // append the first point
    hullSize += 1;

    for (int i = 0; i < numPoints; i++) {
        //remove all points that violate the convex hull property
        while (hullSize >= 2 && //and the sequence of the last two points with the selected point does not make a counter-clockwise turn
               (1LL*(tmpHull[hullSize - 1].getY() - tmpHull[hullSize - 2].getY()) * (l[i].getX() - tmpHull[hullSize - 1].getX())
                - 1LL*(tmpHull[hullSize - 1].getX() - tmpHull[hullSize - 2].getX()) * (l[i].getY() - tmpHull[hullSize - 1].getY()))
               <= 0LL) {
            //Remove the last point from tmpHull
//            tmpHull[hullSize] = nullptr;
            hullSize -= 1;
        }
        //append the last selected point to
        tmpHull[hullSize] = l[i];
        hullSize += 1;
    }
    //append the last point
    tmpHull[hullSize] = *e;
    hullSize += 1;

    //optimise array size and require dynamically allocated array, because of sharing
    Point* hull = new Point[hullSize];
    for (int j = 0; j < hullSize; j++) {
        hull[j] = tmpHull[j];
    }

    delete[] tmpHull;

    //actually use the ConvexHull class and return
    l = hull;
    return hullSize;
//    return new ConvexHull(hull, hullSize);
}


/*
Hull* buildConvexHullPar(Point* l, Ray* crits, int numPoints) {
    //The temporary list that contains the convex hull points (including the artificial points)
    Hull* tmpHull = (Hull*)malloc(sizeof(Hull));
    tmpHull->points = (Point*)malloc((numPoints)*sizeof(Point));
    tmpHull->crits = (Ray**)malloc((numPoints)*sizeof(Ray*));
    tmpHull->size = 0; //Denotes the current size of the convex hull


    Point* s = &l[0];
    Point* e = &l[numPoints-1];

    tmpHull->points[0] = *s; // append the first point
    tmpHull->crits[0] = crits;
    tmpHull->size += 1;

    for (int i = 0; i < numPoints; i++) {
        //remove all points that violate the convex hull property
        while (tmpHull->size >= 2 && //and the sequence of the last two points with the selected point does not make a counter-clockwise turn
               (1LL*(tmpHull->points[tmpHull->size - 1].getY() - tmpHull->points[tmpHull->size - 2].getY()) * (l[i].getX() - tmpHull->points[tmpHull->size - 1].getX())
                - 1LL*(tmpHull->points[tmpHull->size - 1].getX() - tmpHull->points[tmpHull->size - 2].getX()) * (l[i].getY() - tmpHull->points[tmpHull->size - 1].getY()))
               <= 0LL) {
            //Remove the last point from tmpHull
//            tmpHull[tmpHull->size] = nullptr;
            tmpHull->size -= 1;
        }
        //append the last selected point to
        tmpHull->points[tmpHull->size] = l[i];
        tmpHull->crits[tmpHull->size] = &crits[i];
        tmpHull->size += 1;
    }
    //append the last point
//    tmpHull->points[tmpHull->size] = *e;
//    tmpHull->crits[tmpHull->size] = &crits[numPoints - 1];
//    tmpHull->size += 1;

    //optimise array size and require dynamically allocated array, because of sharing
//    Point* hull = new Point[tmpHull->size];
//    for (int j = 0; j < tmpHull->size; j++) {
//        hull[j] = tmpHull[j];
//    }

//    delete[] tmpHull;

    //actually use the ConvexHull class and return
//    l = hull;
//    return tmpHull->size;
    return tmpHull;
//    return new ConvexHull(hull, tmpHull->size);
}
*/
