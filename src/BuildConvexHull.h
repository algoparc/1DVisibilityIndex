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


/** 
 * This class contains the algorithm to build a convex hull based from a list of points. 
 * It uses an adopted version of the algorithm of Adrew's Monotone Chain algorithm.
 * The running time is O(n) (no sorting needed and therefore no O(n log(n)))
 * 
 */

#ifndef BUILDCONVEXHULL_H
#define	BUILDCONVEXHULL_H

//=====Forward Declared Dependencies============================================


//=====Included Dependencies====================================================
//#include "../data/ConvexHull.h"
//#include "../data/Point.h"
#include "convexhull.h"

//=====Actual Class=============================================================

class BuildConvexHull {
    //=====Variables============================================================

    //=====Constructors=========================================================
public:
    /**
     * Constructor for convex hull builder
     */
    //BuildConvexHull();

    /**
     * Copy constructor for convex hull builder
     * 
     * @param orig
     */
    //BuildConvexHull(const BuildConvexHull& orig);

    /**
     * Move constructor for convex hull builder
     * 
     * @param orig
     * @return 
     */
    //BuildConvexHull& operator=(const BuildConvexHull& orig);

    /**
     * Destructor for convex hull builder
     */
    //virtual ~BuildConvexHull();

    //=====Getters and Setters==================================================

    //=====Methods==============================================================
public:
    /**
     * Returns the upper convex hull for the sorted input of points
     * 
     * @param l A sorted array to build the convex hull out of
     * @param numPoints The number of points in the input
     * @return The upper convex hull. Additionally it also adds a first and last element to the hull. 
     *      The first element is below the first point of the convex hull at -1 and symmetrically for the last point. 
     */
//    ConvexHull* buildConvexHull(Point** l, int numPoints);
};

Hull* buildAndSaveHull(Ray* crits, int numPoints);

int buildConvexHull(Point* l, int numPoints);

Hull* buildConvexHullPar(Point* l, Ray* crits, int numPoints);


#endif	/* BUILDCONVEXHULL_H */

