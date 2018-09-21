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


#ifndef DUALRAY_H
#define	DUALRAY_H

//=====Forward Declared Dependencies============================================


//=====Included Dependencies====================================================
//#include "Point.h"
//#include "SliceCell.h"

//=====Actual Class=============================================================

class DualRay {
    //=====Variables============================================================
public:
    /** The a and b from y = a x + b */
    float a, b;

    /** Denotes whether this dual ray is based on a left or right oriented critical ray */
    bool leftOriented;

    float endX, endY;

    //=====Constructors=============================================================
public:
    /**
     * Constructor when the internal values are not yet computed
     * 
     * @param s The starting point in the primal plane is the starting point of the critical ray
     * @param c The corresponding cell of the critical ray
     * @param r The critical ray needed for the dual ray
     */
//    DualRay(CriticalRay* r, SliceCell* c);


    /**
     * Constructor for the computed values
     * 
     * @param x1 The x-coordinate of the starting point of the critical ray
     * @param y1 The y-coordinate of the starting point of the critical ray
     * @param x2 The x-coordinate of the other point of the critical ray
     * @param y2 The y-coordinate of the other point of the critical ray
     * @param c The corresponding cell
     * @param a The a from the formula y = a x + b
     * @param b The b from the formula y = a x + b
     * @param isL Denotes whether it is based on a left or tight oriented critical ray
     * @pre {@code dx != 0}
     * @pre {@code c != nullptr}
     */
//    DualRay(int x1, int y1, int x2, int y2, SliceCell* c, int a, int b, bool isL);
//    DualRay(int x1, int y1, int x2, int y2, int a, int b, bool isL);


    /**
     * Copy constructor
     * 
     * @param orig The original to copy from
     */
//    DualRay(const DualRay& orig);

    /**
     * The move constructor
     * 
     * @param orig The original to move from
     * @return The address to the moved object
     */

    DualRay(float a);
    void setA(float a);

    DualRay& operator=(const DualRay& orig);

    bool operator<(const DualRay& a);

    /**
     * The destructor
     */
    virtual ~DualRay();

    //=====Getters and Setters==================================================
public:
    /**
     * Returns the corresponding cell
     * 
     * @return {@code cell}
     */
//    SliceCell* getCell();

    /**
     * Returns the {@code a} of the equation {@code y = a x + b}
     * 
     * @return {@code a}
     */
    float getA();

    /**
     * Returns the {@code b} of the equation {@code y = a x + b}
     * 
     * @return {@code b}
     */
    float getB();

    /**
     * Returns the x-coordinate of the starting point of the critical ray
     * 
     * @return {@code startX}
     */
//    int getStartX();

    /**
     * Returns the y-coordinate of the starting point of the critical ray
     * 
     * @return {@code startY}
     */
//    int getStartY();

    /**
     * Returns the x-coordinate of the other point of the critical ray
     * 
     * @return {@code otherX}
     */
//    int getOtherX();

    /**
     * Returns the y-coordinate of the other point of the critical ray
     * 
     * @return {@code otherY}
     */
//    int getOtherY();

    /**
     * Returns the numerator of the slope of the critical ray
     * 
     * @return {@code num}
     */
//    int getNum();

    /**
     * Returns the denominator of the slop of the critical ray
     * 
     * @return {@code denom}
     */
//    int getDenom();

    /**
     * Returns whether this dual ray is based on a left or right critical ray.
     * 
     * @return {@code leftOriented}
     */
    bool isLeftOriented();

    float getEndX();
    float getEndY();

    int Overlap;
    int Above;
    int Under;
    DualRay* immBelow;
    int label;

    bool rAbove(DualRay* b);

    //=====Methods==============================================================

    /**
     * Returns whether this dual rays starts before the other. Note that this incorporates the endpoints and is not only based on starting points. 
     * As there are no endpoints we consider that the scope starts at -infinity and continues to infinity. 
     * This ordering is only horizontal (on x)
     * 
     * @param ray The dual ray to test against
     * @return {@code true} if the rays has an endpoint or a start point that is before that of the other
     */
//    bool startsBeforeOther(DualRay* ray);

    /**
     * Returns whether the starring point of the dual ray lies before the starting point of the other dual ray. 
     * In essence it returns whether the starting point of {@code this} lies before the starting point of {@code ray}.
     * Note that it is a strict operator: so when starting points are the same, returns false for both. 
     * This ordering is only horizontal (on x)
     * 
     * @param ray The other ray to test against
     * @pre both have a starting point
     * @return {@code true} if the starting point of {@code this} is before the starting point of {@code ray}, return {@code false} otherwise
     * note that this means that if that have the same starting point we get false
     */
//    bool startingPointBeforeOther(DualRay* ray);
};

#endif	/* DUALRAY_H */

