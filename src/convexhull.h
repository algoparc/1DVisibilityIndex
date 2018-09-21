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


#include "visTypes.h"

class Point {
  public:
  float x;
  float y;

  float getX() {
    return x;
  }
  float getY() {
    return y;
  }
  Point(float inX, float inY) {
    x = inX;
    y = inY;
  }

  Point(){}

};

struct Hull {
  int size;
  Ray** crits;
};

int CreateConvexHull(int size, Point* p);

void updateCriticalRays(Ray* Rcrit, Ray* Lcrit, int step, int threads);

Hull* updateCriticalRaysReturnHull(Ray* Rcrit, Ray* Lcrit, int step, int threads);

Hull* updateCriticalRaysPar(Ray* Rcrit, Ray* Lcrit, Hull* lHull, Hull* rHull, int step, int threads);
