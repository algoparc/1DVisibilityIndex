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


#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


class Ray {
  public:
  float slope;
  float icp;
  float x;
  float y;

  inline float dualSlope() {
    return x;
  }
  inline float dualX() {
    return slope;
  }
  inline float dualY() {
    return -icp;
  }
  inline float dualIcp() {
    return -y;
  }
};

struct Line {
  float slope;
  float icp;
};

struct Peak {
  int relCount;
  int below;
  int above;
  float pi;
};

void vis1D(int* xvals, int* peaks, int* visIdx, int SIZE, int numThreads);

float calcSlope(int x1, int y1, int x2, int y2);

float calcIcp(int x, int y, float slope);
