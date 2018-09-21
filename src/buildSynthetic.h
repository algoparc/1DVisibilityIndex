
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

#include<stdlib.h>
#include<time.h>

void buildRandom(int size, int* peaks) {
	srand(time(NULL));
  for(int i=0; i<size; i++)
    peaks[i] = (int)(rand()%1000000);
}

void buildFlat(int size, int* peaks) {
  for(int i=0; i<size; i++)
	  peaks[i] = 1;

}

void buildParabolic(int size, int* peaks) {
  for(int i=0; i<size; i++)
	  peaks[i] = i*i;

}
