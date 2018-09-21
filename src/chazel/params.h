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


#include<stdint.h>

//#define WORDTYPE uint8_t
//#define WORDSIZE 8
//#define MAX_INT 0xff

//#define WORDTYPE uint16_t
//#define WORDSIZE 16
//#define MAX_INT 0xffff

#define WORDTYPE unsigned int
#define WORDSIZE 32
#define MAX_INT -1

//#define WORDTYPE unsigned long long
//#define WORDSIZE 64
//#define MAX_INT -1

//#define WORDTYPE bool
//#define WORDSIZE 1

#define THREADS 1
//#define QTHREADS 

//#define QNUM 1000000
#define QNUM 1
#define ITERS 1
