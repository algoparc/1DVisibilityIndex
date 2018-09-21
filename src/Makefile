
#Copyright 2016-2018 Ben Karsin, Nodari Sitchinava
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

all: linPar visAoT sweep


visAoT: BuildConvexHull.h BuildConvexHull.cpp AOT.h visAoT.cpp visAoT.h mainAoT.cpp convexhull.cpp convexhull.h
	g++ -Ofast -ggdb -std=c++11 -o visAoT mainAoT.cpp visAoT.cpp convexhull.cpp BuildConvexHull.cpp -fopenmp 

linPar: BuildConvexHull.h BuildConvexHull.cpp AOT.h vis-succinct.cpp vis-succinct.h chazel/params.h linPar.cpp convexhull.cpp convexhull.h chazel/chazelBits.cpp buildSynthetic.h
	g++ -Ofast -ggdb -std=c++11 -o linPar linPar.cpp vis-succinct.cpp convexhull.cpp BuildConvexHull.cpp -fopenmp 

sweep: DualRay.o PTree.o BuildConvexHull.h BuildConvexHull.cpp SimpleRBI.h SimpleRBI.cpp sweep.cpp convexhull.cpp convexhull.h 
	g++ -O3 -std=c++11 -o sweep sweep.cpp SimpleRBI.cpp convexhull.cpp BuildConvexHull.cpp DualRay.o PTree.o 

DualRay.o: data/DualRay.cpp
	g++ -Ofast -std=c++11 -c data/DualRay.cpp

PTree.o: data/PTree.cpp
	g++ -Ofast -std=c++11 -c data/PTree.cpp

clean:
	rm -rf *.o visAoT linPar sweep
