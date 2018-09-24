# 1D Total Visibility-Index Calculation

Developed and implemented by the <a href=http://algoparc.ics.hawaii.edu/>Algorithms and Parallel Computing (AlgoPARC) research group</a> at the University of Hawai'i at Manoa.

This is an implementation of the algorithms and datastructures used and described in the papers:

<a href=https://algoparc.ics.hawaii.edu/pubs/17-alenex.pdf>P. Afshani, M. de Berg, H. Casanova, B. Karsin, C. Lambrechts, N. Sitchinava, and C. Tsirogiannis. An Efficient Algorithm for the 1D Total Visibility-Index Problem. In Proceedings of the 19th Meeting on Algorithm Engineering and Experiments (ALENEX), pages 218-231, 2017</a>.

<a href=https://algoparc.ics.hawaii.edu/pubs/18-jea.pdf>P. Afshani, M. de Berg, H. Casanova, B. Karsin, C. Lambrechts, N. Sitchinava, and C. Tsirogiannis. An Efficient Algorithm for the 1D Total Visibility-Index Problem and its Parallelization. In the Journal of Experimental Algorithms (JEA), Volume 23(2): 2.3:1-2.3:23 (2018)</a>. 


## Algorithms implemented in the code

The code contains three implementations that each solve the 1D total visibility-index problem.  They all use the divide-and-conquer method described in the paper, but each solve the line segment intersection subproblem in a different way.

PLANE SWEEP (sweep):
The 'sweep' method is only implemented sequentially and, overall, performs the worst.  It uses a plane sweep approach to solving the segment intersection problem

ARRAY-OF-TREES (visAoT):
The visAoT method uses the array-of-trees persistent data structure solve the segment intersection problem.  The implementation can be executed in parallel and performs quite well.  

LINEAR PARALLEL (LinPar):
The 'linPar' method uses the succinct datastructure that uses bits to create a persistent data structure that we use to solve the segment intersection problem.  It uses only O(N) space, and, overall, performs better than the other two approaches.

## Using the code

To build the executables for each algorithm, simply go to the ./src/ subdirectory and use:
> make

This requires that a g++ compiler is installed and available.  Additionally, the AoT and LinPar methods require that openMP is installed on the system.  Once built, simply run them with the following usage:

- plane sweep
> ./sweep <size (2^x)>

- array-of-trees
> ./visAoT <size (2^x)> <# threads>

- linPar
> ./linPar <size (2^x)> <# threads>

Each executable creates a random synthetic dataset of the given size and computes the total visibility-index on it.  It repeats this a number of times and outputs the average execution time.  The number of runs is given by the ITERS declaration at the top of each main file.  To change the type of synthetic dataset, simply change the "buildRandom()" function call to one that builds the desired dataset (buildFlat() or buildParabolic()) and re-build the executable.

## Code organization

The code is organized into source code and a couple of sample scripts.  Below is an overview of the source code organization in the 'src/' directory:

- Makefile: builds executables for all three methods.
- sweep.cpp: the main function that runs the planesweep method.
- mainAoT.cpp: the main function that runs the AoT method.
- linPar.cpp: the main function that runs the LinPar method.
- buildSynthetic.h: contains code to create synthetic datasets.
- converxhull.h/.cpp: files to compute the convex hull, used to update critical rays.
- vis-succinct.h/cpp: functions used by the linPar method to compute segment intersection.
- visAoT.h/.cpp: functions used by the AoT method to compute segment intersection.
- SimpleRBI.h/.cpp: functions used by the sweep method to compute segment intersection.
- chazel/: directory containing helper functions used by linPar for bit manipulation and packing.
- data/: directory containing datatypes used by the AoT method datastructures.

There are also some scripts included in the 'scripts/' directory.  They simply run each implementation for different input sizes and number of parallel execution threads.


### Questions

For questions about this code or if you have any issues building or running it, please contact Ben Karsin at karsin@hawaii.edu.
