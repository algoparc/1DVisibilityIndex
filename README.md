- Chazel's working for everything using arrays of booleans!
- Queries work for any time and key values and properly calculates counts

BUILDING STRUCT:

- Don't build prefixSums list yet, need to implement this...
- Don't build lookup tables yet

NEXT STEP:

- Instead of using bits, use ints and shift values into them
  - With this, 32 bits are saved to 1 int value
  - For each int, save a prefix sum count of all bits leading up to it

- Build lookup tables for every 32-bit combination?
- After finding prefix sum of leading ints, mask out necessary bits and use lookup table

EXPERIMENTS:

- Find out if lookup tables are worth it: Better to just scan the last 32 bits (or 16 bits if use short?)
- Compare build & query performance for different configs
- Compare build and query performance vs. AOT (seq. and parallel)

FUTURE STUFF:

- PARALLELIZE!!
- Read papers and think of applications of this, can we do everything AOT can do??
- What about using Hamming weight isntead of lookup table?
- popcnt primitive usable instead?  We can mask out unwanted bits and popcnt!
  - build bitmask as (1 << index) - 1
  - mask out unwanted bits as (mask AND int)
  - Then just use popcnt, hamming weight, or lookup table!

WORK SCHEDULE

TASK                           -     TIME NEEDED

1 Put code on bitbucket          -        DONE

2 Try bits vs. bools (seq)              

 - Build tree with bits       -        2 hr
 - Save block prefix sum      -        1.5 hr
 - total                      -        3.5 hr

3 Basic query with bits      -        3 hr

 - Try popcnt with mask       -        1 hr
 - Try lookup tables          -        1 hr
 - total                        -        5 hr

4 Basic parallel version         -        2 hr

5 Add mergepath for top levels   -        3 hr

6 Integrate into vis-code        -        2 hr

7 Get new results/graphs/etc.    -        4-5 hr

Priorities:

MUST DO:

4, 6, 7 - Total 9-10 hrs


Need to do:

2a, 2b, 3 - Total 6.5 hrs


Would like to:

3a, 3b  -  Total 2 hrs


Maybe if time:

5         - Total 3 hrs (maybe more)

Order of work: 

- DONE build repo 
- DONE Build tree with bits
- DONE Basic query, scan all bits for prefix sum (seq)

CHECKPOINT - Is it worth using bits vs. bools?
YES.  Using popcnt on all, 60-400% faster for query, slightly slower to build though

- DONE Save prefix sums while building
- DONE Use prefix sums for querying

CHECKPOINT - Prefix sums worth it?  bits vs. bools?
YES! Prefix sums reduces Q significantly!  

- Try popcount vs. scan
- Try lookup tables vs. popcnt or scan

Lowering priority to try these, is already quite fast with popcount and prefix sums list

CONTINUE

- DONE Make current fastest version parallel (BASIC)
- DONE First version: calculate prefix sums sequentially as well

CHECKPOINT -  How much faster is paralell tree building?  What about querying?

On 4-core (with HTT), Maximum speedup building tree: 3.23,
max speedup querying: 4.70

IF HAVE TIME

- DONE Further imrpove building parallelism with Mergepath for upper merge levels

CONTINUE 

- DONE Eval parallel speedup on cluster, better than AOT??

Not much speedup on cluster.  Need mergepath and better query parallelization?

Improved build time speedup.  Added mergepath and parallel prefix sum calculation.  Max speedup 10.4 with 16 threads on cluster

CONTINUE

- Integrate into vis code
- Get results on all hardware & datasets


