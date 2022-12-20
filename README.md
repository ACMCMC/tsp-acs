# tsp-acs

Implementation of the **Ant Colony System** for solving the **Travelling Salesman Problem**.

## Algorithm

I used a modified version of the Ant Colony System (ACS). I started implementing that algorithm, and then I modified
it to perform a *2-opt* local search, followed by a *3-opt* local search.
Also, I choose to save the most promising solutions in a list. After the 20%
first part of the run, these solutions are optimized using the local search techniques
mentioned above. This allows to perform more exploration at first, and then, focus
on improving the best solutions. These solutions are the only ones allowed to spread pheromone.
Moreover, I implemented alpha so that it gradually decreases, using a truncated logarithmic function.

## Implementation

I implemented the algorithm using C++, and the highly-performant mathematics library Blaze. For the compilation, I delegated the task to CMake.
I chose to precompute some parameters, such as the heuristics matrix that does not change in the evolution of the algorithm.

## Choice of hyperparameters

I used Bayesian Optimization to choose the hyperparameters for the algorithm.
To get my choice of hyperparameters, I used a Python script, included with my
code, that uses the `bayesian_optimization` library to perform repeated runs of the
algorithm, with the objective of finding the best choice of the hyperparameters

## Results

Average error (with respect to the optimum solution) of 3%.
