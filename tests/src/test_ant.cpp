#include "../../headers/tsp.hpp"

blaze::DynamicMatrix<double> distanceMatrix({{0, 4, 1, 3}, {4, 0, 2, 1}, {1, 2, 0, 5}, {3, 1, 5, 0}});
blaze::DynamicMatrix<double> pheromoneMatrix({{1,1,1,1}, {1,1,1,1}, {1,1,1,1}, {1,1,1,1}});

int testMove()
{
  Node start(1, 0, 0);
  Ant ant(start, {Node(2, 0, 0), Node(3, 0, 0), Node(4, 0, 0)});
  ant.move(distanceMatrix, pheromoneMatrix, 0.5);
  return 0;
}