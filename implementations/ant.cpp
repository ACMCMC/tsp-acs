#include "../headers/tsp.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

long unsigned int Ant::_index = 0;

Ant::Ant(Node& start, std::vector<Node> nodes) : currentNode(start), visitedNodes({start}), unvisitedNodes(nodes)
{
  unvisitedNodes.erase(unvisitedNodes.begin() + start.getInternalId());
  //std::cout << "Ant " << _id << ": Created at node " << start.getHumanId() << std::endl;
}

void Ant::move(blaze::DynamicMatrix<double> &pheromoneMatrix, blaze::DynamicMatrix<double> &distanceMatrix)
{
  double exploitProbability = 0.5;

  // Calculate the probability of moving to each node in the unvisited nodes list
  blaze::DynamicVector<double> probabilities(unvisitedNodes.size());
  double sum = 0;
  for (long unsigned int i = 0; i < unvisitedNodes.size(); i++)
  {
    Node node = unvisitedNodes[i];
    double pheromone = pheromoneMatrix(currentNode.getInternalId(), node.getInternalId());
    double distance = distanceMatrix(currentNode.getInternalId(), node.getInternalId());
    double probability = pow(pheromone, TSPConstants::alpha) * pow(1.0 / distance, TSPConstants::beta);
    probabilities[i] = probability;
    sum += probability;
  }

  // Normalize the probabilities
  probabilities /= sum;

  // Decide whether to exploit or explore based on exploitProbability
  double random = (double)rand() / RAND_MAX;
  if (random < exploitProbability)
  {
    // Exploit
    // Choose the node with the highest probability
    long unsigned int maxIndex = 0;
    double maxProbability = probabilities[0];
    for (long unsigned int i = 1; i < probabilities.size(); i++)
    {
      if (probabilities[i] > maxProbability)
      {
        maxIndex = i;
        maxProbability = probabilities[i];
      }
    }

    // Move to that node
    currentNode = unvisitedNodes[maxIndex];
    visitedNodes.push_back(currentNode);
    unvisitedNodes.erase(unvisitedNodes.begin() + maxIndex);
  }
  else
  {
    // Explore
    // Choose a random node
    double random = (double)rand() / RAND_MAX;
    double sum = 0;
    long unsigned int index = 0;
    for (long unsigned int i = 0; i < probabilities.size(); i++)
    {
      sum += probabilities[i];
      if (sum > random)
      {
        index = i;
        break;
      }
    }

    // Move to that node
    currentNode = unvisitedNodes[index];
    visitedNodes.push_back(currentNode);
    unvisitedNodes.erase(unvisitedNodes.begin() + index);
  }
  //std::cout << "Ant " << _id << ": Moved to node " << currentNode.getHumanId() << std::endl;
}

void Ant::localUpdatePheromone(blaze::DynamicMatrix<double> &pheromoneMatrix)
{
  // Formula: pheromone = (1 - phi) * pheromone + phi * tau0
  long unsigned int node1 = visitedNodes[visitedNodes.size() - 2].getInternalId();
  long unsigned int node2 = visitedNodes[visitedNodes.size() - 1].getInternalId();
  double newPheromoneVal = (1 - TSPConstants::phi) * pheromoneMatrix(node1, node2) + TSPConstants::phi * TSPConstants::tau0;
  pheromoneMatrix(node1, node2) = newPheromoneVal;
  pheromoneMatrix(node2, node1) = newPheromoneVal;
}

double Ant::getTourLength(blaze::DynamicMatrix<double> &distanceMatrix)
{
  double tourLength = 0;
  for (long unsigned int i = 0; i < visitedNodes.size() - 1; i++)
  {
    tourLength += distanceMatrix(visitedNodes[i].getInternalId(), visitedNodes[i + 1].getInternalId());
  }
  return tourLength;
}

std::vector<Node> Ant::getVisitedNodes()
{
  return visitedNodes;
}

blaze::DynamicVector<long unsigned int> Ant::getVisitedNodesAsVector()
{
  blaze::DynamicVector<long unsigned int> visitedNodesAsVector(visitedNodes.size());
  for (long unsigned int i = 0; i < visitedNodes.size(); i++)
  {
    visitedNodesAsVector[i] = visitedNodes[i].getInternalId();
  }
  return visitedNodesAsVector;
}

/**
 * @brief Ant destructor
 */
Ant::~Ant() {
  //std::cout << "Ant " << _id << ": Destroyed" << std::endl;
}