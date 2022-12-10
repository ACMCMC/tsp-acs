#ifndef TSP_HPP
#define TSP_HPP

#include <array>
#include <string>
#include <blaze/Blaze.h>

class TSPConstants {
public:
    static constexpr double rho = 0.1;      // Evaporation rate
    static constexpr double phi = 0.1;      // Pheromone decay rate
    static constexpr double tau0 = 1;       // Initial pheromone
    static constexpr double alpha = 5;      // Importance of pheromone
    static constexpr double beta = 1;       // Importance of distance
};

class Node
{
public:
    Node(long unsigned int id, double x, double y);
    ~Node();
    double getX() const;
    double getY() const;
    long unsigned int getInternalId() const;
    long unsigned int getHumanId() const;
private:
    double x;
    double y;
    long unsigned int id; // Id in the distance matrix, not the real id (which starts at 1)
};

class TSPStatement
{
public:
    void read(const char *filename);
    long unsigned int getDimension() const;
    double getDistance(long unsigned int i, long unsigned int j) const;
    void solve_aco();

private:
    long unsigned int dimension;
    std::string name;
    long unsigned int bestKnown;
    blaze::DynamicMatrix<double> distanceMatrix;
    std::vector<Node> nodes;

    void createDistanceMatrix();
};

/**
 * @brief An Ant is a class that represents an ant in the Ant Colony Optimization algorithm.
 */
class Ant
{
public:
    Ant(Node& start, std::vector<Node> nodes);
    ~Ant();
    void move(blaze::DynamicMatrix<double> &pheromoneMatrix,blaze::DynamicMatrix<double> &distanceMatrix);
    void offlineUpdatePheromone(blaze::DynamicMatrix<double> &pheromoneMatrix, blaze::DynamicMatrix<double> &distanceMatrix);
    void localUpdatePheromone(blaze::DynamicMatrix<double> &pheromoneMatrix);
    double getTourLength(blaze::DynamicMatrix<double> &distanceMatrix);
    std::vector<Node> getVisitedNodes();
    blaze::DynamicVector<long unsigned int> getVisitedNodesAsVector();
private:
    Node currentNode;
    std::vector<Node> visitedNodes;
    std::vector<Node> unvisitedNodes;
    const long unsigned int _id = _index++;
    static long unsigned int _index;
};

#endif