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
    static constexpr double alpha = 1;      // Importance of pheromone
    static constexpr double beta = 1;       // Importance of distance
};

class Node
{
public:
    Node(int id, double x, double y);
    ~Node();
    double getX() const;
    double getY() const;
    int getInternalId() const;
    int getHumanId() const;
private:
    double x;
    double y;
    int id; // Id in the distance matrix, not the real id (which starts at 1)
};

class TSPStatement
{
public:
    void read(const char *filename);
    int getDimension() const;
    double getDistance(int i, int j) const;
    void solve_aco();

private:
    int dimension;
    std::string name;
    int bestKnown;
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
    blaze::DynamicVector<int> getVisitedNodesAsVector();
private:
    Node currentNode;
    std::vector<Node> visitedNodes;
    std::vector<Node> unvisitedNodes;
    const int _id = _index++;
    static int _index;
};

#endif