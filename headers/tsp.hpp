#ifndef TSP_HPP
#define TSP_HPP

#include <array>
#include <string>
#include <blaze/Blaze.h>
#include <chrono>

class TSPConstants {
public:
    static double rho;      // Evaporation rate
    static double phi;      // Pheromone decay rate
    static double tau0;       // Initial pheromone
    static constexpr double minPheromone = 10e-10;       // Minimum pheromone value
    static constexpr double maxPheromone = 1.5;       // Maximum pheromone value
    static long unsigned int nAnts;      // Number of ants
    static double alpha;      // Importance of pheromone
    static double beta;       // Importance of distance
    static constexpr auto timeout = std::chrono::minutes{3}; // Timeout for the algorithm
    static constexpr int localSearchEnabled = 1; // Enable local search
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
    void createDistanceMatrix();
    long unsigned int getBestKnown() const;
    int getBestCost() const;
    std::string getName() const;
    blaze::DynamicVector<long unsigned int, false> getBestPath() const;
    void localSearch3Opt();
    void writeSolution(const char* filename);
    void localSearch2Opt();
private:
    long unsigned int dimension;
    std::string name;
    std::string comment;
    long unsigned int bestKnown;
    blaze::DynamicMatrix<double> distanceMatrix;
    std::vector<Node> nodes;
    blaze::DynamicVector<long unsigned int> bestPath;
    int bestCost = std::numeric_limits<int>::max();
};

/**
 * @brief An Ant is a class that represents an ant in the Ant Colony Optimization algorithm.
 */
class Ant
{
public:
    Ant(Node& start, std::vector<Node> nodes);
    ~Ant();
    void move(blaze::DynamicMatrix<double> &precalculatedTargetMatrix, double exploitProbability);
    void offlineUpdatePheromone(blaze::DynamicMatrix<double> &pheromoneMatrix, blaze::DynamicMatrix<double> &distanceMatrix);
    void localUpdatePheromone(blaze::DynamicMatrix<double> &pheromoneMatrix);
    int getSolutionLength(blaze::DynamicMatrix<double> &distanceMatrix);
    std::vector<Node> getVisitedNodes();
    blaze::DynamicVector<long unsigned int> getSolutionAsVector();
private:
    Node currentNode;
    std::vector<Node> visitedNodes;
    std::vector<Node> unvisitedNodes;
    const long unsigned int _id = _index++;
    static long unsigned int _index;
};

#endif