#include "../headers/tsp.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

auto timeout = std::chrono::minutes{3}; // Timeout for the algorithm

long unsigned int TSPStatement::getDimension() const
{
    return dimension;
}

double TSPStatement::getDistance(long unsigned int i, long unsigned int j) const
{
    return distanceMatrix(i, j);
}

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage, std::string problemName)
{
    // Taken from https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%s: %3d%% [%.*s%*s]", problemName.c_str(), val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

void TSPStatement::createDistanceMatrix()
{
    long unsigned int dimension = getDimension();
    distanceMatrix = blaze::DynamicMatrix<double>(dimension, dimension);
    for (long unsigned int i = 0; i < dimension; i++)
    {
        for (long unsigned int j = 0; j < i; j++)
        {
            // Calculate the Euclidean distance between nodes i and j
            double x1 = nodes[i].getX();
            double y1 = nodes[i].getY();
            double x2 = nodes[j].getX();
            double y2 = nodes[j].getY();
            double distance = round(sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)));
            distanceMatrix(i, j) = distance;
            distanceMatrix(j, i) = distance; // Symmetric matrix
        }
    }
    for (long unsigned int i = 0; i < dimension; i++)
    {
        distanceMatrix(i, i) = 0; // Diagonal is 0
    }
}

void TSPStatement::read(const char *filename)
{
    std::ifstream file(filename);

    if (!file)
    {
        perror("Error opening file");
        return;
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string word;
        if (line.find("NAME") != std::string::npos)
        {
            iss >> word >> word;
            if (word == ":") // if there is a colon, read the next word
            {
                iss >> word;
            }
            name = word;
        }
        if (line.find("BEST_KNOWN") != std::string::npos)
        {
            iss >> word >> word;
            if (word == ":") // if there is a colon, read the next word
            {
                iss >> word;
            }
            bestKnown = std::stoi(word);
        }
        if (line.find("DIMENSION") != std::string::npos)
        {
            iss >> word >> word;
            if (word == ":") // if there is a colon, read the next word
            {
                iss >> word;
            }
            dimension = std::stoi(word);
            // distance.resize(dimension);
            // distance.fill({});
        }
        if (line.find("NODE_COORD_SECTION") != std::string::npos)
        {
            break; // we are done with the header
        }
    }
    // Now, read the distance matrix
    while (std::getline(file, line) && line.find("EOF") == std::string::npos)
    {
        std::istringstream iss(line);
        std::string word;
        iss >> word;
        long unsigned int i = std::stoi(word);
        iss >> word;
        double x = std::stod(word);
        iss >> word;
        double y = std::stod(word);
        Node node(i, x, y);
        nodes.insert(nodes.end(), node);
    }
    std::cout << "Name: " << name << std::endl;
    std::cout << "Best known: " << bestKnown << std::endl;
    std::cout << "Dimension: " << dimension << std::endl;

    createDistanceMatrix();

    std::cout << "Distance matrix: " << distanceMatrix.rows() << "x" << distanceMatrix.columns() << std::endl;

    file.close();
}

void offlineUpdatePheromone(blaze::DynamicMatrix<double> &pheromoneMatrix, blaze::DynamicMatrix<double> &distanceMatrix, blaze::DynamicVector<long unsigned int> &visitedNodes, double cost)
{
    pheromoneMatrix = (1 - TSPConstants::rho) * pheromoneMatrix; // Evaporation in all edges
    // Update the pheromone matrix to sum the edges of the best tour
    long unsigned int node1; // Will be written to in the first iteration
    long unsigned int node2 = visitedNodes[0];
    for (long unsigned int i = 1; i < visitedNodes.size(); i++)
    {
        node1 = node2; // node2 from the previous iteration
        node2 = visitedNodes[i];
        double newAmount = pheromoneMatrix(node1, node2) + TSPConstants::rho / cost;
        //double newAmount = pheromoneMatrix(node1, node2) + 1.0 / cost;
        pheromoneMatrix(node1, node2) = newAmount;
        pheromoneMatrix(node2, node1) = newAmount;
    }
}

void TSPStatement::solve_aco()
{
    // Solution using the Ant Colony Optimization algorithm

    // Parameters
    long unsigned int nAnts = 10;                                             // Number of ants
    auto finish = std::chrono::system_clock::now() + timeout; // Stop after 3 minutes
    double deadlineInSeconds = std::chrono::duration<double>(timeout).count();

    double exploitProbability = 0.2;

    // Initialize pheromone matrix
    long unsigned int dimension = getDimension();
    double tau = TSPConstants::tau0;
    blaze::DynamicMatrix<double> pheromone(dimension, dimension, tau); // All elements are tau0, even the diagonal, but it doesn't matter (it will never be used)

    // Initialize ants
    // Randomly position nAnts artificial ants on nNodes nodes
    std::vector<Ant> ants;
    // Main loop
    for (long unsigned int iteration = 0; std::chrono::system_clock::now() < finish; iteration++)
    {
        if (iteration % 100 == 0) // Perform "expensive" updates every 100 iterations
        {
            // Print progress
            // Only print progress every 100 iterations
            double remainingSeconds = std::chrono::duration<double>(finish - std::chrono::system_clock::now()).count();
            double percentage = 1.0 - (remainingSeconds / deadlineInSeconds);
            printProgress(percentage, name);

            exploitProbability = 0.3 + 0.4 * percentage; // Start with 30% exploitation, then increase to 70%
        }

        // Generate nAnts ants
        ants.clear();
        for (long unsigned int j = 0; j < nAnts; j++)
        {
            // Choose a random node
            long unsigned int node = rand() % dimension;
            Node n = nodes[node];

            // Create an ant on that node
            Ant ant(n, nodes);
            // Add the ant to the ant list
            ants.push_back(Ant(n, nodes));
        }

        // Construct solutions
        // Visit all nodes
        for (long unsigned int j = 0; j < dimension - 1; j++)
        {
            // For each ant
            for (long unsigned int k = 0; k < nAnts; k++)
            {
                // Move the ant
                Ant &ant = ants[k];
                ant.move(pheromone, distanceMatrix, exploitProbability);
                ant.localUpdatePheromone(pheromone);
            }
        }

        // Update best solution
        for (long unsigned int j = 0; j < nAnts; j++)
        {
            Ant &ant = ants[j];
            double cost = ant.getSolutionLength(distanceMatrix);
            if (cost < bestCost)
            {
                bestCost = cost;
                bestPath = ant.getSolutionAsVector();

                std::cout << "\r" << std::endl << "[Iteration " << iteration << "] Best cost: " << bestCost << std::endl;
            }
        }

        // Offline pheromone update
        offlineUpdatePheromone(pheromone, distanceMatrix, bestPath, bestCost);

        // Print pheromone matrix
        // std::cout << "Pheromone matrix: " << std::endl;
        // std::cout << pheromone << std::endl;
    }
}

blaze::DynamicVector<long unsigned int> TSPStatement::getBestPath() const
{
    return bestPath;
}

long unsigned int TSPStatement::getBestKnown() const
{
    return bestKnown;
}

int TSPStatement::getBestCost() const
{
    return bestCost;
}

std::string TSPStatement::getName() const
{
    return name;
}