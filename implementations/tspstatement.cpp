#include "../headers/tsp.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

int TSPStatement::getDimension() const
{
    return dimension;
}

double TSPStatement::getDistance(int i, int j) const
{
    return distanceMatrix(i, j);
}

void TSPStatement::createDistanceMatrix()
{
    int dimension = getDimension();
    distanceMatrix = blaze::DynamicMatrix<double>(dimension, dimension);
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < i; j++)
        {
            // Calculate the Euclidean distance between nodes i and j
            double x1 = nodes[i].getX();
            double y1 = nodes[i].getY();
            double x2 = nodes[j].getX();
            double y2 = nodes[j].getY();
            double distance = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
            distanceMatrix(i, j) = distance;
            distanceMatrix(j, i) = distance; // Symmetric matrix
        }
    }
    for (int i = 0; i < dimension; i++)
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
        int i = std::stoi(word);
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

void TSPStatement::solve_aco()
{
    // Solution using the Ant Colony Optimization algorithm

    // Print distance matrix
    std::cout << "Distance matrix:" << std::endl;
    std::cout << distanceMatrix << std::endl;

    // Parameters
    int nAnts = 2;        // Number of ants
    int nIterations = 5; // Number of iterations

    // Initialize pheromone matrix
    int dimension = getDimension();
    double tau = TSPConstants::tau0;
    blaze::DynamicMatrix<double> pheromone(dimension, dimension, tau); // All elements are tau0, even the diagonal, but it doesn't matter (it will never be used)
    
    // Initialize best solution
    blaze::DynamicVector<int> bestSolution(dimension);
    double bestCost = std::numeric_limits<double>::max();

    // Initialize ants
    // Randomly position nAnts artificial ants on nNodes nodes
    std::vector<Ant> ants;
    // Main loop
    for (int i = 0; i < nIterations; i++)
    {
        std::cout << std::endl << "Iteration " << i << std::endl;
        // Generate nAnts ants
        ants.clear();
        for (int i = 0; i < nAnts; i++)
        {
            // Choose a random node
            int node = rand() % dimension;
            Node n = nodes[node];

            // Create an ant on that node
            Ant ant(n, nodes);
            // Add the ant to the ant list
            ants.push_back(Ant(n, nodes));
        }

        // Construct solutions
        // Visit all nodes
        for (int j = 0; j < dimension - 1; j++)
        {
            // For each ant
            for (int k = 0; k < nAnts; k++)
            {
                // Move the ant
                Ant &ant = ants[k];
                ant.move(pheromone, distanceMatrix);
                ant.localUpdatePheromone(pheromone);
            }
        }

        // Update best solution
        for (int j = 0; j < nAnts; j++)
        {
            Ant &ant = ants[j];
            double cost = ant.getTourLength(distanceMatrix);
            if (cost < bestCost)
            {
                bestCost = cost;
                bestSolution = ant.getVisitedNodesAsVector();
            }
        }

        // Offline pheromone update
        

        std::cout << "Best cost: " << bestCost << std::endl;
        std::cout << "Best solution: " << bestSolution << std::endl;

        // Print pheromone matrix
        std::cout << "Pheromone matrix: " << std::endl;
        std::cout << pheromone << std::endl;
    }
}