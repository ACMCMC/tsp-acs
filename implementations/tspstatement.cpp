#include "../headers/tsp.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#define TIME_NOT

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
        if (line.find("COMMENT") != std::string::npos)
        {
            iss >> word >> word;
            if (word == ":") // if there is a colon, read the next word
            {
                iss >> word;
            }
            comment = word;
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

void offlineUpdatePheromone(blaze::DynamicMatrix<double> &pheromoneMatrix, blaze::DynamicMatrix<double> &distanceMatrix, std::vector<blaze::DynamicVector<long unsigned int>> &visitedNodesVector, std::vector<double> &costs, double updateThreshold)
{
    pheromoneMatrix = (1 - TSPConstants::rho) * pheromoneMatrix; // Evaporation in all edges
    for (int i = 0; i < visitedNodesVector.size(); i++)
    {
        blaze::DynamicVector<long unsigned int> visitedNodes = visitedNodesVector.at(i);
        double cost = costs.at(i);
        if (cost > updateThreshold)
        {
            continue;
        }

        // Update the pheromone matrix to sum the edges of the best tour
        long unsigned int node1; // Will be written to in the first iteration
        long unsigned int node2 = visitedNodes[0];
        for (long unsigned int i = 1; i < visitedNodes.size(); i++)
        {
            node1 = node2; // node2 from the previous iteration
            node2 = visitedNodes[i];

            // double newAmount = pheromoneMatrix(node1, node2) + (TSPConstants::rho * pow(distanceMatrix.rows(), 10)) / cost; // rho * n / cost (temporal fix)
            double newAmount = pheromoneMatrix(node1, node2) + TSPConstants::rho / cost;
            pheromoneMatrix(node1, node2) = newAmount;
            pheromoneMatrix(node2, node1) = newAmount;
        }
    }
}

void TSPStatement::solve_aco()
{
    // Solution using the Ant Colony Optimization algorithm

    // Parameters
    auto timeout = TSPConstants::timeout;                     // Timeout for the algorithm
    auto finish = std::chrono::system_clock::now() + timeout; // Stop after 3 minutes
    double deadlineInSeconds = std::chrono::duration<double>(timeout).count();

    double exploitProbability = 0.5;
    double adaptedAlpha = 0.0;

    int lastLocalSearchImproved = 1; // 0 = no, 1 = yes

    std::vector<blaze::DynamicVector<long unsigned int>> bestPaths;
    std::vector<double> bestCosts;
    std::vector<double> originalACOCosts;
    std::vector<int> optimum;

#ifdef TIME
    std::chrono::duration<double> totalMovementsTime = std::chrono::duration<double>::zero();
    std::chrono::duration<double> totalLocalSearchTime = std::chrono::duration<double>::zero();
#endif

    // Initialize pheromone matrix
    long unsigned int dimension = getDimension();
    double tau0 = TSPConstants::tau0;
    double maxPheromone = TSPConstants::maxPheromone;
    double minPheromone = TSPConstants::minPheromone;
    blaze::DynamicMatrix<double> pheromone(dimension, dimension, tau0); // All elements are tau0, even the diagonal, but it doesn't matter (it will never be used)

    double savePathThreshold = std::numeric_limits<double>::max();
    double minACOLength = std::numeric_limits<double>::max();

    // Utility matrix
    blaze::DynamicMatrix<double> inverseDistanceMatrix = blaze::pow(distanceMatrix, -1.0);                // All elements are 1/distance
    blaze::DynamicMatrix<double> heuristicMatrix = blaze::pow(inverseDistanceMatrix, TSPConstants::beta); // All elements are (1/distance)^beta

    double lastPercentage = 0.0;

    // Initialize ants
    // Randomly position nAnts artificial ants on nNodes nodes
    std::vector<Ant> ants;
    // Main loop
    for (long unsigned int iteration = 0; std::chrono::system_clock::now() < finish; iteration++)
    {
        double remainingSeconds = std::chrono::duration<double>(finish - std::chrono::system_clock::now()).count();
        double percentage = 1.0 - (remainingSeconds / deadlineInSeconds);
        adaptedAlpha = TSPConstants::alpha * std::min(2.0, (-log(percentage)));

        exploitProbability = 0.5 + 0.2 * percentage; // Start with 50% exploitation, then increase to 70%

        if (percentage - lastPercentage >= 0.009)
        {
            // Print progress
            printProgress(percentage, name);
            lastPercentage = percentage;
        }

        // Generate nAnts ants
        ants.clear();
        for (long unsigned int j = 0; j < TSPConstants::nAnts; j++)
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
        // Calculate the target value of the function now, to avoid doing it in the inner loop
        blaze::DynamicMatrix<double> precalculatedTargetMatrix = (blaze::pow(pheromone, adaptedAlpha) % heuristicMatrix);

#ifdef TIME
        auto startMovements = std::chrono::system_clock::now(); // Stop after 3 minutes
#endif
        // Visit all nodes
        for (long unsigned int j = 0; j < dimension - 1; j++)
        {
            // For each ant
            for (long unsigned int k = 0; k < TSPConstants::nAnts; k++)
            {
                // Move the ant
                Ant &ant = ants[k];
                ant.move(precalculatedTargetMatrix, exploitProbability);
                ant.localUpdatePheromone(pheromone);
            }
        }
#ifdef TIME
        auto finishMovements = std::chrono::system_clock::now(); // Stop after 3 minutes
        totalMovementsTime += std::chrono::duration<double>(finishMovements - startMovements);
#endif

        int isThereANewBest = 0; // 0 = no, 1 = yes

        double averagePathLength = 0.0;
        // Update best solution
        for (long unsigned int j = 0; j < TSPConstants::nAnts; j++)
        {
            Ant &ant = ants[j];
            double cost = ant.getSolutionLength(distanceMatrix);
            if (cost < savePathThreshold)
            {
                // Does the same cost already exist?
                int alreadyExists = 0;
                for (long unsigned int k = 0; k < originalACOCosts.size(); k++)
                {
                    if (originalACOCosts[k] == cost)
                    {
                        alreadyExists++;
                        if (alreadyExists >= 2) // Save the same cost (not necessarily same solution) at most 2 times
                        {
                            break;
                        }
                    }
                }

                if (alreadyExists >= 2)
                {
                    continue;
                }
                isThereANewBest = 1;

                bestPaths.push_back(ant.getSolutionAsVector());
                bestCosts.push_back(cost);
                originalACOCosts.push_back(cost);
                optimum.push_back(0);

                //std::cout << "\r" << std::endl
                //          << "[Iteration " << iteration << "] Saving with cost: " << cost << std::endl;

                if (cost < minACOLength) // Update the best solution found by ACO
                {
                    minACOLength = cost;
                }
            }
            averagePathLength += cost;
        }
        averagePathLength /= TSPConstants::nAnts;
        savePathThreshold = std::min(savePathThreshold, (4.0 * minACOLength + averagePathLength) / 5.0); // Middle ground between the average and the best

        // Remove bad solutions
        for (long unsigned int j = 0; j < bestCosts.size(); j++)
        {
            if (bestCosts[j] > savePathThreshold)
            {
                bestPaths.erase(bestPaths.begin() + j);
                bestCosts.erase(bestCosts.begin() + j);
                originalACOCosts.erase(originalACOCosts.begin() + j);
                optimum.erase(optimum.begin() + j);
                j--;
            }
        }

        if (TSPConstants::localSearchEnabled && percentage > 0.2 && (isThereANewBest || lastLocalSearchImproved))
        {
#ifdef TIME
            auto startLocalSearch = std::chrono::system_clock::now(); // Stop after 3 minutes
#endif
            // Search all the best paths and costs
            for (int i = 0; i < bestPaths.size(); i++)
            {
                if (optimum.at(i) == 0)
                {
                    // Local search only if there is a new best solution
                    int localSearch2OptImproved = localSearch2Opt(bestPaths, bestCosts, i);
                    // std::cout << bestPath << std::endl;
                    int localSearch3OptImproved = localSearch3Opt(bestPaths, bestCosts, i);
                    // std::cout << bestPath << std::endl;
                    if (!(localSearch2OptImproved || localSearch3OptImproved))
                    {
                        optimum.at(i) = 1;
                    }
                    // std::cout << "[Iteration " << iteration << "] Best cost after local search (" << i << "): " << bestCosts.at(i) << std::endl;
                }
            }

            lastLocalSearchImproved = 0;
            for (int i = 0; i < optimum.size(); i++)
            {
                if (optimum.at(i) == 0)
                {
                    lastLocalSearchImproved = 1;
                    break;
                }
            }
#ifdef TIME
            auto endLocalSearch = std::chrono::system_clock::now(); // Stop after 3 minutes
            totalLocalSearchTime += std::chrono::duration<double>(endLocalSearch - startLocalSearch);
#endif
        }

        // std::cout << "Pheromone matrix: " << std::endl;
        // std::cout << pheromone << std::endl;
        //  Offline pheromone update
        offlineUpdatePheromone(pheromone, distanceMatrix, bestPaths, bestCosts, savePathThreshold);

        // Clamp pheromone matrix
        pheromone = blaze::clamp(pheromone, minPheromone, maxPheromone);

        // Print pheromone matrix
        // std::cout << "Pheromone matrix: " << std::endl;
        // std::cout << pheromone << std::endl;
    }

    for (int i = 0; i < bestCosts.size(); i++)
    {
        if (bestCosts.at(i) < bestCost)
        {
            bestCost = bestCosts.at(i);
            bestPath = bestPaths.at(i);
        }
    }

#ifdef TIME
    std::cout << "Total movements time: " << totalMovementsTime.count() << " seconds" << std::endl;
    std::cout << "Total local search time: " << totalLocalSearchTime.count() << " seconds" << std::endl;
#endif
}

int TSPStatement::localSearch2Opt(std::vector<blaze::DynamicVector<long unsigned int>> &paths, std::vector<double> &costs, int position)
{
    // 2-opt local search
    int isThereANewBest = 0; // 0 = no, 1 = yes

    blaze::DynamicVector<double, false> path = paths.at(position);

    /*blaze::DynamicVector<long unsigned int> testPath(8);
    testPath[0] = 0;
    testPath[1] = 1;
    testPath[2] = 2;
    testPath[3] = 3;
    testPath[4] = 4;
    testPath[5] = 5;
    testPath[6] = 6;
    testPath[7] = 7;
    int i = 1;
    int j = 5;
            long unsigned int a = i;
            long unsigned int b = i + 1;
            long unsigned int c = j;
            long unsigned int d = j + 1 % testPath.size();
    blaze::DynamicVector<long unsigned int> newPath(8);
                // Copy the first part of the path
                for (long unsigned int l = 0; l <= a; l++)
                {
                    newPath[l] = testPath[l];
                }
                // Copy the second part of the path, from c to b, reversed
                for (long unsigned int l = 0; l <= (c - b); l++)
                {
                    newPath[(a + 1) + l] = testPath[(c - l)];
                }
                // Copy the last part of the path
                for (long unsigned int l = d; l < testPath.size(); l++)
                {
                    newPath[l] = testPath[l];
                }
                    std::cout << "New path: " << newPath << std::endl;*/

    for (long unsigned int i = 0; i < path.size() - 3; i++)
    {
        for (long unsigned int j = i + 2; j < path.size() - 1; j++)
        {
            // Calculate the gain of the swap
            double gain = 0;
            long unsigned int a = i;
            long unsigned int b = i + 1;
            long unsigned int c = j;
            long unsigned int d = j + 1 % path.size();
            gain += distanceMatrix(path[a], path[c]);
            gain += distanceMatrix(path[b], path[d]);
            gain -= distanceMatrix(path[c], path[d]);
            gain -= distanceMatrix(path[a], path[b]);

            // If the gain is negative, perform the swap
            if (gain < 0)
            {
                // Perform the swap: the path becomes bestPath[a] -> bestPath[c] ... (reversed) bestPath[b] -> bestPath[d] ...
                blaze::DynamicVector<long unsigned int> newPath(path.size());
                // Copy the first part of the path
                for (long unsigned int l = 0; l <= a; l++)
                {
                    newPath[l] = path[l];
                }
                // Copy the second part of the path, from c to b, reversed
                for (long unsigned int l = 0; l <= (c - b); l++)
                {
                    newPath[(a + 1) + l] = path[(c - l)];
                }
                // Copy the last part of the path
                for (long unsigned int l = d; l < path.size(); l++)
                {
                    newPath[l] = path[l];
                }

                // Update the best cost
                costs.at(position) += gain;

                // Update the best path
                path = newPath;

                // std::cout << "Improved with 2-opt! Gain: " << gain << std::endl;
                isThereANewBest = 1;
            }
        }
    }

    paths.at(position) = path;

    return isThereANewBest;
};

int TSPStatement::localSearch3Opt(std::vector<blaze::DynamicVector<long unsigned int>> &paths, std::vector<double> &costs, int position)
{
    // 3-opt local search
    int isThereANewBest = 0; // 0 = no, 1 = yes
    blaze::DynamicVector<long unsigned int> path = paths.at(position);

    /*blaze::DynamicVector<long unsigned int> testPath(14);
    testPath[0] = 0;
    testPath[1] = 1;
    testPath[2] = 2;
    testPath[3] = 3;
    testPath[4] = 4;
    testPath[5] = 5;
    testPath[6] = 6;
    testPath[7] = 7;
    testPath[8] = 8;
    testPath[9] = 9;
    testPath[10] = 10;
    testPath[11] = 11;
    testPath[12] = 12;
    testPath[13] = 13;
    int i = 1;
    int j = 7;
    int k = 11;
    blaze::DynamicVector<long unsigned int> newPath(14);
                    for (long unsigned int l = 0; l <= i; l++) {
                        newPath[l] = testPath[l];
                    }
                    // Copy the second part of the path, from j + 1 to k
                    for (long unsigned int l = 0; l <= (k - (j + 1)); l++) {
                        newPath[(i + 1) + l] = testPath[(j + 1) + l];
                    }
                    // Copy the third part of the path, from i + 1 to j
                    for (long unsigned int l = 0; l <= (j - (i + 1)); l++) {
                        newPath[(i + 1) + (k - (j + 1)) + 1 + l] = testPath[(i + 1) + l];
                    }
                    // newPath[k + 1] = bestPath[k + 1];
                    for (long unsigned int l = k + 1; l < testPath.size(); l++) {
                        newPath[l] = testPath[l];
                    }
                    std::cout << "New path: " << newPath << std::endl;*/

    for (long unsigned int i = 0; i < path.size() - 3; i++)
    {
        for (long unsigned int j = i + 2; j < path.size() - 1; j++)
        {
            for (long unsigned int k = j + 2; k < path.size(); k++)
            {
                // Calculate the gain of the swap
                double gain = 0;
                gain += distanceMatrix(path[i], path[j + 1]);
                gain += distanceMatrix(path[k], path[i + 1]);
                gain += distanceMatrix(path[j], path[(k + 1) % path.size()]);
                gain -= distanceMatrix(path[i], path[i + 1]);
                gain -= distanceMatrix(path[j], path[j + 1]);
                gain -= distanceMatrix(path[k], path[(k + 1) % path.size()]);

                // If the gain is negative, perform the swap
                if (gain < 0)
                {
                    // Perform the swap: the path becomes bestPath[i] -> bestPath[j + 1] ... bestPath[k] -> bestPath[i+1] ... bestPath[j] -> bestPath[k+1] ...
                    blaze::DynamicVector<long unsigned int> newPath(path.size());
                    // Copy the first part of the path
                    for (long unsigned int l = 0; l <= i; l++)
                    {
                        newPath[l] = path[l];
                    }
                    // Copy the second part of the path, from j + 1 to k
                    for (long unsigned int l = 0; l <= (k - (j + 1)); l++)
                    {
                        newPath[(i + 1) + l] = path[(j + 1) + l];
                    }
                    // Copy the third part of the path, from i + 1 to j
                    for (long unsigned int l = 0; l <= (j - (i + 1)); l++)
                    {
                        newPath[(i + 1) + (k - (j + 1)) + 1 + l] = path[(i + 1) + l];
                    }
                    // newPath[k + 1] = path[k + 1];
                    for (long unsigned int l = k + 1; l < path.size(); l++)
                    {
                        newPath[l] = path[l];
                    }

                    // Update the best cost
                    costs.at(position) += gain;

                    // Update the best path
                    path = newPath;

                    // std::cout << "Improved with 3-opt! Gain: " << gain << std::endl;
                    isThereANewBest = 1;
                }
            }
        }
    }

    paths.at(position) = path;

    return isThereANewBest;
};

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

void TSPStatement::writeSolution(const char *filename)
{
    double costDiff = (std::floor(getBestCost()) - ((double)getBestKnown())) / ((double)getBestKnown());

    std::ofstream resultsFile;
    // The filename was "./problems/eil76.tsp", so we want to write the results to "./problems/eil76.opt.tour"
    std::string filenameString(filename);
    std::string filenameWithoutExtension = filenameString.substr(0, filenameString.find_last_of("."));
    std::string filenameWithExtension = filenameWithoutExtension + ".opt.tour";
    resultsFile.open(filenameWithExtension, std::ios_base::out); // if no output file is specified, use results.csv
    resultsFile << "NAME : " << getName() << std::endl;
    resultsFile << "COMMENT : " << comment << std::endl;
    resultsFile << "TYPE : TOUR" << std::endl;
    resultsFile << "DIMENSION : " << dimension << std::endl;
    resultsFile << "TOUR_LENGTH : " << getBestCost() << std::endl;
    resultsFile << "ERROR : " << costDiff << std::endl;
    resultsFile << "TOUR_SECTION" << std::endl;
    for (unsigned long int i = 0; i < getBestPath().size() - 1; i++)
    {
        resultsFile << getBestPath()[i] << std::endl;
    }
    resultsFile << "-1" << std::endl;
    resultsFile << "EOF" << std::endl;
    resultsFile.close();
}