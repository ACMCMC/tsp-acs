#include "../headers/tsp.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

int main(int argc, char **argv)
{
    TSPStatement statement;

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::cout << "Reading file: " << argv[1] << std::endl;

    statement.read(argv[1]);

    srand(1); // Set the seed for the random number generator. This is so that the results are reproducible.

    statement.solve_aco();

    std::cout << "Best path found: " << statement.getBestPath() << std::endl;
    std::cout << "Cost difference: " << (statement.getBestCost() - ((double) statement.getBestKnown())) / ((double) statement.getBestKnown()) << std::endl;

    // Append the result to the results file
    std::ofstream resultsFile;
    resultsFile.open("results.csv", std::ios_base::app);
    resultsFile << statement.getName() << "," << std::floor(statement.getBestCost()) << std::endl;
    resultsFile.close();

    return 0;
}