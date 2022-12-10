#include "../headers/tsp.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

double TSPConstants::alpha = 1;
double TSPConstants::beta = 5;

int main(int argc, char **argv)
{
    TSPStatement statement;

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <filename> [<out_filename>] [<alpha>] [<beta>]" << std::endl;
        return 1;
    }

    if (argc > 3)
    {
        std::cout << "Using alpha = " << argv[3] << " and beta = " << argv[4] << std::endl;
        TSPConstants::alpha = std::stod(argv[3]);
        TSPConstants::beta = std::stod(argv[4]);
    } else {
    }

    std::cout << "Reading file: " << argv[1] << std::endl;

    statement.read(argv[1]);

    srand(1); // Set the seed for the random number generator. This is so that the results are reproducible.

    statement.solve_aco();

    std::cout << "Best path found: " << statement.getBestPath() << std::endl;
    std::cout << "Cost difference: " << (statement.getBestCost() - ((double) statement.getBestKnown())) / ((double) statement.getBestKnown()) << std::endl;

    // Append the result to the results file
    std::ofstream resultsFile;
    resultsFile.open((argc > 2 ? argv[2] : "results.csv"), std::ios_base::app); // if no output file is specified, use results.csv
    resultsFile << statement.getName() << "," << std::floor(statement.getBestCost()) << std::endl;
    resultsFile.close();

    return 0;
}