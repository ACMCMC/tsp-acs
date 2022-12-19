#include "../headers/tsp.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

double TSPConstants::alpha = 3.854172313492491;
double TSPConstants::beta = 9.502065547951132;
double TSPConstants::rho = 0.00018384087804787033;      // Evaporation rate
double TSPConstants::phi = 0.0068969337927708005;      // Pheromone decay rate
double TSPConstants::tau0 = 0.023001192342792635;       // Initial pheromone
long unsigned int TSPConstants::nAnts = 10;      // Number of ants

int main(int argc, char **argv)
{
    TSPStatement statement;

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <filename> [<out_filename>] [<alpha>] [<beta>] [rho] [phi] [tau0] [num_ants]" << std::endl;
        return 1;
    }

    if (argc > 3)
    {
        std::cout << "Using alpha = " << argv[3] << " beta = " << argv[4] << " rho = " << argv[5] << " phi = " << argv[6] << " tau0 = " << argv[7] << " nAnts = " << argv[8] << std::endl;
        TSPConstants::alpha = std::stod(argv[3]);
        TSPConstants::beta = std::stod(argv[4]);
        TSPConstants::rho = std::stod(argv[5]);
        TSPConstants::phi = std::stod(argv[6]);
        TSPConstants::tau0 = std::stod(argv[7]);
        TSPConstants::nAnts = std::stoi(argv[8]);
    }

    std::cout << "Reading file: " << argv[1] << std::endl;

    statement.read(argv[1]);

    TSPConstants::nAnts = TSPConstants::nAnts * ceil(((double) statement.getDimension()) / 100.0);

    srand(1); // Set the seed for the random number generator. This is so that the results are reproducible.

    statement.solve_aco();

    std::cout << "Best path found: " << statement.getBestPath() << std::endl;
    double costDiff = (std::floor(statement.getBestCost()) - ((double) statement.getBestKnown())) / ((double) statement.getBestKnown());
    std::cout << "Cost difference: " << costDiff << std::endl;

    // Append the result to the results file
    statement.writeSolution((argc > 2 ? argv[2] : argv[1])); // If the output filename is not specified, use the input filename

    return 0;
}