#include "../headers/tsp.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

double TSPConstants::alpha = 1;
double TSPConstants::beta = 1;
double TSPConstants::rho = 0.05; // Evaporation rate
double TSPConstants::phi = 0.05;  // Pheromone decay rate
double TSPConstants::tau0 = 0.0002;  // Initial pheromone (NOT USED)
double TSPConstants::maxPheromone = 0.1;           // max Pheromone value (NOT USED)
double TSPConstants::minPheromone = 10e-6;           // max Pheromone value
long unsigned int TSPConstants::nAnts = 15;        // Number of ants

int main(int argc, char **argv)
{
    TSPStatement statement;
    double factor_tau = 1.0;

    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <filename> [<out_filename>] [<alpha>] [<beta>] [rho] [phi] [tau0_factor] [num_ants] [min_pheromone]" << std::endl;
        return 1;
    }

    if (argc > 3)
    {
        std::cout << "Using alpha = " << argv[3] << " beta = " << argv[4] << " rho = " << argv[5] << " phi = " << argv[6] << " tau0 = " << argv[7] << " nAnts = " << argv[8] << " min_pheromone = " << argv[9] << std::endl;
        TSPConstants::alpha = std::stod(argv[3]);
        TSPConstants::beta = std::stod(argv[4]);
        TSPConstants::rho = std::stod(argv[5]);
        TSPConstants::phi = std::stod(argv[6]);
        factor_tau = std::stod(argv[7]);
        TSPConstants::nAnts = std::stoi(argv[8]);
        TSPConstants::minPheromone = std::stod(argv[9]);
    }

    std::cout << "Reading file: " << argv[1] << std::endl;

    statement.read(argv[1]);

    TSPConstants::nAnts = TSPConstants::nAnts * ceil(((double)statement.getDimension()) / 100.0);
    TSPConstants::maxPheromone = 1.0 / (TSPConstants::rho * ((double)statement.getBestKnown()));
    TSPConstants::tau0 = pow(TSPConstants::maxPheromone, factor_tau);

    srand(1); // Set the seed for the random number generator. This is so that the results are reproducible.

    statement.solve_aco();

    std::cout << "Best path found: " << statement.getBestCost() << std::endl;
    double costDiff = (std::floor(statement.getBestCost()) - ((double)statement.getBestKnown())) / ((double)statement.getBestKnown());
    std::cout << "Cost difference: " << costDiff << std::endl;

    // Append the result to the results file
    statement.writeSolution((argc > 2 ? argv[2] : argv[1])); // If the output filename is not specified, use the input filename

    return 0;
}