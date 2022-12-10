#include "../headers/tsp.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

int main(int argc, char** argv) {
    TSPStatement statement;

    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    std::cout << "Reading file: " << argv[1] << std::endl;

    statement.read(argv[1]);

    srand(1); // Set the seed for the random number generator. This is so that the results are reproducible.

    statement.solve_aco();

    return 0;
}