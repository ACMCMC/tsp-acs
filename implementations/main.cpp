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
    //std::cout << statement.getDimension() << std::endl;
    //std::cout << statement.getDistance(0, 1) << std::endl;
    return 0;
}