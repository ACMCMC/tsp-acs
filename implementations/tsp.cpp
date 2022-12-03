#include "../headers/tsp.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

/*
NAME: ch130
TYPE: TSP
COMMENT: 130 city problem (Churritz)
DIMENSION: 130
EDGE_WEIGHT_TYPE: EUC_2D
BEST_KNOWN : 6110
NODE_COORD_SECTION
1 334.5909245845 161.7809319139
*/

void TSPStatement::read(const char* filename) {
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("NAME") != std::string::npos) {
            std::istringstream iss(line);
            std::string word;
            iss >> word >> word >> name;
        }
        if (line.find("BEST_KNOWN") != std::string::npos) {
            std::istringstream iss(line);
            std::string word;
            iss >> word >> word >> bestKnown;
        }
        if (line.find("DIMENSION") != std::string::npos) {
            std::istringstream iss(line);
            std::string word;
            iss >> word >> word >> dimension;
            distance.resize(dimension);
            distance.fill({});
        }
        if (line.find("EDGE_WEIGHT_SECTION") != std::string::npos) {
            for (int i = 0; i < dimension; i++) {
                distance[i].resize(dimension);
                for (int j = 0; j < dimension; j++) {
                    file >> distance[i][j];
                }
            }
        }
    }
}