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

TSPStatement::Node::Node(int id, double x, double y)
{
    this->id = id;
    this->x = x;
    this->y = y;
}
double TSPStatement::Node::getX() const
{
    return x;
}
double TSPStatement::Node::getY() const
{
    return y;
}
int TSPStatement::Node::getId() const
{
    return id;
}

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
