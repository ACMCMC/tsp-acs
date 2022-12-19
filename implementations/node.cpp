#include "../headers/tsp.hpp"

#include <iostream>

Node::Node(long unsigned int id, double x, double y)
{
    this->id = id - 1;
    this->x = x;
    this->y = y;
}
double Node::getX() const
{
    return x;
}
double Node::getY() const
{
    return y;
}
long unsigned int Node::getInternalId() const
{
    return id;
}
long unsigned int Node::getHumanId() const
{
    return id + 1;
}

/**
 * @brief Node destructor
 * 
 */
Node::~Node()
{
    //std::cout << "Node " << id << " destroyed" << std::endl;
}