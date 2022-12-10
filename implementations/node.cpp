#include "../headers/tsp.hpp"

#include <iostream>

Node::Node(int id, double x, double y)
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
int Node::getInternalId() const
{
    return id;
}
int Node::getHumanId() const
{
    return id + 1;
}

/**
 * @brief Node destructor
 * 
 */
Node::~Node()
{
    std::cout << "Node " << id << " destroyed" << std::endl;
}