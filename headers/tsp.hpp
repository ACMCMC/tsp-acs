#include <array>
#include <string>
#include <blaze/Blaze.h>

/*class TSPSolution {
    public:
        TSPSolution();
        TSPSolution(const TSPSolution& orig);
        virtual ~TSPSolution();
        void setSolution(const std::array<int>& solution);
        void setCost(double cost);
        double getCost() const;
        std::array<int> getSolution() const;
        void printSolution() const;
    private:
        std::array<int> solution;
        double cost;
};*/

class TSPStatement {
    public:
        void read(const char* filename);
        int getDimension() const;
        double getDistance(int i, int j) const;
    private:
        int dimension;
        std::string name;
        int bestKnown;
        blaze::DynamicMatrix<double> *distanceMatrix;
};