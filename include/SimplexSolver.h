#ifndef SIMPLEX_SOLVER_H
#define SIMPLEX_SOLVER_H

#include <vector>

class SimplexSolver{
private:
    // A 2D matrix to hold the objective function and constraints
    std::vector<std::vector<double>> tableau;

public:
    SimplexSolver(); 
    void Solve();
    void PrintTableau();
};

#endif