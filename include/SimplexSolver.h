#ifndef SIMPLEX_SOLVER_H
#define SIMPLEX_SOLVER_H

#include <vector>
#include <string>

enum class SolutionStatus{
    NOT_SOLVED,
    UNIQUE,
    ALTERNATE,
    ALTERNATE_RAY,
    UNBOUNDED,
    INFEASIBLE
};

class SimplexSolver{
private:
    // A 2D matrix to hold constraints
    std::vector<std::vector<double>> Tableau;

    std::vector<double>Cj; // Objective function coefficients
    std::vector<double>Cb; // Cost of current basic variables
    std::vector<int>Yb; // Indices of current basic variable

    // Status of our solution
    SolutionStatus Status;

    // Number of variables & Number of constraints respectively
    int NumVariables;
    int NumConstraints;
    int OptType; // 1 for Maximization, 2 for Minimization

    // Function to calculate net evaluation
    void CalculateNetEvaluation();

public:
    SimplexSolver();

    // For passing equations to the engine
    void LoadEquations();
    void AddSlackVariables();

    // Core methods
    void SolveSimplex();
    void SolveTwoPhase();

    // Helper tools
    void PrintTableau();
    void PrintResult();
};

#endif