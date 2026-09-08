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
    std::vector<std::vector<double>> Tableau;

    std::vector<double>Cj;
    std::vector<double>Cb;
    std::vector<int>Yb;
    std::vector<int>ConstraintSigns;
    std::vector<double>OriginalCj;
    std::vector<std::string>VariableNames;

    bool IsPhase1;

    SolutionStatus Status;

    int NumVariables;
    int NumConstraints;
    int OptType;

    void CalculateNetEvaluation();

public:
    SimplexSolver();

    void LoadEquations();
    void AddSlackVariables();
    void InitializeTwoPhase();

    void SolveSimplex();
    void SolveTwoPhase();

    void PrintTableau();
    void PrintResult();
};

#endif