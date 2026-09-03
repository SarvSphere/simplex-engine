#include <iostream>
#include "SimplexSolver.h"

using namespace std;

int main(){
    cout<<"--- Linear Programming Engine ---\n";
    
    SimplexSolver engine;

    engine.LoadEquations();
    engine.AddSlackVariables();
    engine.PrintTableau();
    engine.SolveSimplex();
    engine.PrintResult();
    
    return 0;
}