#include <iostream>
#include "SimplexSolver.h"

using namespace std;

int main(){
    cout<<"--- Linear Programming Engine ---\n";
    
    SimplexSolver engine;

    engine.LoadEquations();
    engine.SolveTwoPhase();
    engine.PrintResult();
    
    return 0;
}