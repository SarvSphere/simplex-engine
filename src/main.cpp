#include <iostream>
#include "SimplexSolver.h"

using namespace std;

int main(){
    cout<<"--- Linear Programming Engine ---\n";
    
    SimplexSolver engine;
    engine.Solve();

    return 0;
}