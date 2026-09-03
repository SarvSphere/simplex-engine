#include "SimplexSolver.h"
#include <iostream>
#include <iomanip>

using namespace std;

SimplexSolver::SimplexSolver(){
    Status = SolutionStatus::NOT_SOLVED;
    NumVariables = 0;
    NumConstraints = 0;

    cout<<fixed<<setprecision(2);
}

void SimplexSolver::CalculateNetEvaluation(){
    int totcols = Tableau[0].size();

    for(int i=0;i<totcols-1;i++){
        double Zj = 0.0;
        for(int j=0;j<NumConstraints;j++){
            Zj += Cb[j] * Tableau[j][i];
        }
        Tableau[NumConstraints][i] = Zj - Cj[i];
    }

    double Z = 0.0;
    for(int i=0;i<NumConstraints;i++){
        Z += Cb[i] * Tableau[i][totcols-1];
    }
    Tableau[NumConstraints][totcols-1] = Z;
}

void SimplexSolver::LoadEquations(){
    cout<<"Enter 1 for Maximization Problem\n";
    cout<<"Enter 2 for Minimization Problem\n";
    cout<<"Enter 1 or 2 : ";
    cin>>OptType;

    cout<<"Enter number of decision variables : ";
    cin>>NumVariables;

    cout<<"Enter number of constraints : ";
    cin>>NumConstraints;

    Tableau.assign(NumConstraints + 1,vector<double>(NumVariables + 1,0.0));
    Cj.assign(NumVariables,0.0);

    cout<<"Enter the coefficients of the objective function (Cj) : ";
    for(int j=0;j<NumVariables;j++){
        cin>>Cj[j];
        if(OptType == 2) Cj[j] *= -1;
    }

    cout<<"Enter the constraints (coefficients followed by their RHS value) : \n";
    for(int i=0;i<NumConstraints;i++){
        for(int j=0;j<NumVariables;j++){
            cin>>Tableau[i][j];
        }
        cin>>Tableau[i][NumVariables];
    }

    cout<<"\n[System] Equations loaded successfully.\n";
}

void SimplexSolver::AddSlackVariables(){
    cout<<"[System] Adding slack variables to form initial basis...\n";

    for(int i=0;i<NumConstraints;i++){
        Cj.push_back(0.0);

        for(int j=0;j<NumConstraints;j++){
            double SlackValue = (j == i) ? 1.0 : 0.0;
            Tableau[j].insert(Tableau[j].end() - 1,SlackValue); 
        }

        Tableau[NumConstraints].insert(Tableau[NumConstraints].end() - 1,0.0);

        Yb.push_back(NumVariables + i);
        Cb.push_back(0.0);
    }

    CalculateNetEvaluation();
}

void SimplexSolver::SolveSimplex(){
    cout<<"\n[Engine] Running Standard Simplex Method...\n";

    int TotalCols = Tableau[0].size();

    while(true){
        int PivotCol = -1;
        double MostNegative = 0.0;

        for(int j=0;j<TotalCols-1;j++){
            if(Tableau[NumConstraints][j] < MostNegative){
                MostNegative = Tableau[NumConstraints][j];
                PivotCol = j;
            }
        }

        if(PivotCol == -1){
            if(Status == SolutionStatus::ALTERNATE){
                break;
            }

            vector<bool>IsBasic(TotalCols-1,false);
            for(int i=0;i<NumConstraints;i++){
                IsBasic[Yb[i]] = true;
            }

            int AltCol = -1;
            bool InfiniteSolution = false;
            for(int j=0;j<TotalCols-1;j++){
                if(Tableau[NumConstraints][j] == 0 && !IsBasic[j]){
                    InfiniteSolution = true;
                    for(int i=0;i<NumConstraints;i++){
                        if(Tableau[i][j] > 0){
                            AltCol = j;
                            break;
                        }
                    }
                }
                if(AltCol != -1) break;
            }

            if(AltCol != -1){
                Status = SolutionStatus::ALTERNATE;
                cout<<"Result : Alternate (Infinite) Optimal Solutions exist.\n";
                cout<<"--- First Optimal Solution ---\n";

                double Z = Tableau[NumConstraints][TotalCols-1];
                if(OptType == 2) Z *= -1;
                cout<<"Optimal Z = "<<Z<<"\n";

                vector<double>DecisionVariables(NumVariables,0.0);
                for(int i=0;i<NumConstraints;i++){
                    if(Yb[i] < NumVariables) DecisionVariables[Yb[i]] = Tableau[i][TotalCols-1];
                }

                for(int j=0;j<NumVariables;j++){
                    cout<<"x"<<(j+1)<<" = "<<DecisionVariables[j]<<"\n";
                }

                cout<<"\n[Engine] Pivoting for second basic feasible solution\n";

                PivotCol = AltCol;
            }

            else if(InfiniteSolution){
                Status = SolutionStatus::ALTERNATE_RAY;
                break;
            }

            else{
                Status = SolutionStatus::UNIQUE;
                break;
            }
        }

        int PivotRow = -1;
        double MinRatio = 1e9;

        for(int i=0;i<NumConstraints;i++){
            if(Tableau[i][PivotCol] > 0){
                double Ratio = Tableau[i][TotalCols-1] / Tableau[i][PivotCol];
                if(Ratio < MinRatio){
                    MinRatio = Ratio;
                    PivotRow = i;
                }
            }
        }

        if(PivotRow == -1){
            Status = SolutionStatus::UNBOUNDED;
            break;
        }

        cout<<"\n[Engine] Pivoting on Row "<<(PivotRow+1)<<" & Column "<<(PivotCol+1)<<"\n";

        Yb[PivotRow] = PivotCol;
        Cb[PivotRow] = Cj[PivotCol];

        double PivotElement = Tableau[PivotRow][PivotCol];

        for(int j=0;j<TotalCols;j++){
            Tableau[PivotRow][j] /= PivotElement;
        }

        for(int i=0;i<=NumConstraints;i++){
            if(i != PivotRow){
                double Factor = Tableau[i][PivotCol];
                for(int j=0;j<TotalCols;j++){
                    Tableau[i][j] -= Factor * Tableau[PivotRow][j];
                }
            }
        }

        PrintTableau();
    }
}

void SimplexSolver::SolveTwoPhase(){
    cout<<"[Engine] Running Two-Phase Simplex Method...\n";
}

void SimplexSolver::PrintTableau(){
    cout<<"\nPrinting current Simplex tableau...\n\n";

    cout<<"\tCj\t";
    for(int i=0;i<Cj.size();i++){
        cout<<Cj[i]<<"\t";
    }
    cout<<"\n";

    cout<<"Cb\tYb\t";
    for(int i=0;i<Cj.size();i++){
        cout<<"y"<<(i+1)<<"\t";
    }

    cout<<"Xb\n";
    cout<<"--------------------------------------------------------\n";

    for(int i=0;i<NumConstraints;i++){
        cout<<Cb[i]<<"\ty"<<(Yb[i]+1)<<"\t";
        for(int j=0;j<Tableau[0].size();j++){
            if(j == Tableau[0].size() - 1) cout<<"| ";
            cout<<Tableau[i][j]<<"\t";
        }
        cout<<"\n";
    }

    cout<<"--------------------------------------------------------\n";

    cout<<"Zj-Cj\t\t";
    for(int i=0;i<Tableau[0].size()-1;i++){
        cout<<Tableau[NumConstraints][i]<<"\t";
    }

    cout<<"| Z = "<<Tableau[NumConstraints].back()<<"\n\n";
}

void SimplexSolver::PrintResult(){
    switch(Status){
        case SolutionStatus::NOT_SOLVED:
            cout<<"\n--- Final Conclusion ---\n";
            cout<<"Result : Problem has not been solved yet.\n";
            break;

        case SolutionStatus::UNIQUE:
        {
            cout<<"\n--- Final Conclusion ---\n";
            cout<<"Result : Unique Optimal Solution found.\n";

            int TotalColumns = Tableau[0].size();
            double Z = Tableau[NumConstraints][TotalColumns-1];
            if(OptType == 2) Z *= -1;

            cout<<"Optimal Z = "<<Z<<"\n";
            cout<<"Decision Variables :\n";

            vector<double>DecisionVariables(NumVariables,0.0);

            for(int i=0;i<NumConstraints;i++){
                if(Yb[i] < NumVariables){
                    DecisionVariables[Yb[i]] = Tableau[i][TotalColumns-1];
                }
            }

            for(int j=0;j<NumVariables;j++){
                cout<<"x"<<(j+1)<<" = "<<DecisionVariables[j]<<"\n";
            }

            break;
        }

        case SolutionStatus::ALTERNATE:{
            cout<<"--- Second Optimal Solution ---\n";

            int TotalCols = Tableau[0].size();
            double Z = Tableau[NumConstraints][TotalCols-1];
            if(OptType == 2) Z *= -1;

            cout<<"Optimal Z = "<<Z<<"\n";
            cout<<"Decision Variables :\n";

            vector<double>DecisionVariables(NumVariables,0.0);
            for(int i=0;i<NumConstraints;i++){
                if(Yb[i] < NumVariables){
                    DecisionVariables[Yb[i]] = Tableau[i][TotalCols-1];
                }
            }

            for(int j=0;j<NumVariables;j++){
                cout<<"x"<<(j+1)<<" = "<<DecisionVariables[j]<<"\n";
            }

            break;
        }

        case SolutionStatus::ALTERNATE_RAY:
        {
            cout<<"\n--- Final Conclusion ---\n";
            cout<<"Result : Alternate (Infinite) Optimal Solutions exist but only one Optimal Basic Feasible Solution exists.\n";

            int TotalCols = Tableau[0].size();
            double Z = Tableau[NumConstraints][TotalCols-1];
            if(OptType == 2) Z *= -1;

            cout<<"Optimal Z = "<<Z<<"\n";
            cout<<"Decision Variables :\n";

            vector<double>DecisionVariables(NumVariables,0.0);
            for(int i=0;i<NumConstraints;i++){
                if(Yb[i] < NumVariables){
                    DecisionVariables[Yb[i]] = Tableau[i][TotalCols-1];
                }
            }

            for(int j=0;j<NumVariables;j++){
                cout<<"x"<<(j+1)<<" = "<<DecisionVariables[j]<<"\n";
            }

            break;
        }

        case SolutionStatus::UNBOUNDED:
            cout<<"\n--- Final Conclusion ---\n";
            cout<<"Result : Solution is Unbounded.\n";
            break;

        case SolutionStatus::INFEASIBLE:
            cout<<"\n--- Final Conclusion ---\n";
            cout<<"Result : Problem is Infeasible (No Solution).\n";
            break;
    }
}