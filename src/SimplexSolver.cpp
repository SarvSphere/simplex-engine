#include "SimplexSolver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

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
    VariableNames.clear();
    Yb.clear();
    Cb.clear();
    Status = SolutionStatus::NOT_SOLVED;
    IsPhase1 = false;

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
    ConstraintSigns.assign(NumConstraints,0);

    cout<<"Enter the coefficients of the objective function (Cj) : ";
    for(int j=0;j<NumVariables;j++){
        cin>>Cj[j];
        if(OptType == 2) Cj[j] *= -1;
    }

    cout<<"Enter the constraints (coefficients, then sign [1:<=, 2:>=, 3:=], then RHS) : \n";
    for(int i=0;i<NumConstraints;i++){
        for(int j=0;j<NumVariables;j++){
            cin>>Tableau[i][j];
        }

        cin>>ConstraintSigns[i];
        cin>>Tableau[i][NumVariables];

        if(Tableau[i][NumVariables] < 0){
            for(int j=0;j<=NumVariables;j++){
                Tableau[i][j] *= -1;
            }

            if(ConstraintSigns[i] == 1) ConstraintSigns[i] = 2;
            else if(ConstraintSigns[i] == 2) ConstraintSigns[i] = 1;

            cout<<"[System] Negative RHS in constraint "<<(i+1)<<". Multiplied by -1 and flipped sign.\n";
        }
    }

    for(int j=0;j<NumVariables;j++){
        VariableNames.push_back("x" + to_string(j+1));
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

        VariableNames.push_back("s" + to_string(i+1));
    }

    CalculateNetEvaluation();
}

void SimplexSolver::InitializeTwoPhase(){
    cout<<"\n[System] Initializing Two Phase Method (Phase 1)...\n";

    OriginalCj = Cj;

    for(int j=0;j<NumVariables;j++){
        Cj[j] = 0;
    }

    int sCount = 1;
    int aCount = 1;

    for(int i=0;i<NumConstraints;i++){
        if(ConstraintSigns[i] == 1){
            Cj.push_back(0.0);
            for(int j=0;j<=NumConstraints;j++){
                double val = (j==i) ? 1.0 : 0.0;
                Tableau[j].insert(Tableau[j].end() - 1,val);
            }
            Yb.push_back(Cj.size()-1);
            Cb.push_back(0.0);
            VariableNames.push_back("s" + to_string(sCount++));
        }

        else if(ConstraintSigns[i] == 2){
            Cj.push_back(0.0);
            for(int j=0;j<=NumConstraints;j++){
                double val = (j==i) ? -1.0 : 0.0;
                Tableau[j].insert(Tableau[j].end() - 1,val);
            }
            VariableNames.push_back("s" + to_string(sCount++));
            
            Cj.push_back(-1.0);
            for(int j=0;j<=NumConstraints;j++){
                double val = (j == i) ? 1.0 : 0.0;
                Tableau[j].insert(Tableau[j].end() - 1,val);
            }
            Yb.push_back(Cj.size()-1);
            Cb.push_back(-1.0);
            VariableNames.push_back("a" + to_string(aCount++));
        }

        else if(ConstraintSigns[i] == 3){ 
            Cj.push_back(-1.0);
            for(int j=0;j<=NumConstraints;j++){
                double val = (j == i) ? 1.0 : 0.0;
                Tableau[j].insert(Tableau[j].end() - 1,val);
            }
            Yb.push_back(Cj.size()-1);
            Cb.push_back(-1.0);
            VariableNames.push_back("a" + to_string(aCount++));
        }
    }

    CalculateNetEvaluation();
}

void SimplexSolver::SolveSimplex(){
    cout<<"\n[Engine] Running Standard Simplex Method...\n";

    int TotalCols = Tableau[0].size();

    const double EPSILON = 1e-7;

    while(true){
        int PivotCol = -1;
        double MostNegative = -EPSILON;

        for(int j=0;j<TotalCols-1;j++){
            if(Tableau[NumConstraints][j] < MostNegative){
                MostNegative = Tableau[NumConstraints][j];
                PivotCol = j;
            }
        }

        if(PivotCol == -1){
            if(IsPhase1) return;
            if(Status == SolutionStatus::ALTERNATE){
                break;
            }

            vector<bool>IsBasic(TotalCols-1,false);
            for(int i=0;i<NumConstraints;i++){
                IsBasic[Yb[i]] = true;
            }

            int AltCol = -1;
            bool InfiniteRay = false;
            bool TrueAlternateVertex = false;

            for(int j=0;j<TotalCols-1;j++){
                if(abs(Tableau[NumConstraints][j]) < EPSILON && !IsBasic[j]){
                    InfiniteRay = true;

                    double PreviewMinRatio = 1e9;
                    for(int i=0;i<NumConstraints;i++){
                        if(Tableau[i][j] > EPSILON){
                            double Ratio = Tableau[i][TotalCols-1] / Tableau[i][j];
                            if(Ratio < PreviewMinRatio){
                                PreviewMinRatio = Ratio;
                            }
                        }
                    }

                    if(PreviewMinRatio > EPSILON && PreviewMinRatio != 1e9){
                        AltCol = j;
                        TrueAlternateVertex = true;
                        InfiniteRay = false;
                        break;
                    }
                    else if(PreviewMinRatio <= EPSILON && AltCol == -1){
                        AltCol = j;
                        InfiniteRay = false;
                    }
                }
            }

            if(AltCol != -1 && TrueAlternateVertex){
                Status = SolutionStatus::ALTERNATE;
                cout<<"Result : Alternate (Infinite) Optimal Solutions exist.\n";
                cout<<"--- First Optimal Solution ---\n";

                double Z = Tableau[NumConstraints][TotalCols-1];
                if(OptType == 2) Z *= -1;
                if(abs(Z) < EPSILON) Z = 0.0;
                cout<<"Optimal Z = "<<Z<<"\n";

                vector<double>DecisionVariables(NumVariables,0.0);
                for(int i=0;i<NumConstraints;i++){
                    if(Yb[i] < NumVariables) DecisionVariables[Yb[i]] = Tableau[i][TotalCols-1];
                }

                for(int j=0;j<NumVariables;j++){
                    double val = DecisionVariables[j];
                    if(abs(val) < EPSILON) val = 0.0;
                    cout<<"x"<<(j+1)<<" = "<<val<<"\n";
                }

                cout<<"\n[Engine] Pivoting for second basic feasible solution\n";

                PivotCol = AltCol;
            }

            else if(InfiniteRay){
                Status = SolutionStatus::ALTERNATE_RAY;
                break;
            }

            else{
                Status = SolutionStatus::UNIQUE;
                if(AltCol != -1 && !TrueAlternateVertex){
                    cout<<"\n[System] Note : Degenerate alternate basis detected, but it will give same vertex (Unique Solution).\n";
                }
                break;
            }
        }

        int PivotRow = -1;
        double MinRatio = 1e9;

        for(int i=0;i<NumConstraints;i++){
            if(Tableau[i][PivotCol] > EPSILON){
                double Ratio = Tableau[i][TotalCols-1] / Tableau[i][PivotCol];
                if(Ratio < MinRatio){
                    MinRatio = Ratio;
                    PivotRow = i;
                }
            }
        }

        if(PivotRow == -1){
            if(IsPhase1) return;
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

        for(int i=0;i<=NumConstraints;i++){
            for(int j=0;j<TotalCols;j++){
                if(abs(Tableau[i][j]) < EPSILON){
                    Tableau[i][j] = 0.0;
                }
            }
        }

        PrintTableau();
    }
}

void SimplexSolver::SolveTwoPhase(){
    bool NeedsPhase1 = false;
    for(int i=0;i<NumConstraints;i++){
        if(ConstraintSigns[i] == 2 || ConstraintSigns[i] == 3){
            NeedsPhase1 = true;
            break;
        }
    }

    if(!NeedsPhase1){
        cout<<"\n[System] Given problem can be solved by Standard Simplex Method. No need of Two Phase Simplex Method\n";

        AddSlackVariables();
        PrintTableau();
        SolveSimplex();

        return;
    }

    cout<<"[Engine] Running Two-Phase Simplex Method...\n";

    IsPhase1 = true;
    InitializeTwoPhase();
    PrintTableau();

    SolveSimplex();

    IsPhase1 = false;

    int TotalCols = Tableau[0].size();
    double Phase1Z = Tableau[NumConstraints][TotalCols-1];

    if(Phase1Z < -1e-5){
        cout<<"\n[System] Phase 1 is completed. Artificial variables can not be eliminated.\n";
        Status = SolutionStatus::INFEASIBLE;
        return;
    }

    cout<<"\n[System] Phase 1 is completed. Transitioning to Phase 2...\n";

    for(int j=Cj.size()-1;j>=NumVariables;j--){
        if(Cj[j] == -1.0){
            bool IsBasic = false;

            int BasicRow = -1;

            for(int i=0;i<NumConstraints;i++){
                if(Yb[i] == j){
                    IsBasic = true;
                    BasicRow = i;
                    break;
                }
            }

            if(!IsBasic){
                Cj.erase(Cj.begin()+j);
                VariableNames.erase(VariableNames.begin() + j);
                for(int i=0;i<=NumConstraints;i++){
                    Tableau[i].erase(Tableau[i].begin()+j);
                }

                for(int i=0;i<NumConstraints;i++){
                    if(Yb[i] > j) Yb[i]--;
                }
            }
            else{
                int PivotCol = -1;
                for(int col=0;col<Cj.size();col++){
                    if(abs(Cj[col] - (-1.0)) > 1e-7 && abs(Tableau[BasicRow][col]) > 1e-7){
                        PivotCol = col;
                        break;
                    }
                }

                if(PivotCol != -1){
                    cout<<"[System] Degenerate basic artificial variable found at column "<<(j+1)<<" in row "<<(BasicRow+1)<<".\n";
                    cout<<"[System] Pivoting non-artificial variable "<<VariableNames[PivotCol]<<" in.\n";

                    Yb[BasicRow] = PivotCol;
                    Cb[BasicRow] = Cj[PivotCol]; 

                    double PivotElement = Tableau[BasicRow][PivotCol];

                    for(int col=0;col<Cj.size()+1;col++){
                        Tableau[BasicRow][col] /= PivotElement;
                    }

                    for(int i=0;i<=NumConstraints;i++){
                        if(i != BasicRow){
                            double factor = Tableau[i][PivotCol];
                            for(int col=0;col<Cj.size()+1;col++){
                                Tableau[i][col] -= factor*Tableau[BasicRow][col];
                            }
                        }
                    }

                    for(int i=0;i<=NumConstraints;i++){
                        Tableau[i].erase(Tableau[i].begin() + j);
                    }
                    Cj.erase(Cj.begin() + j);
                    VariableNames.erase(VariableNames.begin() + j);

                    for(int i=0;i<NumConstraints;i++){
                        if(Yb[i] > j) Yb[i]--;
                    }
                }

                else{
                    cout<<"[System] Redundant constraint detected (Degenerate Artificial Variable at column "<<(j+1)<<"). Deleting row "<<(BasicRow+1)<<" from the tableau.\n";
                    
                    Tableau.erase(Tableau.begin()+BasicRow);
                    Yb.erase(Yb.begin()+BasicRow);
                    Cb.erase(Cb.begin()+BasicRow);

                    NumConstraints--;

                    for(int i=0;i<=NumConstraints;i++){
                        Tableau[i].erase(Tableau[i].begin()+j);
                    }

                    Cj.erase(Cj.begin()+j);
                    VariableNames.erase(VariableNames.begin()+j);

                    for(int i=0;i<NumConstraints;i++){
                        if(Yb[i] > j) Yb[i]--;
                    }
                }
            }
        }
    }

    for(int j=Cj.size()-1;j>=NumVariables;j--){
        bool IsPhantom = true;
        for(int i=0;i<NumConstraints;i++){
            if(abs(Tableau[i][j]) > 1e-7){
                IsPhantom = false;
                break;
            }
        }

        if(IsPhantom){
            cout<<"[System] Phantom variable detected ("<<VariableNames[j]<<") from deleted row. Removing this column.\n";
            Cj.erase(Cj.begin()+j);
            VariableNames.erase(VariableNames.begin()+j);
            for(int i=0;i<=NumConstraints;i++){
                Tableau[i].erase(Tableau[i].begin()+j);
            }
            for(int i=0;i<NumConstraints;i++){
                if(Yb[i] > j) Yb[i]--;
            }
        }
    }

    for(int j=0;j<NumVariables;j++){
        Cj[j] = OriginalCj[j];
    }

    for(int i=0;i<NumConstraints;i++){
        Cb[i] = Cj[Yb[i]];
    }

    CalculateNetEvaluation();

    cout<<"[System] Phase 2 Initialised. Running Phase 2 calculations\n";
    Status = SolutionStatus::NOT_SOLVED;
    PrintTableau();

    SolveSimplex();
}

void SimplexSolver::PrintTableau(){
    cout<<"\nPrinting current Simplex tableau...\n\n";

    int W = 18;
    int CbW = 18;
    int YbW = 6;

    cout<<setw(CbW + YbW)<<"Cj";
    for(int i=0;i<Cj.size();i++){
        cout<<setw(W)<<Cj[i];
    }
    cout<<"\n";

    cout<<setw(CbW)<<"Cb"<<setw(YbW)<<"Yb";
    for(int i=0;i<Cj.size();i++){
        cout<<setw(W)<<VariableNames[i];
    }

    cout<<setw(W)<<"Xb\n";

    int LineLen = CbW + YbW + (Cj.size()+1)*W;
    cout<<string(LineLen,'-')<<"\n";

    for(int i=0;i<NumConstraints;i++){
        cout<<setw(CbW)<<Cb[i]<<setw(YbW)<<VariableNames[Yb[i]];
        for(int j=0;j<Tableau[0].size();j++){
            if(j == Tableau[0].size() - 1){
                cout<<" |"<<setw(W-2)<<Tableau[i][j];
            }
            else{
                cout<<setw(W)<<Tableau[i][j];
            }
        }
        cout<<"\n";
    }

    cout<<string(LineLen,'-')<<"\n";

    cout<<setw(CbW + YbW)<<"Zj-Cj";

    for(int i=0;i<Tableau[0].size()-1;i++){
        cout<<setw(W)<<Tableau[NumConstraints][i];
    }

    cout<<" |      Z = "<<Tableau[NumConstraints].back()<<"\n\n";
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