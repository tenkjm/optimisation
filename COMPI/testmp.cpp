#include <iostream>
#include "mpproblem.hpp"
#include "mputils.hpp"
#include "funccnt.hpp"

class F : public COMPI::Functor <double> {
    public:
        
        double func(const double* x) {
            return x[0] * x[1];
        }
}; 

int main() {
    COMPI::MPProblem<double> mpp;
    mpp.mObjectives.push_back(new COMPI::FuncCnt<double>(* new F()));
    const int n = 4, neq = 2, nineq = 3;
    for(int i = 0; i < n; i ++) {
        int v = COMPI::MPProblem<double>::VariableTypes::BOOLEAN;
        mpp.mVarTypes.push_back(v);
    }
    
    for(int i = 0; i < neq; i ++) {
        mpp.mEqConstr.push_back(new F());
    }
    for(int i = 0; i < nineq; i ++) {
        mpp.mEqConstr.push_back(new F());
    }
    
    snowgoose::Box<double> box(n);
    for(int i = 0; i < n; i ++) {
        box.mA[i] = 0;
        box.mB[i] = 1;
    }
    mpp.mBox = &box;
    SG_ASSERT(COMPI::MPUtils::getProblemType(mpp) == COMPI::MPUtils::ProblemTypes::BOXCONSTR | 
                                                     COMPI::MPUtils::ProblemTypes::SINGLEOBJ |
                                                     COMPI::MPUtils::ProblemTypes::CONTINUOUS);
    
    
    return 0;
}
