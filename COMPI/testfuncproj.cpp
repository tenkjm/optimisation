/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream>
#include <common/sgerrcheck.hpp>
#include "funcproj.hpp"

class Dejong : public COMPI::Functor<double> {
    
    double func(const double* x) {
        return x[0] * x[0] + x[1] * x[1];
    }
    
};
int main() {
    
    Dejong dejong;
    double const dir[] = {1, 0};
    double const xinit[] = {0, 1};
    double const t = 1;
    COMPI::FunctorProjector<double> dejongProj(dejong, 2, xinit, dir);
    std::cout << "Function value = " << dejongProj.func(&t) << "\n";
    SG_ASSERT(dejongProj.func(&t) == 2);
    return 0;
}
