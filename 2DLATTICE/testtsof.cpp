#include <iostream>
#include "ppenergy.hpp"
#include "atoms.hpp"
#include "tsofenergy.hpp"
#include "carbontersoff.hpp"

int main(int argc, char* argv[]) {

    lattice::LatticeData data;

    /**
     * Setup lattice data
     */
    data.mNumLayers = 4;
    data.mLength = 16;
    data.mRadius = 3;
    data.mLayersAtoms[0] = lattice::CARBON;
    data.mLayersAtoms[1] = lattice::CARBON;
    data.mLayersAtoms[2] = lattice::CARBON;
    data.mLayersAtoms[3] = lattice::CARBON;

    /**
     * Setup potential
     */
    lattice::TersoffParams tparam;
    lattice::fillCarbonParametersTersoffOriginal(tparam);
    lattice::TersoffUtils tutils(tparam);
    lattice::LatticeUtils lut(data);
    lattice::TersoffEnergy enrg(lut, tutils);

    /**
     * Setup lattice parameters
     */
    double x[] = {1.46, 0, 2.5, 0.73, 1.25, 2.5, 1.46, 1.25, 2.5, 0.73, 0, 2.5};

    double v = enrg.energy(x);

    std::cout << "v = " << v << "\n";

    for(auto i : enrg.getLBounds()) {
        std::cout << i << "\n";
    }
    for(auto i : enrg.getUBounds()) {
        std::cout << i << "\n";
    }
    return 0;
}
