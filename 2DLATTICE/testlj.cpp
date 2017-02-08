#include <iostream>
#include "ppenergy.hpp"
#include "atoms.hpp"
#include "pairpotentials.hpp"

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
    double qcut = data.mRadius * data.mRadius;
    double d = 0.15;
    double qmin = (data.mRadius - d) * (data.mRadius - d);
    double qdelta = qcut - qmin;
    // Lennard Jones
    lattice::PotentialCutter pc(qcut, qdelta, lattice::ljpotent);

    /**
     * Setup energy
     */
    lattice::LatticeUtils lut(data);
    lattice::PairPotentialEnergy enrg(lut, pc);

    /**
     * Setup lattice parameters
     */
    double x[] = {1, 0, 1, 1, 0.5, 1, 1, 0, 1, 1, 0.5, 1};

    double v = enrg.energy(x);
    
    std::cout << "v = " << v << "\n";

    return 0;
}
