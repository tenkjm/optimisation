/* 
 * File:   energy.hpp
 * Author: medved
 *
 * Created on October 23, 2015, 5:01 PM
 */

#ifndef ENERGY_HPP
#define	ENERGY_HPP

namespace lattice {

    /**
     * Abstract class of the energy of the material piece
     */
    class Energy {
    public:

        /**
         * Compute the energy of the material piece
         * @param x parameters, defining the layers height, displacements and interatomic distances
         * @return the energy
         */
        virtual double energy(const double* x) = 0;
        
        /**
         * Sets fixed atoms regime on or off. Fixed atoms means that the energy is computed for a set of atoms rather than a piece of material
         * @param onoff on - true, off - false
         */
        virtual void setFixedAtoms(bool onoff) = 0;
    };
}

#endif	/* ENERGY_HPP */

