/* 
 * File:   2dlatticedata.hpp
 * Author: medved
 *
 * Created on February 12, 2016, 3:49 PM
 */

#ifndef LATTICEDATA_HPP
#define	LATTICEDATA_HPP

#include <vector>
#include "latlimits.hpp"

namespace lattice {

    /**
     * Data for 2-dimensional lattice
     */
    struct LatticeData {
        /**
         * Number of modeling layers
         */
        int mNumLayers;

        /**
         * The maximal radius of atomic interaction, beyond this radius the energy is zero
         */
        double mRadius;

        /**
         * The length of the material piece to model
         */
        double mLength;

        /**
         * Layer's atoms, each atom is identified by an integral number
         */
        std::vector<int> mLayersAtoms;

    };

}
#endif	/* LATTICEDATA_HPP */

