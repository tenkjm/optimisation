/* 
 * File:   energy.hpp
 * Author: medved
 *
 * Energy computations
 * 
 * Created on August 24, 2015, 12:48 PM
 */

#ifndef PAIRPOTENTIALENERGY_HPP
#define	PAIRPOTENTIALENERGY_HPP

#include <math.h>
#include <functional>
#include "energy.hpp"
#include "latticeutils.hpp"


namespace lattice {

    /**
     * Computes pair potential energy
     */
    class PairPotentialEnergy : public Energy {
    public:

        /**
         * Constructor
         * @param model reference to the model
         * @param potent reference to the potential functional
         * the potential functional returns the energy of two atoms (given by their kind) interaction and the square of Euclidian distance between them
         */
        PairPotentialEnergy(const LatticeUtils& model, std::function <double (int, int, double) > potent) :
        mLatticeModel(model), mPotent(potent), mFixedAtoms(false), mOnceComputed(false) {
            mLBounds.resize(model.getData().mNumLayers);
            mUBounds.resize(model.getData().mNumLayers);
        }

        /**
         * Compute the energy of the material piece
         * @param x parameters, defining the layers height, displacements and interatomic distances
         * @return the energy
         */
        double energy(const double* x) {
            double E = 0;
            int i;
            for (i = 0; i < mLatticeModel.getData().mNumLayers; i++) {
                E += layerEnergy(i, x);
            }
            return E;
        }

        /**
         * Sets fixed atoms regime on or off. Fixed atoms means that the energy is computed for a set of atoms rather than a piece of material
         * @param onoff on - true, off - false
         */
        void setFixedAtoms(bool onoff) {
            mFixedAtoms = onoff;
            mOnceComputed = false;
        }

    private:

        /**
         * Computes interaction energy for a layer
         * @param i layer's number
         * @param x layer's data
         * @return energy value
         */

        double layerEnergy(int i, const double* x) {
            double E = 0;
            if (!mFixedAtoms)
                mLatticeModel.computeBounds(x, mLBounds, mUBounds);
            else if (!mOnceComputed) {
                mOnceComputed = true;
                mLatticeModel.computeBounds(x, mLBounds, mUBounds);
            }
            for (int j = mLBounds[i]; j <= mUBounds[i]; j++) {
                E += atomEnergy(i, j, x);
            }
            return E;
        }

        /**
         * Computes interaction energy for an atom (interactions with all neighbours are accounted)
         * @param i atoms' layer
         * @param j atoms' number in a layer
         * @param x layer's data
         * @return energy value
         */
        double atomEnergy(int i, int j, const double* x) const {
            double v = 0;

            /**
             * Computes interaction energy for two atoms
             */
            auto inter = [&] (int i2, int j2) {

                if (!((i == i2) && (j == j2))) {
                    double q = mLatticeModel.getSqrDistance(i, j, i2, j2, x);
                    int a1 = mLatticeModel.getData().mLayersAtoms[mLatticeModel.getReferenceLayer(i)];
                    int a2 = mLatticeModel.getData().mLayersAtoms[mLatticeModel.getReferenceLayer(i2)];
                    v += mPotent(a1, a2, q);
                }
            };

            mLatticeModel.traverseLattice(i, j, inter, x);
            return v;
        }


        std::function< double (int, int, double) > mPotent;
        const LatticeUtils mLatticeModel;
        bool mFixedAtoms;
        bool mOnceComputed;
        std::vector< int > mLBounds;
        std::vector< int > mUBounds;
    };

}

#endif	/* PPENERGY_HPP */

