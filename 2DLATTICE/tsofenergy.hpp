/* 
 * File:   tsofenergy.hpp
 * Author: medved
 *
 * Created on October 24, 2015, 11:59 AM
 */

#ifndef TSOFENERGY_HPP
#define	TSOFENERGY_HPP

#include <math.h>
#include <functional>
#include <vector>
#include "energy.hpp"
#include "latticeutils.hpp"
#include "tersoffparams.hpp"
#include "tersoffutils.hpp"


namespace lattice {

    /**
     * Lattice energy by tersoff potential
     */
    class TersoffEnergy : public Energy {
    public:

        /**
         * Constructor
         * @param model reference to the model   
         */
        TersoffEnergy(const LatticeUtils& utils, const TersoffUtils& tutils) :
        mLatticeUtils(utils), mTutils(tutils), mFixedAtoms(false), mOnceComputed(false) {
            mLBounds.resize(utils.getData().mNumLayers);
            mUBounds.resize(utils.getData().mNumLayers);
        }

        /**
         * Compute the energy of the material piece
         * @param x parameters, defining the layers height, displacements and interatomic distances
         * @return the energy
         */
        double energy(const double* x)  {
            double E = 0;
            for (int i = 0; i < mLatticeUtils.getData().mNumLayers; i++) {
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
        
        /**
         * Retrieve lower bounds
         * @return vector of lower bounds
         */
        const std::vector<int>& getLBounds() const {
            return mLBounds;
        }

        /**
         * Retrieve lower bounds
         * @return vector of lower bounds
         */
        const std::vector<int>& getUBounds() const {
            return mUBounds;
        }

    private:

 

        /**
         * Computes interaction energy for a layer
         * @param i layer's number
         * @param x layer's data
         * @return energy value
         */
        double layerEnergy(int i, const double* x)  {
            double E = 0;
            if (!mFixedAtoms)
                mLatticeUtils.computeBounds(x, mLBounds, mUBounds);
            else if(!mOnceComputed) {
                mOnceComputed = true;
                mLatticeUtils.computeBounds(x, mLBounds, mUBounds);
            }
            for (int j = mLBounds[i]; j <= mUBounds[i]; j++) {
                E += atomEnergy(i, j, x);
            }
            return E;
        }

        /**
         * Computes interaction energy for an atom 
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
                    double q = mLatticeUtils.getSqrDistance(i, j, i2, j2, x);
                    double r = sqrt(q);
                    double fc = mTutils.cutoff(r);
                    if (fc > 0) {
                        double zij = 0;

                        /**
                         * Aux function to compute angular term zij
                         */
                        auto ang = [&] (int i3, int j3) {                            
                            bool isij = (i == i3) && (j == j3);
                            bool isi2j2 = (i2 == i3) && (j2 == j3);
                            if (!(isij || isi2j2)) {
                                double q3 = mLatticeUtils.getSqrDistance(i, j, i3, j3, x);
                                double r3 = sqrt(q3);
                                double fc3 = mTutils.cutoff(r3);
                                if (fc3 > 0) {
                                    double o = mTutils.computeOmega(r, r3);
                                    double scalv = mLatticeUtils.getScalarMult(i, j, i2, j2, i3, j3, x);
                                    double cosv = scalv / (r * r3);
                                    double u = mTutils.computeG(cosv);
                                    zij += fc3 * u * o;
                                }
                            }
                        };
                        mLatticeUtils.traverseLattice(i, j, ang, x);
                        double bij = mTutils.computeBij(zij);
                        v += fc * (mTutils.VR(r) - bij * mTutils.VA(r));
                    }
                }
            };

            mLatticeUtils.traverseLattice(i, j, inter, x);
            return v;
        }

        const TersoffUtils mTutils;
        const LatticeUtils mLatticeUtils;
        bool mFixedAtoms;
        bool mOnceComputed;

        std::vector< int > mLBounds;
        std::vector< int > mUBounds;
    };

}


#endif	/* TSOFENERGY_HPP */

