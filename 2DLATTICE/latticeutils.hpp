/* 
 * File:   2dlatticeutils.hpp
 * Author: medved
 *
 * Created on February 12, 2016, 4:00 PM
 */

#ifndef LATTICEUTILS_HPP
#define	LATTICEUTILS_HPP

#include <math.h>
#include <functional>
#include <vector>
#include "latticedata.hpp"


namespace lattice {

    /**
     * Utilities for working with lattice data
     * Lattice parameters are stored as follows: h0 d0 s0 h1 d1 s1 .... hn-1 dn-1 sn-1
     * where hi - the difference between layers i and i - 1, di - displacement of the first atom in i-th row, si - stride in the i-th row
     * h0 is the distance between last layer and the 0-th layer in the next periodic group
     */
    class LatticeUtils {
    public:

        /**
         * Constructor
         */
        LatticeUtils(const LatticeData& data) : mData(data) {
        }

        /**
         * Computes layers height by its number
         * @param i the number of the layer
         * @param x the layer's parameters
         * @return layer's height from the 0 level (may be positive, zero or negative)
         */
        double getLayerHeight(int i, const double* x) const {
            double h = 0;
            if (i > 0) {
                for (int j = 1; j <= i; j++) {
                    int J = getReferenceLayer(j);
                    h += x[3 * J];
                }
            } else if (i < 0) {
                for (int j = 0; j > i; j--) {
                    int J = getReferenceLayer(j);
                    h -= x[3 * J];
                }
            }
            return h;
        }

        /**
         * Computes the reference layer in model data for an arbitrary row
         * taking into account periodicity
         * @param i row number
         * @return reference layer
         */
        int getReferenceLayer(int i) const {
            int k = mData.mNumLayers;
            int m = floor((double) i / (double) k);
            return i - m * k;
        }

        /**
         * Computes the displacement of the first atom and the interatomic distance in a layer
         * @param l layers number
         * @param d displacment
         * @param x layers data
         * @param s stride
         */
        void getDisplacementAndStride(int l, double& d, double& s, const double* x) const {
            int r = getReferenceLayer(l);
            d = x[3 * r + 1];
            s = x[3 * r + 2];
        }

        /**
         * Retrieve atoms 'x' coordinate (positive or negative offset from zero)
         * @param i atoms layer
         * @param j atoms number
         * @param x layers data
         * @return offset
         */
        double getOffset(int i, int j, const double* x) const {
            double d, s;
            getDisplacementAndStride(i, d, s, x);
            return d + j * s;
        }

        /**
         * Computes the square of euclidian distance between two atoms given by their 
         * layer number and number in a layer
         * @param i1 first atoms layer
         * @param j1 first atoms number
         * @param i2 second atoms layer
         * @param j2 second atoms number
         * @return distance
         */
        double getSqrDistance(int i1, int j1, int i2, int j2, const double* x) const {
            double y1 = getLayerHeight(i1, x);
            double x1 = getOffset(i1, j1, x);
            double y2 = getLayerHeight(i2, x);
            double x2 = getOffset(i2, j2, x);
            return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
        }

        /**
         * Computes scalar multiple of vectors r12 and r13
         * @param i1 first atoms layer
         * @param j1 first atoms number
         * @param i2 second atoms layer
         * @param j2 second atoms number
         * @param i3 third atoms layer
         * @param j3 third atoms number
         * @param x
         * @return 
         */
        double getScalarMult(int i1, int j1, int i2, int j2, int i3, int j3, const double* x) const {
            double y1 = getLayerHeight(i1, x);
            double x1 = getOffset(i1, j1, x);
            double y2 = getLayerHeight(i2, x);
            double x2 = getOffset(i2, j2, x);
            double y3 = getLayerHeight(i3, x);
            double x3 = getOffset(i3, j3, x);
            double scalv = (x2 - x1) * (x3 - x1) + (y2 - y1) * (y3 - y1);
            return scalv;
        }

        /**
         * Traverses the neighbourhood (including atom itself) of a given atom defined by the model cut radius
         * @param i atom's layer
         * @param j atom's position in a layer
         * @param compf functor to call for each atom
         */
        void traverseLattice(int i, int j, std::function <void (int, int) > compf, const double* x) const {
            double myx = getOffset(i, j, x);
            double myy = getLayerHeight(i, x);
            /**
             * Computes interaction energy with the layer
             */
            auto lenerg = [&] (int l) {
                double d;
                double s;
                getDisplacementAndStride(l, d, s, x);
                int tl = ceil((myx - mData.mRadius - d) / s);
                int tu = floor((myx + mData.mRadius - d) / s);
                double u = 0;
                for (int t = tl; t <= tu; t++) {
                    compf(l, t);
                }
                return u;
            };


            /**
             * Add energy of 'my' layer
             */
            lenerg(i);
            /**
             * Pass up              
             */
            for (int l = i + 1;; l++) {
                double h = getLayerHeight(l, x);
                if (h - myy > mData.mRadius)
                    break;
                lenerg(l);
            }
            /**
             * Add pass down
             */
            for (int l = i - 1;; l--) {
                double h = getLayerHeight(l, x);
                if (myy - h > mData.mRadius)
                    break;
                lenerg(l);
            }
        };

        /**
         * Computes bounds of atom indeces based on the rectangle shape
         * @param x lattice parameters
         * @param lbounds lower indices
         * @param rbounds upper indices
         */
        void computeBounds(const double* x, std::vector<int> &lbounds, std::vector<int> &rbounds) const {
            for (int i = 0; i < mData.mNumLayers; i++) {
                int a = 0, b = -1;
                int j = -1;
                while (getOffset(i, j, x) >= 0) {
                    a = j;
                    j--;
                }
                j = 0;
                while (getOffset(i, j, x) <= mData.mLength) {
                    b = j;
                    j++;
                }
                lbounds[i] = a;
                rbounds[i] = b;
            }
        }

        /**
         * Retrieve reference to lattice data
         * @return reference to lattice data
         */
        const LatticeData& getData() const {
            return mData;
        }

    private:

        const LatticeData mData;

    };
}

#endif	/* LATTICEUTILS_HPP */

