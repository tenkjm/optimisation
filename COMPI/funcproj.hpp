/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   funcproj.hpp
 * Author: Mikhail Posypkin <mposypkin at gmail.com>
 *
 * Created on 4 мая 2016 г., 12:40
 */

#ifndef FUNCPROJ_HPP
#define FUNCPROJ_HPP

#include <common/smartptr.hpp>
#include <common/vec.hpp>
#include "functor.hpp"

namespace COMPI {

    /**
     * Function that is a projector of the given function 
     * to the given direction
     */
    template <class FT> class FunctorProjector : public Functor<FT> {
    public:

        /**
         * Constructor
         * @param f projected functor
         * @param n number of variables int the project functor
         * @param xinit the initial point 
         * @param direction projection direction
         */
        FunctorProjector(Functor <FT>& f, int n, const FT* xinit, const FT* direction) :
        mF(f), mN(n) {
            mDir.alloc(n);
            mXinit.alloc(n);
            mX.alloc(n);
            snowgoose::VecUtils::vecCopy(n, direction, (FT*) mDir);
            snowgoose::VecUtils::vecCopy(n, xinit, (FT*) mXinit);
        }

        FT func(const FT* x) {
            FT t = *x;
            snowgoose::VecUtils::vecSaxpy(mN, (FT*)mXinit, (FT*)mDir, t, (FT*)mX);
            return mF.func(mX);
        }

    private:
        Functor <FT>& mF;
        const int mN;
        snowgoose::SmartArrayPtr<FT> mDir;
        snowgoose::SmartArrayPtr<FT> mXinit;
        snowgoose::SmartArrayPtr<FT> mX;
    };

}

#endif /* FUNCPROJ_HPP */

