/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   funccnt.hpp
 * Author: mikhail
 *
 * Created on 26 февраля 2016 г., 11:10
 */

#ifndef FUNCCNT_HPP
#define FUNCCNT_HPP

#include "functor.hpp"

namespace COMPI {

    /**
     * The wrapper around functor that computes 
     * the number of function, gradient and hessian evaluations
     */
    template <class FT> class FuncCnt : public Functor<FT> {
    public:

        struct Counters {
            /**
             * Number of function calls
             */
            unsigned int mFuncCalls;
            /**
             * Number of gradient calls
             */
            unsigned int mGradCalls;
            /**
             * Number of Hessian calls
             */
            unsigned int mHessCalls;
        };

        /**
         * Counters
         */
        Counters mCounters;

        /**
         * Constructor
         * @param f wrapper functor
         */
        FuncCnt(Functor<FT>& f) : mF(f) {
            reset();   
        }

        FT func(const FT* x) {
            mCounters.mFuncCalls++;
            FT v = mF.func(x);
            return v;
        }

        void grad(const FT* x, FT* g) {
            mCounters.mGradCalls ++;
            mF.grad(x, g);
        }

        void hess(const FT* x, FT* H) {
            mCounters.mHessCalls ++;
            mF.hess(x, H);            
        }
        
        /**
         * Reset counters
         */
        void reset() {
            mCounters.mFuncCalls = 0;
            mCounters.mGradCalls = 0;
            mCounters.mHessCalls = 0;
        }

        
    private:
        Functor <FT>& mF;
    };

}

#endif /* FUNCCNT_HPP */

