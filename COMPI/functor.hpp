/* 
 * File:   functor.hpp
 * Author: medved
 *
 * Created on February 17, 2016, 3:51 PM
 */

#ifndef FUNCTOR_HPP
#define	FUNCTOR_HPP

#include <common/sgerrcheck.hpp>


namespace COMPI {

    /**
     * Abstract class for storing functional objectives and constraints
     */
    template <class FT> class Functor {
    public:

        /**
         * Compute function
         * @param x argument
         * @return function value
         *
         */
        virtual FT func(const FT* x) = 0;

        /**
         * Compute gradient (optional)        
         * @param x argument
         * @param g gradient
         */
        virtual void grad(const FT* x, FT* g) {
            SG_ERROR_REPORT("Gradient not implemented");
        }

        /**
         * Compute Hessian (optional)
         * @param x point
         * @param H hessian 
         */
        virtual void hess(const FT* x, FT* H) {
            SG_ERROR_REPORT("Hessian not implemented");
        }
    };

}

#endif	/* FUNCTOR_HPP */

