/* 
 * File:   mpproblem.hpp
 * Author: medved
 *
 * Created on February 17, 2016, 2:09 PM
 */

#ifndef MPPROBLEM_HPP
#define MPPROBLEM_HPP

#include <vector>
#include <cstddef>
#include <string>
#include <box/box.hpp>
#include "functor.hpp"

namespace COMPI {

    /**
     * Generic interface to mathematical programming problem
     */
    template <class FT> class MPProblem {
    public:

        MPProblem() : mBox(nullptr) {
        }

        /**
         * Information about the problem
         * @return Problem description
         */
        virtual std::string about() const {
            return "No problem description provided";
        }

        /**
         * Objective functions (maybe more than one if the problem is multicriterial)
         */
        std::vector< Functor<FT>* > mObjectives;

        /**
         * Inequality constraints in the form g(x) <= 0
         */
        std::vector< Functor<FT>* > mIneqConstr;

        /**
         * Equality constraints in the form g(x) = 0
         */
        std::vector< Functor<FT>* > mEqConstr;

        /**
         * Bounding box
         */
        snowgoose::Box<FT>* mBox;

#if 0

        /**
         * Types of variables
         */
        struct VariableTypes {
            /**
             * Normal continuos variable
             */
            static const unsigned int GENERIC = 0;

            /**
             * Integral
             */
            static const unsigned int INTEGRAL = 1;

            /**
             * Boolean
             */
            static const unsigned int BOOLEAN = 2;
        };
#endif

        enum VariableTypes {
            GENERIC,
            INTEGRAL,
            BOOLEAN
        };

        /**
         * Characteristics of variables
         */
        //std::vector<unsigned int> mVarTypes;
        std::vector<VariableTypes> mVarTypes;
    };

}

#endif /* MPPROBLEM_HPP */

