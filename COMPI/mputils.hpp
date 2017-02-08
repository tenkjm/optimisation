/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   mputil.hpp
 * Author: mikhail
 *
 * Created on 19 февраля 2016 г., 13:56
 */

#ifndef MPUTILS_HPP
#define MPUTILS_HPP

#include <box/boxutils.hpp>
#include "mpproblem.hpp"
#include "funcproj.hpp"

namespace COMPI {

    /**
     * Useful utilities for MPProblem
     */
    class MPUtils {
    public:

        /**
         * MP problem classes
         */
        struct ProblemTypes {
            /**
             * Generic
             */
            static const unsigned int GENERIC = 0;
            /**
             * Set if only box constraint presents
             */
            static const unsigned int BOXCONSTR = 1;
            /**
             * Set if problem is single objective
             */
            static const unsigned int SINGLEOBJ = 1 << 1;
            /**
             * Set if all variables are continuous
             */
            static const unsigned int CONTINUOUS = 1 << 2;

        };

        /**
         * Evaluates the problem type
         * @param prob problem
         * @return problem type
         */
        template <class FT> static unsigned int getProblemType(const MPProblem<FT>& prob) {
            unsigned int rv = ProblemTypes::GENERIC;
            if (prob.mObjectives.size() == 1) {
                rv |= ProblemTypes::SINGLEOBJ;
            }
            if (prob.mEqConstr.empty() && prob.mIneqConstr.empty()) {
                rv |= ProblemTypes::BOXCONSTR;
            }
            bool cont = true;
            for (auto o : prob.mVarTypes) {
                if (o != MPProblem<FT>::VariableTypes::GENERIC) {
                    cont = false;
                    break;
                }
            }
            if (cont)
                rv |= ProblemTypes::CONTINUOUS;
            return rv;
        }

        /**
         * Checks feasibility
         * @param prob problem
         * @param x point
         * @return true if x is feasible, false otherwise
         */
        template <class FT> static bool isFeasible(const MPProblem<FT>& prob, FT* x) {
            if (!snowgoose::BoxUtils::isIn(x, *(prob.mBox)))
                return false;
            for (auto c : prob.mEqConstr) {
                if (c->func(x) != 0)
                    return false;
            }
            for (auto c : prob.mIneqConstr) {
                if (c->func(x) > 0)
                    return false;
            }
            return true;
        }

        /**
         * Constructs problems projection to the given direction and point
         * @param prob problem to project
         * @param xinit initial point
         * @param direction direction of projection
         * @return resulting problem
         */
        template <class FT> MPProblem<FT>& getProblemProjection(const MPProblem<FT>& prob, const FT *xinit, const FT *direction) {
            MPProblem<FT> *problem = new MPProblem<FT>();
            const int n = prob.mBox->mDim;
            auto project = [&](const std::vector< Functor<FT>*>& oldvec, std::vector< Functor<FT>*>& newvec) {
                for (auto f : oldvec) {
                    Functor<FT>* pfunc = new FunctorProjector<FT>(*f, n, xinit, direction);
                    newvec.push_back(pfunc);
                }
            };
            project(prob.mObjectives, problem->mObjectives);
            project(prob.mIneqConstr, problem->mIneqConstr);
            project(prob.mEqConstr, problem->mEqConstr);
            // PROJECTION OF THE BOX SHOULD BE ADDED
            return *problem;
        }
    };
}

#endif /* MPUTIL_HPP */

