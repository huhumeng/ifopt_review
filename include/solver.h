#pragma once

#include "problem.h"

#include <memory>

namespace ifopt {

/**
 * @defgroup Solvers
 * @brief Interfaces to IPOPT and SNOPT to solve the optimization problem.
 *
 * These are included in the folders: @ref ifopt_ipopt/ and @ref ifopt_snopt/.
 */

/**
 * @brief Solver interface implemented by IPOPT and SNOPT.
 *
 * @ingroup Solvers
 */
class Solver{
public:
    typedef std::shared_ptr<Solver> Ptr;

    virtual ~Solver () = default;

    /** @brief  Uses a specific solver (IPOPT, SNOPT) to solve the NLP.
        * @param [in/out]  nlp  The nonlinear programming problem.
        */
    virtual void Solve(Problem& nlp) = 0;
};

} /* namespace ifopt */