#pragma once

#include "composite.h"

namespace ifopt{

/**
 * @brief  A container holding a set of related optimization variables.
 *
 * This is a single set of variables representing a single concept, e.g
 * "spline coefficients" or "step durations".
 *
 * @ingroup ProblemFormulation
 * @sa Component
 */

class VariableSet : public Component{
public:

    /**
     * @brief Creates a set of variables representing a single concept.
     * @param n_var  Number of variables.
     * @param name   What the variables represent to (e.g. "spline coefficients").
     */
    VariableSet(int n_var, const std::string& name);
    virtual ~VariableSet() = default;

    // doesn't exist for variables, generated run-time error when used.
    Jacobian GetJacobian() const override final
    {
        throw std::runtime_error("not implemented for variables");
    }

};
}