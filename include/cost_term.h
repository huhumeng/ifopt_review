#pragma once

#include "constraint_set.h"

namespace ifopt{

/**
 * @brief A container holding a single cost term.
 *
 * This container builds a scalar cost term from the values of the variables.
 * This can be seen as a constraint with only one row and no bounds.
 *
 * @ingroup ProblemFormulation
 * @sa Component
 */
class CostTerm : public ConstraintSet {
public:

    CostTerm(const std::string& name);
    virtual ~CostTerm() = default;

private:
    /**
     * @brief  Returns the scalar cost term calculated from the @c variables.
     */
    virtual double GetCost() const = 0;

    public:
    /**
     * @brief  Wrapper function that converts double to Eigen::VectorXd.
     */
    VectorXd GetValues() const final;

    /**
     * @brief  Returns infinite bounds (e.g. no bounds).
     */
    VecBound GetBounds() const final;

    /**
     * Cost term printout slightly different from variables/constraints.
     */
    void Print(double tol, int& index) const final;
};


}