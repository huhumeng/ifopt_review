#pragma once

#include "composite.h"

namespace ifopt{

/**
 * @brief A container holding a set of related constraints.
 *
 * This container holds constraints representing a single concept, e.g.
 * @c n constraints keeping a foot inside its range of motion. Each of the
 * @c n rows is given by:
 * lower_bound < g(x) < upper_bound
 *
 * These constraint sets are later then stitched together to form the overall
 * problem.
 *
 * @ingroup ProblemFormulation
 * @sa Component
 */
class ConstraintSet : public Component{
public:
    typedef std::shared_ptr<ConstraintSet> Ptr;

    using VariablesPtr = Composite::Ptr;

    /**
     * @brief Creates constraints on the variables @c x.
     * @param n_constraints  The number of constraints.
     * @param name  What these constraints represent.
     */
    ConstraintSet(int n_constraints, const std::string& name);
    virtual ~ConstraintSet() = default;

    /**
     * @brief Connects the constraint with the optimization variables.
     * @param x  A pointer to the current values of the optimization variables.
     *
     * The optimization variable values are necessary for calculating constraint
     * violations and Jacobians.
     */
    void LinkWithVariables(const VariablesPtr& x);

    /**
     * @brief  The matrix of derivatives for these constraints and variables.
     *
     * Assuming @c n constraints and @c m variables, the returned Jacobian
     * has dimensions n x m. Every row represents the derivatives of a single
     * constraint, whereas every column refers to a single optimization variable.
     *
     * This function only combines the user-defined jacobians from
     * FillJacobianBlock().
     */
    Jacobian GetJacobian() const final;

protected:

    /**
     * @brief Read access to the value of the optimization variables.
     *
     * This must be used to formulate the constraint violation and Jacobian.
     */
    const VariablesPtr GetVariables() const{
        return variables_; 
    }


private:
    /**
     * @brief Set individual Jacobians corresponding to each decision variable set.
     * @param var_set  Set of variables the current Jacobian block belongs to.
     * @param jac_block  Columns of the overall Jacobian affected by var_set.
     *
     * A convenience function so the user does not have to worry about the
     * ordering of variable sets. All that is required is that the user knows
     * the internal ordering of variables in each individual set and provides
     * the Jacobian of the constraints w.r.t. this set (starting at column 0).
     * GetJacobian() then inserts these columns at the correct position in the
     * overall Jacobian.
     *
     * If the constraint doen't depend on a @c var_set, this function should
     * simply do nothing.
     */
    virtual void FillJacobianBlock(std::string var_set, Jacobian& jac_block) const = 0;
    VariablesPtr variables_;

     /**
     * @brief  Initialize quantities that depend on the optimization variables.
     * @param x  A pointer to the initial values of the optimization variables.
     *
     * Sometimes the number of constraints depends on the variable representation,
     * or shorthands to specific variable sets want to be saved for quicker
     * access later. This function can be overwritten for that.
     */
    virtual void InitVariableDependedQuantities(const VariablesPtr& x_init){}

    // doesn't exist for constraints, generated run-time error when used
    void SetVariables(const VectorXd& x) final{ 
        assert(false); 
    }

};
}