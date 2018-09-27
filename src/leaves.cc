#include "variable_set.h"
#include "constraint_set.h"
#include "cost_term.h"

#include <iostream>
#include <iomanip>

namespace ifopt{

VariableSet::VariableSet(int n_var, const std::string& name)
    : Component(n_var, name)
{}

ConstraintSet::ConstraintSet(int row_count, const std::string& name)
    : Component(row_count, name)
{}

ConstraintSet::Jacobian ConstraintSet::GetJacobian() const
{
    Jacobian jacobian(GetRows(), variables_->GetRows());

    int col = 0;
    Jacobian jac;
    std::vector< Eigen::Triplet<double> > triplet_list;

    for (const auto& vars : variables_->GetComponents()) {
        int n = vars->GetRows();
        jac.resize(GetRows(), n);

        FillJacobianBlock(vars->GetName(), jac);
        // reserve space for the new elements to reduce memory allocation
        triplet_list.reserve(triplet_list.size()+jac.nonZeros());

        // create triplets for the derivative at the correct position in the overall Jacobian
        for (int k=0; k<jac.outerSize(); ++k)
        for (Jacobian::InnerIterator it(jac,k); it; ++it)
            triplet_list.push_back(Eigen::Triplet<double>(it.row(), col+it.col(), it.value()));
        col += n;
    }

    // transform triplet vector into sparse matrix
    jacobian.setFromTriplets(triplet_list.begin(), triplet_list.end());
    return jacobian;
}

void ConstraintSet::LinkWithVariables(const VariablesPtr& x)
{
    variables_ = x;
    InitVariableDependedQuantities(x);
}

CostTerm::CostTerm (const std::string& name) :ConstraintSet(1, name)
{}

CostTerm::VectorXd CostTerm::GetValues() const
{
    VectorXd cost(1);
    cost(0) = GetCost();
    return cost;
}

CostTerm::VecBound CostTerm::GetBounds() const
{
    return VecBound(GetRows(), NoBound);
}

void CostTerm::Print (double tol, int& index) const
{
    // only one scalar cost value
    double cost = GetValues()(0);

    std::cout.precision(2);
    std::cout << std::fixed
                << std::left
                << std::setw(30) << GetName()
                << std::right
                << std::setw(4) << GetRows()
                << std::setw(9) << index
                << std::setfill ('.')
                << std::setw(7) << index+GetRows()-1
                << std::setfill (' ')
                << std::setw(12) << cost
                << std::endl;
}

}