#include "problem.h"

#include <iostream>
#include <iomanip>

namespace ifopt{

Problem::Problem()
    :constraints_("constraint-sets", false),
     costs_("cost-terms", true)
{
    variables_ = std::make_shared<Composite>("variable-sets", false);
}

void Problem::AddVariableSet(VariableSet::Ptr variable_set)
{
    variables_->AddComponent(variable_set);
}

void Problem::AddConstraintSet(ConstraintSet::Ptr constraint_set)
{
    constraint_set->LinkWithVariables(variables_);
    constraints_.AddComponent(constraint_set);
}

void Problem::AddCostSet(CostTerm::Ptr cost_set)
{
    cost_set->LinkWithVariables(variables_);
    costs_.AddComponent(cost_set);
}

int Problem::GetNumberOfOptimizationVariables () const
{
    return variables_->GetRows();
}

Problem::VecBound Problem::GetBoundsOnOptimizationVariables () const
{
    return variables_->GetBounds();
}

Problem::VectorXd Problem::GetVariableValues () const
{
    return variables_->GetValues();
}

void Problem::SetVariables (const double* x)
{
    variables_->SetVariables(ConvertToEigen(x));
}

double Problem::EvaluateCostFunction (const double* x)
{
    VectorXd g = VectorXd::Zero(1);
    if (HasCostTerms()) {
        SetVariables(x);
        g = costs_.GetValues();
    }
    return g(0);
}

Problem::VectorXd Problem::EvaluateCostFunctionGradient (const double* x)
{
    Jacobian jac = Jacobian(1,GetNumberOfOptimizationVariables());
    if (HasCostTerms()) {
        SetVariables(x);
        jac = costs_.GetJacobian();
    }

    return jac.row(0).transpose();
}

Problem::VecBound Problem::GetBoundsOnConstraints () const
{
    return constraints_.GetBounds();
}

int Problem::GetNumberOfConstraints () const
{
    return GetBoundsOnConstraints().size();
}

Problem::VectorXd Problem::EvaluateConstraints (const double* x)
{
    SetVariables(x);
    return constraints_.GetValues();
}

bool Problem::HasCostTerms () const
{
    return costs_.GetRows()>0;
}

void Problem::EvalNonzerosOfJacobian (const double* x, double* values)
{
    SetVariables(x);
    Jacobian jac = GetJacobianOfConstraints();

    jac.makeCompressed(); // so the valuePtr() is dense and accurate
    std::copy(jac.valuePtr(), jac.valuePtr() + jac.nonZeros(), values);
}

Problem::Jacobian Problem::GetJacobianOfConstraints () const
{
    return constraints_.GetJacobian();
}

void Problem::SaveCurrent()
{
    x_prev.push_back(variables_->GetValues());
}

Composite::Ptr Problem::GetOptVariables () const
{
    return variables_;
}

void Problem::SetOptVariables (int iter)
{
    variables_->SetVariables(x_prev.at(iter));
}

void Problem::SetOptVariablesFinal ()
{
    variables_->SetVariables(x_prev.at(GetIterationCount()-1));
}

void Problem::PrintCurrent() const
{
    using namespace std;
    cout << "\n"
        << "************************************************************\n"
        << "    IFOPT - Interface to Nonlinear Optimizers (v2.0)\n"
        << "                \u00a9 Alexander W. Winkler\n"
        << "           https://github.com/ethz-adrl/ifopt\n"
        << "************************************************************"
        << "\n"
        << "Legend:\n"
        << "c - number of variables, constraints or cost terms" << std::endl
        << "i - indices of this set in overall problem" << std::endl
        << "v - number of [violated variable- or constraint-bounds] or [cost term value]"
        << "\n\n"
        << std::right
        << std::setw(33) << ""
        << std::setw(5)  << "c  "
        << std::setw(16) << "i    "
        << std::setw(11) << "v "
        << std::left
        << "\n";

    variables_->PrintAll();
    constraints_.PrintAll();
    costs_.PrintAll();
}

Problem::VectorXd Problem::ConvertToEigen(const double* x) const
{
    return Eigen::Map<const VectorXd>(x,GetNumberOfOptimizationVariables());
}

} /* namespace opt */
