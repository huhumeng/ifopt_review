#include "composite.h"

#include <iostream>
#include <iomanip>

namespace ifopt{

Component::Component (int num_rows, const std::string& name)
{
    num_rows_ = num_rows;
    name_ = name;
}

int Component::GetRows () const
{
    return num_rows_;
}

void Component::SetRows (int num_rows)
{
    num_rows_ = num_rows;
}

std::string Component::GetName () const
{
    return name_;
}

void Component::Print (double tol, int& index) const
{
     // calculate squared bound violation
    VectorXd x = GetValues();
    VecBound bounds = GetBounds();

    std::vector<int> viol_idx;
    for (uint i=0; i<bounds.size(); ++i) {
        double lower = bounds.at(i).lower_;
        double upper = bounds.at(i).upper_;
        double val = x(i);
        if (val < lower-tol || upper+tol < val)
        viol_idx.push_back(i); // constraint out of bounds
    }

    std::string black = "\033[0m";
    std::string red   = "\033[31m";
    std::string color = viol_idx.empty()? black : red;

    std::cout.precision(2);
    std::cout << std::fixed
                << std::left
                << std::setw(30) << name_
                << std::right
                << std::setw(4) << num_rows_
                << std::setw(9) << index
                << std::setfill ('.')
                << std::setw(7) << index+num_rows_-1
                << std::setfill (' ')
                << color
                << std::setw(12) << viol_idx.size()
                << black
                << std::endl;

    index += num_rows_;
}

Composite::Composite (const std::string& name, bool is_cost) :Component(0, name)
{
    is_cost_ = is_cost;
}

void Composite::AddComponent (const Component::Ptr& c)
{
    // at this point the number of rows must be specified.
    assert(c->GetRows() != kSpecifyLater);

    components_.push_back(c);

    if (is_cost_)
        SetRows(1);
    else
        SetRows(GetRows()+ c->GetRows());
}

void Composite::ClearComponents ()
{
    components_.clear();
    SetRows(0);
}

const Component::Ptr Composite::GetComponent (std::string name) const
{
    for (const auto& c : components_)
        if (c->GetName() == name)
        return c;

    assert(false); // component with name doesn't exist, abort program
    return Component::Ptr();
}

Composite::VectorXd Composite::GetValues() const
{
    VectorXd g_all = VectorXd::Zero(GetRows());

    int row = 0;
    for (const auto& c : components_) {
        int n_rows = c->GetRows();
        VectorXd g = c->GetValues();
        g_all.middleRows(row, n_rows) += g;

        if (!is_cost_)
        row += n_rows;
    }
    return g_all;
}

void Composite::SetVariables(const VectorXd& x)
{
    int row = 0;
    for (auto& c : components_) {
        int n_rows = c->GetRows();
        c->SetVariables(x.middleRows(row,n_rows));
        row += n_rows;
    }
}

Composite::Jacobian Composite::GetJacobian() const
{
    int n_var = components_.front()->GetJacobian().cols();
    Jacobian jacobian(GetRows(), n_var);

    int row = 0;
    std::vector<Eigen::Triplet<double>> triplet_list;

    for(const auto& c : components_){
        const Jacobian& jac = c->GetJacobian();
        triplet_list.reserve(triplet_list.size()+jac.nonZeros());

        for (int k=0; k<jac.outerSize(); ++k)
        for (Jacobian::InnerIterator it(jac,k); it; ++it)
            triplet_list.push_back(Eigen::Triplet<double>(row+it.row(), it.col(), it.value()));

        if (!is_cost_)
        row += c->GetRows();
    }

    jacobian.setFromTriplets(triplet_list.begin(), triplet_list.end());
    return jacobian;
}

Composite::VecBound Composite::GetBounds() const
{
    VecBound bounds_;
    for (const auto& c : components_) {
        VecBound b = c->GetBounds();
        bounds_.insert(bounds_.end(), b.begin(), b.end());
    }

    return bounds_;
}

const Composite::ComponentVec Composite::GetComponents() const
{
    return components_;
}

void Composite::PrintAll() const
{
    int index = 0;
    double tol = 0.001; ///< tolerance when printing out constraint/bound violation.

    std::cout << GetName() << ":\n";
    for (auto c : components_) {
        std::cout << "   "; // indent components
        c->Print(tol, index);
    }
    std::cout << std::endl;
}


}