#include "state.h"

namespace towr{

State::State(int dim, int n_derivatives)
{
    values_ = std::vector<VectorXd>(n_derivatives, VectorXd::Zero(dim));
}

const Eigen::VectorXd State::at(Dx deriv) const
{
    return values_.at(deriv);
}

Eigen::VectorXd& State::at(Dx deriv)
{
    return values_.at(deriv);
}

const Eigen::VectorXd State::p() const
{
    return at(kPos);
}

const Eigen::VectorXd State::v () const
{
    return at(kVel);
}

const Eigen::VectorXd State::a () const
{
    return at(kAcc);
}

} // namespace towr