#pragma once

#include "problem.h"
#include "solver.h"


namespace Ipopt {
class IpoptApplication;
}

namespace ifopt {

/**
 * @brief An interface to IPOPT, fully hiding its implementation.
 *
 * To set specific options, see:
 * https://www.coin-or.org/Ipopt/documentation/node40.html
 *
 * @ingroup Solvers
 */
class IpoptSolver : public Solver{
    
public:
    typedef std::shared_ptr<IpoptSolver> Ptr;

    IpoptSolver();
    virtual ~IpoptSolver() = default;

    /** @brief  Creates an IpoptAdapter and solves the NLP.
        * @param [in/out]  nlp  The specific problem.
        */
    void Solve(Problem& nlp) override;

    /** Options for the IPOPT solver. A complete list can be found here:
        * https://www.coin-or.org/Ipopt/documentation/node40.html
        */
    void SetOption(const std::string& name, const std::string& value);
    void SetOption(const std::string& name, int value);
    void SetOption(const std::string& name, double value);

private:
    std::shared_ptr<Ipopt::IpoptApplication> ipopt_app_;
};

} /* namespace ifopt */
