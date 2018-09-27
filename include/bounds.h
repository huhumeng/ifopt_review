#pragma once

namespace ifopt{

/**
 * @brief Upper and lower bound for optimization variables and constraints.
 */
struct Bounds{

    Bounds(double lower = 0.0, double upper = 0.0)
    {
        lower_ = lower;
        upper_ = upper;
    }

    double lower_;
    double upper_;

    void operator+=(double scalar)
    {
        lower_ += scalar;
        upper_ += scalar;
    }

    void operator-=(double scalar)
    {
        lower_ -= scalar;
        upper_ -= scalar;
    }

};

// settings this as signals infinity for IPOPT/SNOPT solvers
static const double inf = 1.0e20;

static const Bounds NoBound          = Bounds(-inf, +inf);
static const Bounds BoundZero        = Bounds( 0.0,  0.0);
static const Bounds BoundGreaterZero = Bounds( 0.0, +inf);
static const Bounds BoundSmallerZero = Bounds(-inf,  0.0);

}