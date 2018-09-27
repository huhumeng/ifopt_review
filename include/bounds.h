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

    Bounds operator+(double scalar){
        return Bounds(lower_ + scalar, upper_ + scalar);
    }

    Bounds operator-(double scalar){
        return Bounds(lower_ - scalar, upper_ - scalar);
    }

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

// 无边界约束
static const Bounds NoBound          = Bounds(-inf, +inf);

// 0等式约束
static const Bounds BoundZero        = Bounds( 0.0,  0.0);

// 不等式约束
static const Bounds BoundGreaterZero = Bounds( 0.0, +inf);
static const Bounds BoundSmallerZero = Bounds(-inf,  0.0);

}