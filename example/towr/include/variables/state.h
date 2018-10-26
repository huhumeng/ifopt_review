#pragma once

#include <vector>

#include <Eigen/Dense>

namespace towr{

///< the values or derivative. For motions e.g. position, velocity, ...
enum Dx { kPos=0, kVel, kAcc, kJerk };

/**
 * @brief Stores at state comprised of values and higher-order derivatives.
 *
 * This state can represent a motion state with position, velocity and
 * accelerations, but also a force-profiles with forces, force-derivatives etc.
 */

class State{
public:

    using VectorXd = Eigen::VectorXd;

    /**
     * @brief Constructs a state object.
     *
     * @param dim  The number of dimensions this state has (e.g. x,y,z).
     * @param n_derivatives  The number of derivatives. In control a state
     *                       is usually made up of two (positions and velocities.
     */
    explicit State(int dim, int n_derivatives);
    virtual ~State() = default;


    /**
     * @brief   Read the state value or it's derivatives by index.
     * @param   deriv  Index for that specific derivative (pos=0, vel=1, acc=2).
     * @return  Read-only n-dimensional position, velocity or acceleration.
     */
    const VectorXd at(Dx deriv) const;

    /**
     * @brief   Read or write a specific state derivative by index.
     * @param   deriv  Index for that specific derivative (pos=0, vel=1, acc=2).
     * @return  Read/write n-dimensional position, velocity or acceleration.
     */
    VectorXd& at(Dx deriv);

    /**
     * @brief read access to the zero-derivative of the state, e.g. position.
     */
    const VectorXd p() const;

    /**
     * @brief read access to the first-derivative of the state, e.g. velocity.
     */
    const VectorXd v() const;

    /**
     * @brief read access to the second-derivative of the state, e.g. acceleration.
     */
    const VectorXd a() const;

private:
    std::vector<VectorXd> values_; ///< e.g. position, velocity and acceleration, ...
};

/**
 * @brief A node represents the state of a trajectory at a specific time.
 *
 * Given a set of nodes, cubic polynomials can be used to smoothly interpolate
 * between them. Therefore, if optimal node values have been found, the
 * continuous trajectory for that spline can be reconstructed.
 *
 * In this framework a node only has position and velocity values, no
 * acceleration.
 */
class Node : public State {
public:

    static const int n_derivatives = 2; ///< value and first derivative.

    /**
     * @brief Constructs a @a dim - dimensional node (default zero-dimensional).
     */
    explicit Node(int dim = 0) : State(dim, n_derivatives) {};
    virtual ~Node() = default;
};

/**
 * @brief Can represent the 6Degree-of-Freedom floating base of a robot.
 */
struct BaseState {
    
    BaseState(): lin(3), ang(3) {}

    Node lin; ///< linear position x,y,z and velocities.
    Node ang; ///< angular euler roll, pitch, yaw and rates.
};

}