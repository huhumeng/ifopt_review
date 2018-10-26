#pragma once

namespace towr{

class NodesVariables;

/**
 * @brief Base class to receive up-to-date values of the NodeVariables.
 *
 * This class registers with the node variables and everytime the positions or
 * velocities of a node change, the subject updates this class by calling the
 * UpdatePolynomials() method.
 *
 * Used by spline.h
 *
 * This class implements the observer pattern:
 * https://sourcemaking.com/design_patterns/observer
 */

class NodesObserver {
public:
    
    using NodeSubjectPtr = NodesVariables*; // observer shouldn't own subject

    /**
     * @brief Registers this observer with the subject class to receive updates.
     * @param node_values  The subject holding the Hermite node values.
     */
    NodesObserver(NodeSubjectPtr node_values);
    virtual ~NodesObserver() = default;

    /**
     * @brief Callback method called every time the subject changes.
     */
    virtual void UpdateNodes() = 0;

protected:
    NodeSubjectPtr node_values_;
};

} /* namespace towr */