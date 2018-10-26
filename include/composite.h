#pragma once

#include "bounds.h"

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace ifopt{

/**
 * @brief Interface representing either Variable, Cost or Constraint.
 *
 * Variables, costs and constraints can all be fit into the same interface
 * (Component). For example, each has a "value", which is either the actual
 * value of the variables, the constraint value g or the cost.
 * This representation provides one common interface
 * ("smallest common denominator") that can be contain either individual
 * variables/costs/constraints or a Composite of these. This pattern
 * takes care of stacking variables, ordering Jacobians and providing bounds
 * for the complete problem without duplicating code. For more information on
 * the composite pattern visit https://sourcemaking.com/design_patterns/composite
 */

class Component{
public:
    typedef std::shared_ptr<Component> Ptr;

    using Jacobian = Eigen::SparseMatrix<double, Eigen::RowMajor>; // Sparse matrix存储了非零元素与元素的索引
    using VectorXd = Eigen::VectorXd;
    using VecBound = std::vector<Bounds>;

    /**
     * @brief  Creates a component.
     * @param  num_rows  The number of rows of this components.
     * @param  name  The identifier for this component.
     *
     * The number of rows @c num_rows can represent either
     * @li number of variables in this variables set        用于表示变量的数目
     * @li number of constraints in this constraint set     用于表示约束的数目
     * @li 1 if this component represents a Cost.           表示cost时候设为1
     */

    Component(int num_rows, const std::string& name);
    virtual ~Component() = default;
    
    /**
     * @brief  Returns the "values" of whatever this component represents.
     *
     * @li For Variable this represents the actual optimization values. 取回变量的值
     * @li For Constraint this represents the constraint value g.       表示约束的值
     * @li For Cost this represents the cost value.                     表示cost的值
     */
    virtual VectorXd GetValues() const = 0;

     /**
     * @brief  Returns the "bounds" of this component.
     *
     * @li For Variable these are the upper and lower variable bound.  变量的上下界
     * @li For Constraint this represents the constraint bounds.       约束的上下界
     * @li For Cost these done't exists (set to infinity).             cost不应包含此属性
     */
    virtual VecBound GetBounds() const = 0;

    /**
     * @brief  Sets the optimization variables from an Eigen vector.
     *
     * This is only done for Variable, where these are set from the current 设置变量的值, 只对变量有效
     * values of the @ref solvers.
     */
    virtual void SetVariables(const VectorXd& x) = 0;

    /**
     * @brief  Returns derivatives of each row w.r.t. the variables
     *
     * @li For Constraint this is a matrix with one row per constraint. 约束对变量的雅克比
     * @li For a Cost this is a row vector (gradient transpose).        cost对变量的雅克比
     * @li Not sensible for Variable.                                   对变量是没有用的
     */
    virtual Jacobian GetJacobian() const = 0;

    /**
     * @brief Returns the number of rows of this component.
     */
    int GetRows() const;

    /**
     * @brief Returns the name (id) of this component.
     */
    std::string GetName() const;

    /**
     * @brief Prints the relevant information (name, rows, values) of this component.  打印相关的信息
     * @param tolerance  When to flag constraint/bound violation.                        
     * @param index_start  Of this specific variables-, constraint- or cost set.
     */
    virtual void Print(double tolerance, int& index_start) const;

    /**
     * @brief Sets the number of rows of this component.
     *
     * @attention This should correctly be done through constructor call, only
     * delay this by using @c kSpecifyLater if you have good reason.
     */
    void SetRows(int num_rows);
    
    static const int kSpecifyLater = -1; // 稍后指定rows
    
private:
    int num_rows_ = kSpecifyLater;
    std::string name_;
};

/**
 * @brief A collection of components which is treated as another Component.
 *
 * This class follows the Component interface as well, but doesn't actually
 * do any evaluation, but only stitches together the results of the
 * components it is holding. This is where multiple sets of variables,
 * constraints or costs are ordered and combined.
 *
 * See Component and Composite Pattern for more information.
 * 
 * 将各个Component拼接在一起
 */

class Composite : public Component{
public:
    typedef std::shared_ptr<Composite> Ptr;

    using ComponentVec = std::vector<Component::Ptr>;

    /**
     * @brief  Creates a Composite holding either variables, costs or constraints.
     * @param  is_cost  True if this class holds cost terms, false for all others.
     *
     * Constraints and variables append individual constraint sets and Jacobian
     * rows below one another, whereas costs terms are all accumulated to a
     * scalar value/a single Jacobian row.
     */

    Composite(const std::string& name, bool is_cost);
    virtual ~Composite() = default;

    // see Component for documentation
    VectorXd GetValues   () const override;
    Jacobian GetJacobian () const override;
    VecBound GetBounds   () const override;
    void SetVariables(const VectorXd& x) override;
    void PrintAll() const;

    /**
     * @brief  Access generic component with the specified name.
     * @param  name  The name given to the component.
     * @return A generic pointer of that component.
     */
    const Component::Ptr GetComponent(std::string name) const;

    /**
     * @brief  Access type-casted component with the specified name.
     * @param  name  The name given to the component.
     * @tparam T  Type of component.
     * @return A type-casted pointer possibly providing addtional functionality.
     */
    template<typename T> 
    std::shared_ptr<T> GetComponent(const std::string& name) const;

    /**
     * @brief Adds a component to this composite.
     */
    void AddComponent(const Component::Ptr&);

    /**
     * @brief Removes all component from this composite.
     */
    void ClearComponents();

    /**
     * @brief Returns read access to the components.
     */
    const ComponentVec GetComponents() const;

private:

    ComponentVec components_;
    bool is_cost_;    
};

// implementation of template functions
// A type-casted pointer possibly providing addtional functionality. 
template<typename T>
std::shared_ptr<T> Composite::GetComponent(const std::string& name) const
{
    Component::Ptr c = GetComponent(name);
    return std::dynamic_pointer_cast<T>(c);
}

}