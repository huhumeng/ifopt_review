#include "nodes_observer.h"
#include "nodes_variables.h"

namespace towr {

NodesObserver::NodesObserver(NodeSubjectPtr subject)
{
    node_values_ = subject;

    // register observer to subject so this class always up-to-date
    subject->AddObserver(this);
};

}