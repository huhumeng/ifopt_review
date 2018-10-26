#include "phase_durations_observer.h"
#include "phase_durations.h"

namespace towr{

PhaseDurationsObserver::PhaseDurationsObserver (PhaseDurationsSubjectPtr subject)
{
    phase_durations_ = subject;

    // register this observer to subject so this class always up-to-date
    subject->AddObserver(this);
}

} // namespace towr