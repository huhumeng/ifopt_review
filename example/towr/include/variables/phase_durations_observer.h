#pragma once

namespace towr{

class PhaseDurations;

/**
 * @brief Base class to receive up-to-date values of the ContactSchedule.
 *
 * This class registers with the contact schedule and everytime those
 * durations change, the contact schedule updates this class by calling the
 * UpdatePhaseDurations() method.
 *
 * Used by spline.h
 *
 * This class implements the observer pattern:
 * https://sourcemaking.com/design_patterns/observer
 */

class PhaseDurationsObserver {
public:
    
    using PhaseDurationsSubjectPtr = PhaseDurations*; // observer shouldn't own subject

    PhaseDurationsObserver() = default;

    /**
     * @brief Registers this observer with the subject class to receive updates.
     * @param phase_durations  A pointer to the hase durations subject.
     */
    PhaseDurationsObserver(PhaseDurationsSubjectPtr phase_durations);
    virtual ~PhaseDurationsObserver() = default;

    /**
     * @brief Callback method called every time the subject changes.
     */
    virtual void UpdatePolynomialDurations() = 0;

protected:
    PhaseDurationsSubjectPtr phase_durations_;
};

} /* namespace towr */