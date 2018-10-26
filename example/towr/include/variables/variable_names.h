#pragma once

#include <string>

namespace towr{

/**
 * @brief Identifiers (names) used for variables in the optimization problem.
 *
 * @ingroup Variables
 */
namespace id{

static const std::string base_lin_nodes    = "base-lin";
static const std::string base_ang_nodes    = "base-ang";
static const std::string ee_motion_nodes   = "ee-motion_";
static const std::string ee_force_nodes    = "ee-force_";
static const std::string contact_schedule  = "ee-schedule";

static std::string EEMotionNodes(uint ee)
{
    return ee_motion_nodes + std::to_string(ee);
}

static std::string EEForceNodes(uint ee)
{
    return ee_force_nodes + std::to_string(ee);
}

static std::string EESchedule(uint ee)
{
    return contact_schedule + std::to_string(ee);
}

} // namespace id
} // namespace towr