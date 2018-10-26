# towr cmake framework

include_directories(${PROJECT_SOURCE_DIR}/example/towr/include/constraints)
include_directories(${PROJECT_SOURCE_DIR}/example/towr/include/variables)


add_library(towr_core SHARED
    src/state.cc
    src/polynomial.cc
    src/spline.cc
    src/node_observer.cc
    src/nodes_variables.cc
    src/phase_durations_observer.cc
    src/phase_durations.cc
    
)