cmake_minimum_required(VERSION 3.15.0)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
project(sqp LANGUAGES C)

include_directories(smpsReader)
include_directories(solverUtilities)

add_executable(sqp
  stochasticQP_main.c
  dualSolve.c
  smpsReader/smps.c
  smpsReader/rvgen.c
  smpsReader/prob.c
  fullSolve.c
  Source.c
  algo.c
  partSolve.c
  solverUtilities/solver_gurobi.c
  solverUtilities/utilities.c)

target_link_libraries(sqp PRIVATE m)
target_link_libraries(sqp PRIVATE gurobi80)

