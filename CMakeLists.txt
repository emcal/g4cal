# g4cal: standalone PbWO4 crystal-calorimeter Geant4 simulation (calsim) and an
# ilreco-based island reconstruction (calreco). Two independent halves, each
# behind an option so you can build the sim without ilreco or vice versa.
#
# Configure:
#   cmake -B build -DCMAKE_PREFIX_PATH=/path/to/geant4
#   cmake --build build -j
#
# ilreco: found via find_package(ilreco) if installed, otherwise fetched from
# git (see G4CAL_ILRECO_GIT / G4CAL_ILRECO_TAG below). The context API is
# required, so the ref must be new enough (github emcal/ilreco master).
cmake_minimum_required(VERSION 3.16)
project(g4cal CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(G4CAL_BUILD_SIM  "Build the Geant4 simulation calsim (needs Geant4)"        ON)
option(G4CAL_BUILD_RECO "Build the ilreco-based reconstruction calreco (needs ilreco)" ON)

# ---------------- calsim (Geant4) ----------------
if(G4CAL_BUILD_SIM)
  find_package(Geant4 REQUIRED)
  file(GLOB CALSIM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/calsim/*.cc)
  add_executable(calsim ${CALSIM_SOURCES})
  target_link_libraries(calsim PRIVATE ${Geant4_LIBRARIES})
  target_include_directories(calsim PRIVATE ${Geant4_INCLUDE_DIRS})
endif()

# ---------------- ilreco + calreco ----------------
if(G4CAL_BUILD_RECO)
  find_package(ilreco CONFIG QUIET)
  if(ilreco_FOUND)
    set(G4CAL_ILRECO_TARGET ilreco::ilreco)
    message(STATUS "g4cal: using installed ilreco")
  else()
    include(FetchContent)
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../ilreco/.git)
      set(_ilreco_default ${CMAKE_CURRENT_SOURCE_DIR}/../ilreco)  # sibling checkout, offline-friendly
    else()
      set(_ilreco_default https://github.com/emcal/ilreco)
    endif()
    set(G4CAL_ILRECO_GIT ${_ilreco_default} CACHE STRING "ilreco git repo or local path for FetchContent")
    set(G4CAL_ILRECO_TAG master CACHE STRING "ilreco git ref (must expose the context API)")
    message(STATUS "g4cal: fetching ilreco from ${G4CAL_ILRECO_GIT}@${G4CAL_ILRECO_TAG}")
    # ilreco's root CMakeLists gates its tests behind PROJECT_IS_TOP_LEVEL, so an
    # embedded build pulls in just the library target (no Catch2 fetch).
    FetchContent_Declare(ilreco
      GIT_REPOSITORY ${G4CAL_ILRECO_GIT}
      GIT_TAG        ${G4CAL_ILRECO_TAG})
    FetchContent_MakeAvailable(ilreco)
    set(G4CAL_ILRECO_TARGET ilreco::ilreco)
  endif()

  # Thin adapter over ilreco's context API (see src/ilreco_adapter/adapter.h).
  add_library(calreco_lib STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/ilreco_adapter/adapter.c)
  target_include_directories(calreco_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/ilreco_adapter)
  target_link_libraries(calreco_lib PUBLIC ${G4CAL_ILRECO_TARGET} m)

  # calreco: array size is runtime configuration (--nx/--ny); the calreco_NxN
  # binaries are identical and kept only so production scripts' paths stay valid.
  enable_testing()
  foreach(N 3 10 100)
    add_executable(calreco_${N}x${N} ${CMAKE_CURRENT_SOURCE_DIR}/src/calreco/reco_main.cc)
    target_link_libraries(calreco_${N}x${N} PRIVATE calreco_lib)

    add_executable(test_adapter_${N} ${CMAKE_CURRENT_SOURCE_DIR}/src/tests/test_adapter.cc)
    target_compile_definitions(test_adapter_${N} PRIVATE TEST_GRID_N=${N})
    target_link_libraries(test_adapter_${N} PRIVATE calreco_lib)
    add_test(NAME adapter_${N} COMMAND test_adapter_${N}
             ${CMAKE_CURRENT_SOURCE_DIR}/data/prof_pwo.dat)
  endforeach()
endif()
