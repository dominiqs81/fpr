# Common boilerplate
cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE INTERNAL "" FORCE)

# Project Name
project(fixproprep)

# Prevent in-source builds
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory and run CMake from there")
endif()

# Determine if this is built as a subproject (using add_subdirectory)
# or if it is the master project.
set(MASTER_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(MASTER_PROJECT ON)
  message(STATUS "${PROJECT_NAME} is MASTER PROJECT")
  message(STATUS "CMake version: ${CMAKE_VERSION}")
endif()

if (MASTER_PROJECT)
  find_package(Git QUIET)
  # Setup git support
  if (GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/GitUtils.cmake)
    GitUtilsSetup()
  endif()
endif()

# main library
add_library(${PROJECT_NAME} STATIC src/vectorlist.cpp src/mip.cpp src/repair.cpp src/heurmgr.cpp
  src/cliquetable.cpp src/cliquecover.cpp src/impltable.cpp
  src/propagation.cpp src/table_propagators.cpp src/linear_propagator.cpp)
add_dependencies(${PROJECT_NAME} AlwaysCheckGit)

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(${PROJECT_NAME} PUBLIC Utils::Lib fmt::fmt)

add_library(FPR::Lib ALIAS ${PROJECT_NAME})

# Define fpr executable
add_executable(fpr src/main.cpp)

target_include_directories(fpr PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(fpr FPR::Lib Utils::Lib fmt::fmt)

# Define tryfpr executable
add_executable(tryfpr src/trial_main.cpp)

target_include_directories(tryfpr PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(tryfpr FPR::Lib Utils::Lib fmt::fmt)

# Define solchecker executable
add_executable(solchecker src/solchecker.cpp)

target_include_directories(solchecker PRIVATE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)
target_link_libraries(solchecker FPR::Lib Utils::Lib fmt::fmt)

# Add subprojects if this is master project
if (MASTER_PROJECT)
  add_subdirectory(extern/fmt)
  add_subdirectory(extern/utils)
endif()
