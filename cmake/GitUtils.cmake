# paths management
# we have to deal with two cases
# 1) file included by CMakeLists.txt, where all the usual variables are defined
# 2) file executed as a standalone script by make
if (NOT DEFINED my_source_dir)
    set(my_source_dir ${CMAKE_CURRENT_SOURCE_DIR})
endif ()

if (NOT DEFINED my_build_dir)
    set(my_build_dir ${CMAKE_CURRENT_BINARY_DIR})
endif ()

set(pre_configure_file ${my_source_dir}/include/${PROJECT_NAME}/version.h.in)
set(post_configure_file ${my_build_dir}/include/${PROJECT_NAME}/version.h)

# This can only happen in script mode
if (NOT GIT_FOUND)
  find_package(Git QUIET)
endif()

# Get current hash and branch from Git
function(CheckGitVersion)
  # Get the current working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${my_source_dir}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --always --dirty
    WORKING_DIRECTORY ${my_source_dir}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  configure_file(
    ${pre_configure_file}
    ${post_configure_file}
  )
endfunction()

# Update Git submodules (if any)
function(UpdateSubmodules)
  message(STATUS "Submodules update")
  execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init
                  WORKING_DIRECTORY ${my_source_dir}
                  RESULT_VARIABLE GIT_SUBMOD_RESULT)
  if(NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(FATAL_ERROR "git submodule update failed")
  endif()
endfunction()


function(GitUtilsSetup)
    UpdateSubmodules()

    # Phony target to keep git hash always up to date
    add_custom_target(AlwaysCheckGit COMMAND ${CMAKE_COMMAND}
        -DRUN_CHECK_GIT_VERSION=1
        -DPROJECT_NAME=${PROJECT_NAME}
        -Dmy_source_dir=${my_source_dir}
        -Dmy_build_dir=${my_build_dir}
        -P ${my_source_dir}/cmake/GitUtils.cmake
        BYPRODUCTS ${post_configure_file}
        )

    CheckGitVersion()
endfunction()


# This is used to run this function from an external cmake process.
if (RUN_CHECK_GIT_VERSION)
    CheckGitVersion()
endif ()
