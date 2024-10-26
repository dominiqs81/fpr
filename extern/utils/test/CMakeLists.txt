# Little helper function to create tests
function(add_my_test name)
  add_executable(${name} ${name}.cpp ${ARGN})
  target_include_directories(test_floats PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  )
  target_link_libraries(${name} PUBLIC utils)
  add_test(NAME ${name} COMMAND ${name} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)
endfunction()

add_my_test(test_floats)
add_my_test(test_fileconfig)
add_my_test(test_index_partition)
add_my_test(test_maths)
add_my_test(test_numbers)
add_my_test(test_sorting)
add_my_test(test_consolelog)
