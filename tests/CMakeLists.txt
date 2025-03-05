add_subdirectory(unit)
add_subdirectory(regression)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests)

add_executable("test_unit" unit/main.cpp)
target_link_libraries("test_unit" PRIVATE TestingUnit ${Boost_LIBRARIES})
target_link_libraries("test_unit" PUBLIC ${LIBRARIES} UMR)

enable_testing()
add_test(NAME test_unit COMMAND test_unit)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/regression)
file(GLOB REGSRC regression/*/test.cpp)

foreach(filename ${REGSRC})

  get_filename_component(stem ${filename} DIRECTORY)
  get_filename_component(regression_test ${stem} NAME)
  set(libname "_${regression_test}")

  add_executable(${regression_test} ${filename})
  target_link_libraries(${regression_test} PRIVATE ${Boost_LIBRARIES})
  target_link_libraries(${regression_test} PUBLIC ${LIBRARIES})
  target_link_libraries(${regression_test} PUBLIC TestingRegression UMR)
  target_link_libraries(${regression_test} PUBLIC UMR AlgLib jsoncpp_lib ${VTK_LIBRARIES})

  # Note that automating executables and testing like this won't allow for prioritizing
  add_test(NAME ${regression_test} COMMAND ${regression_test} WORKING_DIRECTORY ${stem})

endforeach(filename ${REGSRC})
