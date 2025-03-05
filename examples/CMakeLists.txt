file(GLOB LIB_SOURCES *.cpp)
file(GLOB CSV *.csv)
file(COPY ${CSV} DESTINATION ${CMAKE_BINARY_DIR}/bin/examples)
# add_library(Examples ${LIB_SOURCES} ${LIB_HEADERS})
# target_link_libraries(Examples PUBLIC "${LIBRARIES}")
# target_link_libraries(Examples PUBLIC "${VTK_LIBRARIES}")
# target_link_libraries(Examples PUBLIC "${Boost_LIBRARIES}")

