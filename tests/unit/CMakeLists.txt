file (GLOB headers CONFIGURE_DEPENDS "*.hpp")
file (GLOB source CONFIGURE_DEPENDS "*.cpp")
add_library(TestingUnit SHARED ${source} ${headers})
target_link_libraries(TestingUnit PUBLIC ${LIBRARIES} AlgLib)
