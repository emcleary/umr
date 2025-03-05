if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/alglib-cpp")
  message("ALGLIB exists, skipping download")
else()
  message("Downloading ALGLIB")
  file(DOWNLOAD "https://www.alglib.net/translator/re/alglib-4.04.0.cpp.gpl.zip"
    SHOW_PROGRESS "${CMAKE_CURRENT_SOURCE_DIR}/alglib-cpp.zip")
  file(ARCHIVE_EXTRACT
    INPUT "${CMAKE_CURRENT_SOURCE_DIR}/alglib-cpp.zip"
    DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

add_definitions(-Wno-strict-aliasing)
add_definitions(-Wno-maybe-uninitialized)
add_definitions(-Wno-volatile)
add_definitions(-Wno-unused-but-set-variable)

add_library(AlgLib STATIC
  alglib-cpp/src/alglibinternal.cpp
  alglib-cpp/src/alglibinternal.h
  alglib-cpp/src/alglibmisc.cpp
  alglib-cpp/src/alglibmisc.h
  alglib-cpp/src/ap.cpp
  alglib-cpp/src/ap.h
  alglib-cpp/src/integration.cpp
  alglib-cpp/src/integration.h
  alglib-cpp/src/interpolation.cpp
  alglib-cpp/src/interpolation.h
  alglib-cpp/src/linalg.cpp
  alglib-cpp/src/linalg.h
  alglib-cpp/src/optimization.cpp
  alglib-cpp/src/optimization.h
  alglib-cpp/src/solvers.cpp
  alglib-cpp/src/solvers.h
  alglib-cpp/src/specialfunctions.cpp
  alglib-cpp/src/specialfunctions.h
)
