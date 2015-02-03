# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()


# This project is designed to be built outside the Insight source tree.
project(IsotropicAnomalousDiffusionImageFilter)

# Find ITK.
find_package(ITK REQUIRED)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR
            "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
  ENDIF(ITK_FOUND)

add_executable(IsotropicAnomalousDiffusionImageFilter
Main.cpp
IsotropicAnomalousDiffusionImageFilter.h
IsotropicAnomalousDiffusionImageFilter.hxx
)


target_link_libraries(IsotropicAnomalousDiffusionImageFilter ${ITK_LIBRARIES} )
