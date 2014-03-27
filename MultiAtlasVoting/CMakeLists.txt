CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(MultiAtlasVoting)
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

add_executable(atlasvoting sprMultiAtlasVoting.cxx )
target_link_libraries(atlasvoting ${ITK_LIBRARIES})
