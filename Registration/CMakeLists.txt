CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(Registration)
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

add_executable(registering Registration.cxx )
target_link_libraries(registering ${ITK_LIBRARIES})

add_executable(DeformableRegistrationlastversion DeformableRegistrationlastversion.cxx )
target_link_libraries(DeformableRegistrationlastversion ${ITK_LIBRARIES})
