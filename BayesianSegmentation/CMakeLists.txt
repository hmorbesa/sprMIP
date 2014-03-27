CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(BayesianSegmentation)
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

include_directories(../MRIntensityNormalization)

add_executable(bayes_train sprBayesTrain.cxx )
target_link_libraries(bayes_train ${ITK_LIBRARIES})
