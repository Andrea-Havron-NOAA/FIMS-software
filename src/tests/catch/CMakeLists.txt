cmake_minimum_required(VERSION 3.14)
project(FIMS-software)

# Catch requires at least C++11
set(CMAKE_CXX_STANDARD 11)

include(FetchContent)
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG devel
)

FetchContent_MakeAvailable(catch)

list(APPEND CMAKE_MODULE_PATH ${catch_SOURCE_DIR}/extras)

enable_testing()

add_executable(
  likelihoods_test
  catch-likelihoods.cpp
)

target_include_directories(
  likelihoods_test PUBLIC
  ${CMAKE_SOURCE_DIR}/../
)

target_link_libraries(
  likelihoods_test
  Catch2::Catch2WithMain
)

add_executable(
  model_test
  catch-model.cpp
)

target_include_directories(
  model_test PUBLIC
  ${CMAKE_SOURCE_DIR}/../
)

target_link_libraries(
  model_test
  Catch2::Catch2WithMain
)

target_compile_definitions(
  model_test PUBLIC
  STD_LIB
)

include(Catch)
catch_discover_tests(likelihoods_test)
catch_discover_tests(model_test)