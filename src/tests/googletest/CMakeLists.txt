cmake_minimum_required(VERSION 3.14)
project(FIMS-software)

# GoogleTest requires at least C++11
set(CMAKE_CXX_STANDARD 11)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.11.0.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  likelihoods_test
  googletest-likelihoods.cpp
)

target_include_directories(
  likelihoods_test PUBLIC
  ${CMAKE_SOURCE_DIR}/../
)

target_link_libraries(
  likelihoods_test
  gtest_main
)

add_executable(
  model_test
  googletest-model.cpp
)

target_include_directories(
  model_test PUBLIC
  ${CMAKE_SOURCE_DIR}/../
)

target_link_libraries(
  model_test
  gtest_main
)

target_compile_definitions(
  model_test PUBLIC
  STD_LIB
)

include(GoogleTest)
gtest_discover_tests(likelihoods_test)
gtest_discover_tests(model_test)