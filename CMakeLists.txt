cmake_minimum_required(VERSION 3.1)
get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})
set(CMAKE_BUILD_TYPE debug)
#set(CMAKE_BUILD_TYPE release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  polyscope
  GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(polyscope)

include(FetchContent)
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(fmt)

FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
set(EIGEN_BUILD_DOC OFF)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
set( OFF)
FetchContent_MakeAvailable(Eigen)

FetchContent_Declare(
  libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG v2.4.0
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
FetchContent_MakeAvailable(libigl)

find_package(CGAL REQUIRED COMPONENTS Core)
# find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# Add your project files
include_directories("./src/")
file(GLOB SRCFILES3D ./src/*.cpp)
file(GLOB HEADERS *.h ./src/*.h ./lib/*.hpp ./lib/*.h)

add_executable(${PROJECT_NAME} incremental.cpp ${SRCFILES} ${SRCFILES3D} ${HEADERS})
target_link_libraries(${PROJECT_NAME} CGAL::CGAL Eigen3::Eigen polyscope fmt::fmt igl::core)
