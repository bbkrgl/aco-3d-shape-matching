cmake_minimum_required(VERSION 3.16)
project(3daco)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Libigl
option(LIBIGL_GLFW "Build target igl::glfw" ON)
include(libigl)

include(dlib)

# Add your project files
file(GLOB SRC_FILES src/*.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw dlib::dlib)
