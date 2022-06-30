cmake_minimum_required(VERSION 3.16)
project(3daco)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#set (CMAKE_CXX_STANDARD 20)
# Libigl
option(LIBIGL_GLFW "Build target igl::glfw" ON)
include(libigl)

#include(FetchContent)
#FetchContent_Declare(
#  or-tools
#  GIT_REPOSITORY https://github.com/google/or-tools.git
#  GIT_TAG        master
#)

# After the following call, the CMake targets defined by or-tools
# will be defined and available to the rest of the build
#FetchContent_MakeAvailable(or-tools)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.45.0) 

if(Boost_FOUND)
    	include_directories(${Boost_INCLUDE_DIRS})
	# Add your project files
	file(GLOB SRC_FILES src/*.cpp)
	add_executable(${PROJECT_NAME} ${SRC_FILES})
	target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw ${Boost_LIBRARIES})
endif()
