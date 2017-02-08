cmake_minimum_required(VERSION 2.8)
project(benchmarks)

set(CMAKE_CXX_COMPILER g++-6)
set(CMAKE_CXX_FLAGS "-O3 -g -std=c++14 -Wall -Wextra -Wundef -Wshadow -pedantic -fopenmp")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

add_definitions(-DOMP)
file(GLOB SRC
	"src/*.h"
	"src/*.hpp"
	)

add_executable(evo ${SRC} src/mainevo.cpp)
add_executable(dbg ${SRC} src/maindbg.cpp)
add_executable(flappyKeyboard ${SRC} src/mainflappykeyboard.cpp)
target_link_libraries(flappyKeyboard ncurses)
target_link_libraries(dbg ncurses)

