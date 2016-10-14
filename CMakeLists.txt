cmake_minimum_required(VERSION 2.8)

project("counting")

set(GCC_COMPILE_FLAGS "--std=c++11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")


set(SOURCES 
    Counter.cpp
  )

set(HEADERS 
    Counter.h
  )

include_directories(/home/roar/root/root-6.06.06/include)
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/lib/Timer)
link_directories(/home/roar/root/root-6.06.06/lib)

add_executable(count main.cpp ${SOURCES} ${HEADERS})
target_link_libraries(count dl Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Gui pthread m)
