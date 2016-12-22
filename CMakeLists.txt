cmake_minimum_required(VERSION 2.8)

project("counting")

set(GCC_COMPILE_FLAGS "-Wall --std=c++11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")


set(SOURCES 
    src/Counter.cpp
  )

set(HEADERS 
    src/Counter.h
  )

include_directories(/home/roar/div_progs/root/include)
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/lib/Timer)
link_directories(/home/roar/div_progs/root/lib)

add_executable(count src/main.cpp ${SOURCES} ${HEADERS})
target_link_libraries(count Core RIO Hist)
# (count dl Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread Gui pthread m)
