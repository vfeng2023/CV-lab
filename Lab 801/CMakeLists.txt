cmake_minimum_required(VERSION 2.8)
project( l081 )
find_package( OpenCV REQUIRED )
add_executable( l081 l081.cpp )
target_link_libraries( l081 ${OpenCV_LIBS} )

