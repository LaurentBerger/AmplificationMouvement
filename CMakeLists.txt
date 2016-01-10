cmake_minimum_required(VERSION 2.8)

PROJECT (AmplificationMouvement)

find_package(OpenCV  REQUIRED)


file(GLOB AmplificationMouvement_SRCS
    "*.h"
    "*.cpp")
ADD_EXECUTABLE (AmplificationMouvement ${AmplificationMouvement_SRCS})



if (OpenCV_FOUND)
	include_directories( ${OpenCV_INCLUDE_DIRS} )
	target_link_libraries( AmplificationMouvement ${OpenCV_LIBS} )
else (OpenCV_FOUND)
message("PB->OPENCV = ${OpenCV_INCLUDE_DIRS}")
endif (OpenCV_FOUND)

