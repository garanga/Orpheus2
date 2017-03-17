







#set(SOURCE_EXE main.cpp)

#file(GLOB CPPS "*.cpp")

#set(SOURCE_LIBS )

#set(${PROJECT_NAME}_BIN_DIR bin)
#set(${PROJECT_NAME}_SRC_DIR src)
#set(${PROJECT_NAME}_INCLUDE_DIR include)

#set (CMAKE_INCLUDE_PATH ${PROJECT_SOURCE_DIR}/${${PROJECT_NAME}_INCLUDE_DIR})
#set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/${${PROJECT_NAME}_BIN_DIR})

#set (CMAKE_VERBOSE_MAKEFILE ON)
#set (CMAKE_BUILD_TYPE Debug)

#set (CMAKE_BUILD_TYPE Debug)

#if ( CMAKE_COMPILER_IS_GNUCXX )
#    set(MY_CXX_FLAGS "-std=c++0x -Wall -fmessage-length=0")
#    set(CMAKE_CXX_FLAGS_DEBUG "-O3 ${MY_CXX_FLAGS}")
#else ()
#    message (FATAL_ERROR ${ERR_NO_UNIX})
#endif ()
#
#
#
#add_subdirectory(src)
#
#







#include_directories("include")

#add_executable(Orpheus2 main.cpp)





#message(STATUS ${CMAKE_MODULE_PATH})


#option(WITH_DEBUG_MODE "Build with debug mode" ON)

#if (NOT UNIX )
#    message (FATAL_ERROR "Not Unix!")
#endif ()

#message ("Starting cmake")
#
#message(STATUS ${PROJECT} "1")
#message(STATUS "1")
#




#set ("${PROJECT}_BINARY_DIR" bin)
#set ("${PROJECT}_SOURCE_DIR" src)
#set ("${PROJECT}_LIB_DIR" lib)

#set (CMAKE_INCLUDE_PATH ${${PROJECT}_SOURCE_DIR})
#set (CMAKE_LIBRARY_PATH ${${PROJECT}_LIB_DIR})
#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/${${PROJECT}_BINARY_DIR})e
#set (CMAKE_VERBOSE_MAKEFILE ON)
#set (CMAKE_BUILD_TYPE Debug)
#
#set (ERR_NO_UNIX "Cannot build on non Unix systems")
#
#if ( WITH_DEBUG_MODE )
#     ADD_DEFINITIONS( -DMY_DEBUG_MODE=1)
#endif()
#
#if ( CMAKE_COMPILER_IS_GNUCXX )
#    set(MY_CXX_FLAGS "-Wall -std=c++0x -fmessage-length=0 -v -I/usr/local/include/eigen3")
#    set(CMAKE_CXX_FLAGS "-O0 ${MY_CXX_FLAGS}")
#    # я отключил настройку разных флагов для релиза и отладки. Пока что не нужно.
#    #set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
#    #set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
#else
#    message (FATAL_ERROR ${ERR_NO_UNIX})
#endif () 



