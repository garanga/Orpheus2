cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

find_package( Eigen3 REQUIRED )
include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCE_ ${CMAKE_CURRENT_SOURCE_DIR}/Mdb.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Model.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Point.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Line.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Body.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Part.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/MaterialLib/Material.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/MaterialLib/Isotropic.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ElementLib/ElementType.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ElementLib/P4.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Node.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Element.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Step.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/InitialStep.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/StaticStep.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Constraint.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/DisplacementConstraint.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Load.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/ConcentratedLoad.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Job.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Solver.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/LinearSolver.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/Odb.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/OdbStep.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/OdbFrame.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/FieldOutput.cpp
)

#if (TRIANGLE)
#    list (APPEND SOURCE_ ${CMAKE_CURRENT_SOURCE_DIR}/triangle.o)
#endif ()

set (CMAKE_VERBOSE_MAKEFILE ON)

if (LIBS)
    
    if (STATIC_VERSION)
        set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/static)
        add_library(Orpheus2_ STATIC ${SOURCE_})
        
        add_executable (${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/examples/main.cpp)
        target_link_libraries(${PROJECT_NAME} Orpheus2_)
        
        
    else ()
        set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/shared)
        add_library(Orpheus2_ SHARED ${SOURCE_})
#        add_executable (${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/examples/main.cpp)
#        target_link_libraries(${PROJECT_NAME} Orpheus2_)
    endif ()
    
else()
    
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY, ${PROJECT_SOURCE_DIR}/bin)
    
    list(APPEND SOURCE_ ${PROJECT_SOURCE_DIR}/examples/main.cpp)
    
    add_executable (${PROJECT_NAME} ${SOURCE_})
    
endif ()











#
#
#
#
#
#
#find_package( Eigen3 REQUIRED )
#
#
#
#file(GLOB CPPS "*.cpp")
#
#message(${PROJECT_SOURCE_DIR}/${${PROJECT_NAME}_INCLUDE_DIR}/MaterialLib)
#
#
#include_directories(${PROJECT_SOURCE_DIR}/${${PROJECT_NAME}_INCLUDE_DIR})
#include_directories(${PROJECT_SOURCE_DIR}/${${PROJECT_NAME}_INCLUDE_DIR}/ElementLib)
#include_directories(${PROJECT_SOURCE_DIR}/${${PROJECT_NAME}_INCLUDE_DIR}/MaterialLib)
#include_directories(${EIGEN3_INCLUDE_DIR})
#
#message(STATUS ${INCLUDE_DIRECTORIES})
#message(${EIGEN3_INCLUDE_DIR})
#
#add_executable (Orpheus2 ${CPPS})


#set(${PROJECT_NAME}_VERSION_MAJOR 1)
#set(${PROJECT_NAME}_VERSION_MINOR 0)
#set(${PROJECT_NAME}_VERSION_PATCH 0)
#set(${PROJECT_NAME}_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

#message(STATUS ${${PROJECT}_VERSION})


#set (MAIN_SOURCES main.cpp)



#add_executable(Orpheus2 main.cpp)


#message(${PROJECT_NAME})

#message(${PROJECT_NAME})


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