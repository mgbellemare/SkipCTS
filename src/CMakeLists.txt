PROJECT(SkipCTS)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

# remove deprecated warnings from MSVC
IF (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    ADD_DEFINITIONS(/D_CRT_SECURE_NO_WARNINGS)
	ADD_DEFINITIONS(/D_SCL_SECURE_NO_WARNINGS)
ENDIF()

# x86_64 Linux cannot use static libs inside of a dynamic lib.
IF (NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    SET(Boost_USE_STATIC_LIBS ON)
ENDIF()

# Compiler flags
MESSAGE(STATUS "Test GCC")
IF (CMAKE_COMPILER_IS_GNUCXX)
    MESSAGE(STATUS "Add the optimization option -O2")
    SET(CMAKE_CXX_FLAGS "-O2")        ## Optimize
ENDIF()    

# The NO_BOOST_SYSTEM_PATHS flag enables us to specify an alternate path for the root of our
# Boost installation root. This is done via the option BOOST_ROOT.
IF (NO_BOOST_SYSTEM_PATHS)
    SET(Boost_NO_SYSTEM_PATHS ON)
ELSE()
    SET(Boost_NO_SYSTEM_PATHS OFF)
ENDIF()

SET(Boost_USE_STATIC_RUNTIME OFF)
# We require a number of Boost packages for I/O and command-line options.
FIND_PACKAGE(Boost 1.55.0 REQUIRED COMPONENTS program_options filesystem system)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
MESSAGE(STATUS "** Boost Include dirs: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "** Boost Library dirs: ${Boost_LIBRARY_DIRS}")
MESSAGE(STATUS "** Boost Libraries: ${Boost_LIBRARIES}")

FILE(GLOB cppfiles *.cpp *.hpp *.h)
ADD_EXECUTABLE(skipcts ${cppfiles})

TARGET_LINK_LIBRARIES(skipcts ${Boost_LIBRARIES}) 

