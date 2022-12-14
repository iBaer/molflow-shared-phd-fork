cmake_minimum_required(VERSION 3.12.2 FATAL_ERROR)

##################### Variables ############################
# Change if you want modify path or other values           #
############################################################

# Project
set(PROJECT_NAME fmtlib)

# Folders files
#set(CPP_DIR_1 ./molflow_wrapper)
#set(HEADER_DIR_1 ./molflow_wrapper)

################# CMake Project ############################
# The main options of project                              #
############################################################

project(${PROJECT_NAME} C CXX)

############### Files & Targets ############################
# Files of project and target to build                     #
############################################################

# Source Files

# Download and unpack fmt at configure time
configure_file(CMakeLists.txt.in fmt-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fmt-download )
if(result)
    message(FATAL_ERROR "CMake step for nfd failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fmt-download )
if(result)
    message(FATAL_ERROR "Build step for fmt failed: ${result}")
endif()

# Add fmt directly to our build. This defines
# the fmt target.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/fmt-src
        ${CMAKE_CURRENT_BINARY_DIR}/fmt-build
        EXCLUDE_FROM_ALL
        )

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
#[[if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${fmtlib_BINARY_DIR}/fmt-src/include")
endif()]]
include_directories(SYSTEM ${CMAKE_CURRENT_BINARY_DIR}/fmt-src/include)

# Get the current working branch
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest commit hash
execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

#add_subdirectory(nativefiledialog-extended)

#add_library(${PROJECT_NAME} STATIC ${WRAPPER_SRC})

add_library(${PROJECT_NAME}_src INTERFACE)
#add_library(${PROJECT_NAME}_src STATIC ${fmtlib_SOURCE_DIR})
target_include_directories(${PROJECT_NAME}_src INTERFACE include/)
target_include_directories(${PROJECT_NAME}_src INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/fmt-src/include/)
#target_link_libraries(${PROJECT_NAME} PUBLIC fmt)
#add_library(${PROJECT_NAME}_src STATIC ${fmt_SRC})
#add_dependencies(${PROJECT_NAME} fmt)

#target_include_directories(${PROJECT_NAME} PRIVATE ${HEADER_DIR_1})
#target_link_libraries(${PROJECT_NAME} PRIVATE fmt)
#target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_src)

######################### Flags ############################
# Defines Flags for Windows and Linux                      #
############################################################

#[[IF (WIN32)

ELSEIF(APPLE)

    target_link_libraries(${PROJECT_NAME}_src
            "-framework AppKit"
            )
ELSE()
    # Use the package PkgConfig to detect GTK+ headers/library files
    FIND_PACKAGE(PkgConfig REQUIRED)
    PKG_CHECK_MODULES(GTK3 REQUIRED gtk+-3.0)

    # Setup CMake to use GTK+, tell the compiler where to look for headers
    # and to the linker where to look for libraries
    INCLUDE_DIRECTORIES(${GTK3_INCLUDE_DIRS})
    LINK_DIRECTORIES(${GTK3_LIBRARY_DIRS})

    # Add other flags to the compiler
    ADD_DEFINITIONS(${GTK3_CFLAGS_OTHER})
    target_link_libraries(${PROJECT_NAME}_src ${GTK3_LIBRARIES})
ENDIF()]]

# Preprocessor definitions
#[[
if(CMAKE_BUILD_TYPE MATCHES Debug|RelWithDebInfo)
    target_compile_definitions(${PROJECT_NAME} PRIVATE)
    if(MSVC)
        target_compile_options(${PROJECT_NAME} PRIVATE  /W3 /MDd /Od /EHsc)
        #for .pdb debugging files
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi")
        set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /DEBUG:FULL /OPT:REF /OPT:ICF")
    endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_definitions(${PROJECT_NAME} PRIVATE)
    if(MSVC)
        target_compile_options(${PROJECT_NAME} PRIVATE  /W3 /GL /Oi /Gy /EHsc)
    endif()
endif()]]
