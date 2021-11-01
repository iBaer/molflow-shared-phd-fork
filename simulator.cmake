################# CMake Project ############################
# The main options of project                              #
############################################################

set(PROJECT_NAME shared_simulator)
project(${PROJECT_NAME} CXX)

# Add library to build.
add_library(${PROJECT_NAME} STATIC
        ${SIM_FILES}
        ${HELPER_FILES}
        )

target_include_directories(${PROJECT_NAME} PRIVATE
        ${HEADER_DIR_1}
        ${HEADER_DIR_2}
        ${HEADER_DIR_3}
        ${HEADER_DIR_4}
        )

target_include_directories(${PROJECT_NAME} PUBLIC
        ${HEADER_DIR_5}
        ${HEADER_DIR_6}
        ${HEADER_DIR_7}
        )

include(SetOpenMP.cmake)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Detected OpenMP version: ${OpenMP_CXX_VERSION}")
    #target_include_directories(${PROJECT_NAME} PRIVATE OpenMP_INCLUDE_DIR)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

if(NOT MSVC)
    find_package(GSL REQUIRED)
    target_include_directories(${PROJECT_NAME} PUBLIC ${GSL_INCLUDE_DIRS})

    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)

    target_link_libraries(${PROJECT_NAME} PUBLIC ${GSL_LIBRARIES})
    target_link_libraries(${PROJECT_NAME} PUBLIC Threads::Threads)

    #for shared memory
    find_library(LIBRT rt)
    if(LIBRT)
        target_link_libraries(${PROJECT_NAME} PUBLIC ${LIBRT})
    endif()

    #[[# Your-external "mylib", add GLOBAL if the imported library is located in directories above the current.
    if (NOT TARGET libzip)
        add_library( libzip STATIC IMPORTED GLOBAL)
    endif()
    # You can define two import-locations: one for debug and one for release.
    get_filename_component(ABS_LINK_DIR_1 "${LINK_DIR_1}" REALPATH)

    # other static libraries need to be built with similar settings for clang (here: libc++)
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set_target_properties( libzip PROPERTIES IMPORTED_LOCATION ${ABS_LINK_DIR_1}/libzip_clang.a )
    else()
        set_target_properties( libzip PROPERTIES IMPORTED_LOCATION ${ABS_LINK_DIR_1}/libzip_gcc.a )
    endif()]]
    target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ziplib/Source)
    target_link_libraries(${PROJECT_NAME} PUBLIC ziplib)
    #target_include_directories(${PROJECT_NAME} INTERFACE fmtlib_src)


    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_link_libraries(${PROJECT_NAME} PUBLIC c++fs)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
        #don´t add anything for filesystem
    else()
        target_link_libraries(${PROJECT_NAME} PUBLIC stdc++fs)
    endif()

endif(NOT MSVC)

target_link_libraries(${PROJECT_NAME} PUBLIC fmtlib_src) # header include
target_link_libraries(${PROJECT_NAME} PUBLIC fmt)
######################### Flags ############################
# Defines Flags for Windows and Linux                      #
############################################################

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

# Preprocessor definitions
#[[
if(CMAKE_BUILD_TYPE MATCHES Debug|RelWithDebInfo)
    target_compile_definitions(${PROJECT_NAME} PRIVATE)
    if(MSVC)
        target_compile_options(${PROJECT_NAME} PRIVATE /MDd /Od /EHsc)
        #for .pdb debugging files
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi")
        set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /DEBUG:FULL /OPT:REF /OPT:ICF")
    endif()
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_definitions(${PROJECT_NAME} PRIVATE)
    if(MSVC)
        target_compile_options(${PROJECT_NAME} PRIVATE /GL /Oi /Gy /EHsc)
    endif()
endif()
]]