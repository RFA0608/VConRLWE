cmake_minimum_required(VERSION 3.16)
project(main)

find_package(SEAL REQUIRED)

find_path(GMPXX_INCLUDE_DIR NAMES gmpxx.h)
find_path(GMP_INCLUDE_DIR NAMES gmp.h)

find_library(GMPXX_LIB NAMES gmpxx)
find_library(GMP_LIB NAMES gmp)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_LIB GMPXX_LIB GMPXX_INCLUDE_DIR)

if(SEAL_FOUND AND GMP_FOUND)
	include_directories(${SEAL_INCLUDE_DIRS})
	
	add_executable(main main.cpp)
    target_include_directories(main PUBLIC ${GMPXX_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
    target_link_libraries(main PUBLIC 
        SEAL::seal 
        ${GMPXX_LIB} 
        ${GMP_LIB}
    )
else()
	message(FATAL_ERROR "Microsoft SEAL library not found.")
endif()