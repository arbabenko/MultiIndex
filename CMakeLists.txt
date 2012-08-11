##########################################################################################
# CMake build script for NearestSearch.
#
##########################################################################################

#let all libraries be static, not shared
OPTION(BUILD_SHARED_LIBS "Build shared libraries (DLLs)." OFF)

##########################################################################################
# lets start describing our project.
project (NearestSearch CXX C)
cmake_minimum_required(VERSION 2.6)

IF (UNIX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
ENDIF (UNIX)

IF (NOT DEFINED CMAKE_CURRENT_LIST_DIR)
    SET(CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF (NOT DEFINED CMAKE_CURRENT_LIST_DIR)

MESSAGE("current dir: ${CMAKE_CURRENT_LIST_DIR}")
SET (Source_Path ${CMAKE_CURRENT_LIST_DIR})

MESSAGE ("PROJECT_BINARY_DIR " ${PROJECT_BINARY_DIR})

# ADD BOOSTDIR
SET(BOOST_DIR C:/Users/arbabenko/Soft/boost_1_47_0)
SET(BOOST_LIB C:/Users/arbabenko/Soft/boost_1_47_0/lib/x64/lib)

# ADD BLAS
SET(BLAS_DIR "C:/Program Files (x86)/Intel/Composer XE 2011 SP1/mkl/include")
SET(BLAS_LIB "C:/Program Files (x86)/Intel/Composer XE 2011 SP1/mkl/lib/intel64"
             "C:/Program Files (x86)/Intel/Composer XE 2011 SP1/compiler/lib/intel64")

##########################################################################################
INCLUDE (${Source_Path}/NearestSearch.cmake NO_POLICY_SCOPE)
INCLUDE_DIRECTORIES(${BOOST_DIR})
INCLUDE_DIRECTORIES(${BLAS_DIR})
LINK_DIRECTORIES(${BOOST_LIB})
LINK_DIRECTORIES(${BLAS_LIB})
INCLUDE_DIRECTORIES(${Source_Path})
ADD_LIBRARY (nearest_search_lib STATIC ${NEAREST_ALL_CC})

##########################################################################################
INCLUDE (${Source_Path}/IndexerLauncher.cmake NO_POLICY_SCOPE)
INCLUDE_DIRECTORIES(${BOOST_DIR})
INCLUDE_DIRECTORIES(${BLAS_DIR})
LINK_DIRECTORIES(${BOOST_LIB})
LINK_DIRECTORIES(${BLAS_LIB})
INCLUDE_DIRECTORIES(${Source_Path})
ADD_EXECUTABLE (indexer_launcher ${IndexerLauncher})
TARGET_LINK_LIBRARIES (indexer_launcher nearest_search_lib)
TARGET_LINK_LIBRARIES (indexer_launcher  mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib)
  
##########################################################################################
INCLUDE (${Source_Path}/SearcherTester.cmake NO_POLICY_SCOPE)
INCLUDE_DIRECTORIES(${BOOST_DIR})
INCLUDE_DIRECTORIES(${BLAS_DIR})
LINK_DIRECTORIES(${BOOST_LIB})
LINK_DIRECTORIES(${BLAS_LIB})
INCLUDE_DIRECTORIES(${Source_Path})
ADD_EXECUTABLE (searcher_tester ${SearchTester})
TARGET_LINK_LIBRARIES (searcher_tester nearest_search_lib)
TARGET_LINK_LIBRARIES (searcher_tester mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib)


