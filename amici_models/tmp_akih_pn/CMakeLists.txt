cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0060)
  cmake_policy(SET CMP0060 NEW)
endif(POLICY CMP0060)
if(POLICY CMP0065)
  cmake_policy(SET CMP0065 NEW)
endif(POLICY CMP0065)
if(POLICY CMP0074)
  # Use package_ROOT environment variables
  cmake_policy(SET CMP0074 NEW)
endif(POLICY CMP0074)

project(tmp_akih_pn)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCXXCompilerFlag)
set(MY_CXX_FLAGS -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable)
foreach(FLAG ${MY_CXX_FLAGS})
    unset(CUR_FLAG_SUPPORTED CACHE)
    CHECK_CXX_COMPILER_FLAG(${FLAG} CUR_FLAG_SUPPORTED)
    if(${CUR_FLAG_SUPPORTED})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
endforeach(FLAG)

find_package(Amici 0.11.6 HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build)
message(STATUS "Found AMICI ${Amici_DIR}")

set(MODEL_DIR ${CMAKE_CURRENT_LIST_DIR})

set(SRC_LIST_LIB tmp_akih_pn_J.cpp 
tmp_akih_pn_JB.cpp 
tmp_akih_pn_JDiag.cpp 
tmp_akih_pn_JSparse.cpp 
tmp_akih_pn_JSparseB.cpp 
tmp_akih_pn_Jy.cpp 
tmp_akih_pn_dJydsigmay.cpp 
tmp_akih_pn_dJydy.cpp 
tmp_akih_pn_dwdp.cpp 
tmp_akih_pn_dwdx.cpp 
tmp_akih_pn_spline_constructors.cpp 
tmp_akih_pn_dspline_valuesdp.cpp 
tmp_akih_pn_dspline_slopesdp.cpp 
tmp_akih_pn_dxdotdw.cpp 
tmp_akih_pn_dxdotdp_explicit.cpp 
tmp_akih_pn_dydx.cpp 
tmp_akih_pn_dydp.cpp 
tmp_akih_pn_dsigmaydp.cpp 
tmp_akih_pn_sigmay.cpp 
tmp_akih_pn_w.cpp 
tmp_akih_pn_x0.cpp 
tmp_akih_pn_x0_fixedParameters.cpp 
tmp_akih_pn_sx0.cpp 
tmp_akih_pn_sx0_fixedParameters.cpp 
tmp_akih_pn_xdot.cpp 
tmp_akih_pn_y.cpp 
tmp_akih_pn_x_rdata.cpp 
tmp_akih_pn_total_cl.cpp 
tmp_akih_pn_x_solver.cpp 
tmp_akih_pn_JSparse_colptrs.cpp
tmp_akih_pn_JSparse_rowvals.cpp 
tmp_akih_pn_JSparseB_colptrs.cpp
tmp_akih_pn_JSparseB_rowvals.cpp 
tmp_akih_pn_dJydy_colptrs.cpp
tmp_akih_pn_dJydy_rowvals.cpp 
tmp_akih_pn_dwdp_colptrs.cpp
tmp_akih_pn_dwdp_rowvals.cpp 
tmp_akih_pn_dwdx_colptrs.cpp
tmp_akih_pn_dwdx_rowvals.cpp 
tmp_akih_pn_dxdotdw_colptrs.cpp
tmp_akih_pn_dxdotdw_rowvals.cpp 
tmp_akih_pn_dxdotdp_explicit_colptrs.cpp
tmp_akih_pn_dxdotdp_explicit_rowvals.cpp 
tmp_akih_pn_dxdotdp_implicit_colptrs.cpp
tmp_akih_pn_dxdotdp_implicit_rowvals.cpp 
tmp_akih_pn.cpp
${MODEL_DIR}/wrapfunctions.cpp
)

add_library(${PROJECT_NAME} ${SRC_LIST_LIB})
add_library(model ALIAS ${PROJECT_NAME})

# This option can be helpful when using the Intel compiler and compilation of
# wrapfunctions.cpp fails due to insufficient memory.
option(ENABLE_WRAPFUNCTIONS_OPTIMIZATIONS "Enable compiler optimizations for wrapfunctions.cpp?" ON)
if(NOT ENABLE_WRAPFUNCTIONS_OPTIMIZATIONS)
    set_source_files_properties(wrapfunctions.cpp PROPERTIES COMPILE_FLAGS -O0)
endif()

target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(${PROJECT_NAME}
    PUBLIC Upstream::amici
)

set(SRC_LIST_EXE main.cpp)

add_executable(simulate_${PROJECT_NAME} ${SRC_LIST_EXE})

target_link_libraries(simulate_${PROJECT_NAME} ${PROJECT_NAME})

if($ENV{ENABLE_GCOV_COVERAGE})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 --coverage")
endif()

## SWIG
option(ENABLE_SWIG "Build swig/python library?" ON)
if(ENABLE_SWIG)
    if(NOT(${CMAKE_VERSION} VERSION_LESS 3.8))
        add_subdirectory(swig)
    else()
        message(WARNING "Unable to build SWIG interface, upgrade CMake to >=3.8.")
    endif()
endif()


# <Export cmake configuration>
include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}Targets
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
export(EXPORT ${PROJECT_NAME}Targets FILE ${PROJECT_NAME}Config.cmake
    NAMESPACE Upstream::
    )
# </Export cmake configuration>

