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

project(vaccination_piecewise)

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

set(SRC_LIST_LIB vaccination_piecewise_J.cpp 
vaccination_piecewise_JB.cpp 
vaccination_piecewise_JDiag.cpp 
vaccination_piecewise_JSparse.cpp 
vaccination_piecewise_JSparseB.cpp 
vaccination_piecewise_Jy.cpp 
vaccination_piecewise_dJydsigmay.cpp 
vaccination_piecewise_dJydy.cpp 
vaccination_piecewise_dwdp.cpp 
vaccination_piecewise_dwdx.cpp 
vaccination_piecewise_spline_constructors.cpp 
vaccination_piecewise_dspline_valuesdp.cpp 
vaccination_piecewise_dspline_slopesdp.cpp 
vaccination_piecewise_dxdotdw.cpp 
vaccination_piecewise_dxdotdp_explicit.cpp 
vaccination_piecewise_dydx.cpp 
vaccination_piecewise_dydp.cpp 
vaccination_piecewise_dsigmaydp.cpp 
vaccination_piecewise_sigmay.cpp 
vaccination_piecewise_w.cpp 
vaccination_piecewise_x0.cpp 
vaccination_piecewise_x0_fixedParameters.cpp 
vaccination_piecewise_sx0.cpp 
vaccination_piecewise_sx0_fixedParameters.cpp 
vaccination_piecewise_xdot.cpp 
vaccination_piecewise_y.cpp 
vaccination_piecewise_x_rdata.cpp 
vaccination_piecewise_total_cl.cpp 
vaccination_piecewise_x_solver.cpp 
vaccination_piecewise_JSparse_colptrs.cpp
vaccination_piecewise_JSparse_rowvals.cpp 
vaccination_piecewise_JSparseB_colptrs.cpp
vaccination_piecewise_JSparseB_rowvals.cpp 
vaccination_piecewise_dJydy_colptrs.cpp
vaccination_piecewise_dJydy_rowvals.cpp 
vaccination_piecewise_dwdp_colptrs.cpp
vaccination_piecewise_dwdp_rowvals.cpp 
vaccination_piecewise_dwdx_colptrs.cpp
vaccination_piecewise_dwdx_rowvals.cpp 
vaccination_piecewise_dxdotdw_colptrs.cpp
vaccination_piecewise_dxdotdw_rowvals.cpp 
vaccination_piecewise_dxdotdp_explicit_colptrs.cpp
vaccination_piecewise_dxdotdp_explicit_rowvals.cpp 
vaccination_piecewise_dxdotdp_implicit_colptrs.cpp
vaccination_piecewise_dxdotdp_implicit_rowvals.cpp 
vaccination_piecewise.cpp
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

