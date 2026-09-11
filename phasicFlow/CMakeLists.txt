cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# set the project name and version
project(phasicFlow VERSION 1.0 )

set(CMAKE_CXX_STANDARD 20 CACHE STRING "" FORCE)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# install into the OpenFOAM user platforms tree ($FOAM_USER_APPBIN, $FOAM_USER_LIBBIN)
if(NOT DEFINED pFlow_INSTALL_PREFIX)
	if(DEFINED ENV{WM_PROJECT_USER_DIR} AND DEFINED ENV{WM_OPTIONS})
		set(pFlow_INSTALL_PREFIX "$ENV{WM_PROJECT_USER_DIR}/platforms/$ENV{WM_OPTIONS}")
	else()
		message(FATAL_ERROR
			"WM_PROJECT_USER_DIR / WM_OPTIONS are not set. Source the OpenFOAM "
			"environment first (e.g. 'source /usr/lib/openfoam/openfoam2512/etc/bashrc'), "
			"or configure with -DpFlow_INSTALL_PREFIX=<path>.")
	endif()
endif()
set(CMAKE_INSTALL_PREFIX "${pFlow_INSTALL_PREFIX}" CACHE PATH "Install path of phasicFlow" FORCE)

# one lib directory (GNUInstallDirs would send Kokkos to lib64 on RHEL)
set(CMAKE_INSTALL_LIBDIR     lib     CACHE PATH "" FORCE)
set(CMAKE_INSTALL_BINDIR     bin     CACHE PATH "" FORCE)
set(CMAKE_INSTALL_INCLUDEDIR include CACHE PATH "" FORCE)

# relative RPATH, so the installed tree resolves itself without LD_LIBRARY_PATH
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "\$ORIGIN;\$ORIGIN/../lib")
set(CMAKE_BUILD_TYPE Release CACHE STRING "build type")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build using shared libraries" FORCE)
mark_as_advanced(FORCE var BUILD_SHARED_LIBS)

message(STATUS "Install prefix is:" ${CMAKE_INSTALL_PREFIX})

include(cmake/globals.cmake)

option(pFlow_STD_Parallel_Alg "Use TTB std parallel algorithms" ON)
option(pFlow_Build_Serial "Build phasicFlow and backends for serial execution" OFF)
option(pFlow_Build_OpenMP "Build phasicFlow and backends for OpenMP execution" OFF)
option(pFlow_Build_Cuda   "Build phasicFlow and backends for Cuda execution" OFF)
option(pFlow_Build_Double "Build phasicFlow with double precision floating-oint variables" ON)
option(pFlow_Build_MPI    "Build for MPI parallelization. This will enable multi-gpu run, CPU run on clusters (distributed memory machine). Use this combination Cuda+MPI, OpenMP + MPI or Serial+MPI " OFF)

# the Kokkos backend must be chosen before preReq.cmake adds the kokkos subdirectory
if(pFlow_Build_Serial)
	set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Serial execution" FORCE)
	set(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "OpenMP execution" FORCE)
	set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "Cuda execution" FORCE)
	set(Kokkos_ENABLE_CUDA_LAMBDA OFF CACHE BOOL "Cuda execution" FORCE)
	set(Kokkos_ENABLE_CUDA_CONSTEXPR OFF CACHE BOOL "Enable constexpr on cuda code" FORCE)
elseif(pFlow_Build_OpenMP )
	set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Serial execution" FORCE)
	set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "OpenMP execution" FORCE)
	set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "Cuda execution" FORCE)
	set(Kokkos_ENABLE_CUDA_LAMBDA OFF CACHE BOOL "Cuda execution" FORCE)
	set(Kokkos_DEFAULT_HOST_PARALLEL_EXECUTION_SPACE Serial CACHE STRING "" FORCE)
	set(Kokkos_DEFAULT_DEVICE_PARALLEL_EXECUTION_SPACE OpenMP CACHE STRING "" FORCE)
	set(Kokkos_ENABLE_CUDA_CONSTEXPR OFF CACHE BOOL "Enable constexpr on cuda code" FORCE)
elseif(pFlow_Build_Cuda)
	set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Serial execution" FORCE)
	set(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "OpenMP execution" FORCE)
	set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Cuda execution" FORCE)
	set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "Cuda execution" FORCE)
	set(Kokkos_ENABLE_CUDA_CONSTEXPR ON CACHE BOOL "Enable constexpr on cuda code" FORCE)
endif()

#for installing the required packages
include(cmake/preReq.cmake)

if(pFlow_Build_MPI)
	find_package(MPI REQUIRED)
endif()

include(cmake/makeLibraryGlobals.cmake)
include(cmake/makeExecutableGlobals.cmake)

configure_file(phasicFlowConfig.H.in phasicFlowConfig.H)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#add a global include directory 
include_directories(src/setHelpers src/demComponent "${PROJECT_BINARY_DIR}")

add_subdirectory(src)

add_subdirectory(solvers)

add_subdirectory(utilities)

add_subdirectory(DEMSystems)

install(FILES "${PROJECT_BINARY_DIR}/phasicFlowConfig.H"
  DESTINATION include)

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${phasicFlow_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${phasicFlow_VERSION_MINOR}")
include(CPack)

