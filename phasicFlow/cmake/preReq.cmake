# TBB is only the libstdc++ parallel-STL backend: use it if found, never install it
if(pFlow_STD_Parallel_Alg)
    find_library(TBB_LIBRARY NAMES tbb HINTS ENV LIBRARY_PATH ENV LD_LIBRARY_PATH)
    if(TBB_LIBRARY)
        message(STATUS "Found TBB: ${TBB_LIBRARY}")
    else()
        message(WARNING
            "TBB not found -> building with pFlow_STD_Parallel_Alg=OFF (serial std:: "
            "algorithms). This is the expected configuration on the cluster. "
            "To force it on, load a TBB module or set LIBRARY_PATH/CPATH.")
        set(pFlow_STD_Parallel_Alg OFF CACHE BOOL "Use TBB std parallel algorithms" FORCE)
    endif()
endif()


# Kokkos source: -DKokkos_Source_DIR, the thirdParty/kokkos submodule (4.4.01),
# $Kokkos_DIR, $HOME/Kokkos/kokkos, else FetchContent
if(DEFINED Kokkos_Source_DIR AND EXISTS "${Kokkos_Source_DIR}/CMakeLists.txt")
    message(STATUS "Kokkos: using -DKokkos_Source_DIR")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/thirdParty/kokkos/CMakeLists.txt")
    set(Kokkos_Source_DIR "${CMAKE_CURRENT_SOURCE_DIR}/thirdParty/kokkos")
    message(STATUS "Kokkos: using the thirdParty/kokkos submodule")
elseif(DEFINED ENV{Kokkos_DIR} AND EXISTS "$ENV{Kokkos_DIR}/CMakeLists.txt")
    set(Kokkos_Source_DIR "$ENV{Kokkos_DIR}")
    message(STATUS "Kokkos: using $Kokkos_DIR from the environment")
elseif(EXISTS "$ENV{HOME}/Kokkos/kokkos/CMakeLists.txt")
    set(Kokkos_Source_DIR "$ENV{HOME}/Kokkos/kokkos")
    message(STATUS "Kokkos: using $HOME/Kokkos/kokkos")
else()
    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.30")
        cmake_policy(SET CMP0169 OLD)
    endif()

    include(FetchContent)
    FetchContent_Declare(
        kokkos
        GIT_REPOSITORY https://github.com/kokkos/kokkos.git
        GIT_TAG 4.4.01
    )

    FetchContent_GetProperties(kokkos)
    if(NOT kokkos_POPULATED)
        message(STATUS "No local Kokkos found. Downloading Kokkos 4.4.01 ...")
        FetchContent_Populate(kokkos)
        set(Kokkos_Source_DIR ${kokkos_SOURCE_DIR})
    endif()
endif()

message(STATUS "Kokkos source directory is ${Kokkos_Source_DIR}")
add_subdirectory(${Kokkos_Source_DIR} ./kokkos)
