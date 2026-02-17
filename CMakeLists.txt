cmake_minimum_required(VERSION 3.20)
project(bmercator LANGUAGES Fortran CXX)

option(BMERCATOR_ENABLE_OPENMP "Enable OpenMP acceleration when available" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -c -fpic -Wall -O3")

include_directories(include)

add_library(trapecis_Sd STATIC include/trapecis_Sd.f)
# add_library(trapecis_Sd_expected STATIC include/trapecis_Sd_expected.f)

set(source_dir "${PROJECT_SOURCE_DIR}/src/")
set(source_files "${source_dir}/embeddingSD_unix.cpp")

set(include_dir "${PROJECT_SOURCE_DIR}/include/")
file(GLOB include_files "${include_dir}/*.hpp")

add_executable(bmercator ${source_files} ${include_files})

target_compile_options(bmercator PRIVATE
  $<$<CONFIG:Debug>:-Wall;-Wextra;-Wshadow;-pedantic>
  $<$<CONFIG:Release>:-O2>
)

find_package(Threads REQUIRED)
target_link_libraries(bmercator PRIVATE trapecis_Sd Threads::Threads)

if(BMERCATOR_ENABLE_OPENMP)
  if(APPLE AND NOT DEFINED OpenMP_ROOT)
    foreach(_libomp_prefix /opt/homebrew/opt/libomp /usr/local/opt/libomp)
      if(EXISTS "${_libomp_prefix}")
        set(OpenMP_ROOT "${_libomp_prefix}")
        message(STATUS "Trying Homebrew libomp at ${OpenMP_ROOT}")
        break()
      endif()
    endforeach()
  endif()

  find_package(OpenMP QUIET COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(bmercator PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled for C++ target.")
  else()
    message(STATUS "OpenMP not found. Building without OpenMP.")
  endif()
else()
  message(STATUS "OpenMP disabled by BMERCATOR_ENABLE_OPENMP=OFF.")
endif()
