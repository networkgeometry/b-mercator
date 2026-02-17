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

if(APPLE)
  # Ensure dyld can resolve Fortran/OpenMP runtimes even when they use @rpath IDs.
  set(_bmercator_rpaths
    /opt/homebrew/opt/gcc/lib/gcc/current
    /usr/local/opt/gcc/lib/gcc/current
    /opt/homebrew/opt/libomp/lib
    /usr/local/opt/libomp/lib
  )

  if(CMAKE_Fortran_COMPILER)
    execute_process(
      COMMAND "${CMAKE_Fortran_COMPILER}" -print-file-name=libgfortran.dylib
      OUTPUT_VARIABLE _bmercator_libgfortran
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    if(_bmercator_libgfortran AND NOT _bmercator_libgfortran STREQUAL "libgfortran.dylib")
      get_filename_component(_bmercator_libgfortran_dir "${_bmercator_libgfortran}" DIRECTORY)
      list(APPEND _bmercator_rpaths "${_bmercator_libgfortran_dir}")
    endif()

    execute_process(
      COMMAND "${CMAKE_Fortran_COMPILER}" -print-file-name=libquadmath.dylib
      OUTPUT_VARIABLE _bmercator_libquadmath
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    if(_bmercator_libquadmath AND NOT _bmercator_libquadmath STREQUAL "libquadmath.dylib")
      get_filename_component(_bmercator_libquadmath_dir "${_bmercator_libquadmath}" DIRECTORY)
      list(APPEND _bmercator_rpaths "${_bmercator_libquadmath_dir}")
    endif()
  endif()

  if(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES)
    list(APPEND _bmercator_rpaths ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
  endif()

  if(DEFINED OpenMP_ROOT AND IS_DIRECTORY "${OpenMP_ROOT}/lib")
    list(APPEND _bmercator_rpaths "${OpenMP_ROOT}/lib")
  endif()

  set(_bmercator_existing_rpaths "")
  foreach(_path IN LISTS _bmercator_rpaths)
    if(IS_DIRECTORY "${_path}")
      list(APPEND _bmercator_existing_rpaths "${_path}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES _bmercator_existing_rpaths)

  if(_bmercator_existing_rpaths)
    set_property(TARGET bmercator PROPERTY BUILD_RPATH "${_bmercator_existing_rpaths}")
    set_property(TARGET bmercator PROPERTY INSTALL_RPATH "${_bmercator_existing_rpaths}")
    set_property(TARGET bmercator PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
    message(STATUS "Configured macOS runtime rpaths: ${_bmercator_existing_rpaths}")
  endif()
endif()
