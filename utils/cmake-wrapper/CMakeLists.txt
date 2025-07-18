# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss Contact
# information: icon-model.org See AUTHORS.TXT for a list of authors See
# LICENSES/ for license information SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

# ~~~
# This is a wrapper file that simplifies the integration of CMake-based packages
# into the ICON build system. The path to the source directory of a wrapped
# package must be provided as the ICON_DEP_SOURCE_DIR cache variable. The
# wrapper supports the following features:
#   1. Generates a .gitignore file in the build directory to ignore all files in
#      it, which is useful when there is no control over the .gitignore file of
#      a git submodule and it is not possible to ignore artefacts produced by
#      the wrapper. See the option ICON_GENERATE_GITIGNORE.
#   2. Checks whether the compilers in use belong to the NEC family (CMake does
#      not support NEC compilers yet and usually mistakes them for GNU) and
#      applies several workarouns:
#        - reruns cmake_determine_compiler_abi with adjusted arguments to get
#          the correct values for the CMAKE_<LANG>_IMPLICIT_LINK_DIRECTORIES and
#          CMAKE_<LANG>_IMPLICIT_LINK_LIBRARIES variables;
#        - sets CMAKE_Fortran_COMPILE_OPTIONS_PREPROCESS_ON to the NEC-specific
#          value.
#      See the cache variable ICON_NEC_WORKAROUND_LANGUAGES, which must be set
#      to the list of languages that the workaround should be applied for.
#   3. Provides the ICON build system with the extra linker arguments required
#      for proper linking to the targets of the wrapped package. See the cache
#      variable ICON_DEP_LIBRARIES, which must be set to the list of targets
#      produced by the package that ICON depends on.
#   4. Ensures that the 'test' target exists even if the test suite of the
#      wrapped package is disabled. This way, the build system if ICON does not
#      have to track what packages are configured with ENABLE_TESTING=OFF and
#      can safely call `make test` regardless.
# ~~~

# Do not impose any additional version constraints:
cmake_minimum_required(
  VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION})

# Do not impose any additional language constraints:
project(IconWrapper LANGUAGES NONE)

option(ICON_GENERATE_GITIGNORE
       "Generate .gitignore file in the root build directory to ignore it" OFF)

if(ICON_GENERATE_GITIGNORE)
  set(gitignore_file "${PROJECT_BINARY_DIR}/.gitignore")
  if(NOT EXISTS "${gitignore_file}")
    file(WRITE "${gitignore_file}" "*")
  endif()
endif()

set(ICON_NEC_WORKAROUND_LANGUAGES
    ""
    CACHE
      STRING
      "List of languages to enable and help CMake recognize NEC compilers and their implicit linker flags"
)
set_property(CACHE ICON_NEC_WORKAROUND_LANGUAGES PROPERTY STRINGS
                                                          "C;CXX;Fortran")

if(ICON_NEC_WORKAROUND_LANGUAGES)
  enable_language(${ICON_NEC_WORKAROUND_LANGUAGES})
  unset(nec_languages)
  # The users might be confused with the NEC-related messages when they do not
  # use NEC compilers:
  set(CMAKE_REQUIRED_QUIET 1)
  # CMake is expected to mistake the NEC compilers for the GNU ones:
  if("C" IN_LIST ICON_NEC_WORKAROUND_LANGUAGES AND CMAKE_C_COMPILER_ID STREQUAL
                                                   "GNU")
    include(CheckSymbolExists)
    check_symbol_exists(__NEC__ "" C_is_NEC)
    if(C_is_NEC)
      set(abi_C_file "${CMAKE_ROOT}/Modules/CMakeCCompilerABI.c")
      list(APPEND nec_languages "C")
    endif()
  endif()
  if("CXX" IN_LIST ICON_NEC_WORKAROUND_LANGUAGES AND CMAKE_CXX_COMPILER_ID
                                                     STREQUAL "GNU")
    include(CheckCXXSymbolExists)
    check_cxx_symbol_exists(__NEC__ "" CXX_is_NEC)
    if(CXX_is_NEC)
      set(abi_CXX_file "${CMAKE_ROOT}/Modules/CMakeCXXCompilerABI.cpp")
      list(APPEND nec_languages "CXX")
    endif()
  endif()
  if("Fortran" IN_LIST ICON_NEC_WORKAROUND_LANGUAGES
     AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    file(
      WRITE "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.F90"
      "      program main
      implicit none
#ifdef __NEC__
      integer a
#else
      choke me
#endif
#ifndef __NEC__
      choke me
#else
      integer b
#endif
      a = 4
      b = 2
      end
")
    try_compile(
      Fortran_is_NEC "${PROJECT_BINARY_DIR}"
      "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.F90")
    if(Fortran_is_NEC)
      if(NOT CMAKE_Fortran_COMPILE_OPTIONS_PREPROCESS_ON)
        set(CMAKE_Fortran_COMPILE_OPTIONS_PREPROCESS_ON -fpp)
      endif()
      set(abi_Fortran_file "${CMAKE_ROOT}/Modules/CMakeFortranCompilerABI.F")
      list(APPEND nec_languages "Fortran")
    endif()
  endif()
  # Re-run the checks for internal compiler flags:
  if(nec_languages)
    message(STATUS "Applying the NEC workarounds")
    include(CMakeDetermineCompilerABI)
    # NEC compilers use 'nld' instead of 'ld' for linking:
    set(CMAKE_LINK_STARTFILE "nld")
    foreach(language IN LISTS nec_languages)
      message(STATUS "The real ${language} compiler ID is NEC")
      unset(CMAKE_${language}_ABI_COMPILED)
      cmake_determine_compiler_abi(${language} "${abi_${language}_file}")
    endforeach()
  endif()
endif()

# Set the path to the root source directory of the bundled library:
set(ICON_DEP_SOURCE_DIR
    ""
    CACHE PATH "Path to the source directory of the ICON dependency")
if(NOT ICON_DEP_SOURCE_DIR)
  message(
    FATAL_ERROR
      "The path to the root source directory of the bundled library ICON_DEP_SOURCE_DIR is missing"
  )
endif()

# Set the path to the root build directory of the bundled library:
set(dep_build_dir "_icon")

add_subdirectory("${ICON_DEP_SOURCE_DIR}" "${dep_build_dir}")

# Set the list of libraries (targets) that ICON depends on:
set(ICON_DEP_LIBRARIES
    ""
    CACHE STRING "List of libraries (targets) that ICON needs to be linked to")

if(ICON_DEP_LIBRARIES)
  # If the list of dependencies is not empty, we activate the linker flag
  # extraction mechanism. First, all languages that are enabled for the bundled
  # library must be enabled for the current project:
  get_property(dep_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
  enable_language(${dep_languages})

  # Make sure that Fortran is enabled:
  if(NOT "Fortran" IN_LIST dep_languages)
    enable_language(Fortran)
  endif()

  # Create a dummy executable (we have to set a source file to avoid complains
  # from CMake, so we give it this file):
  add_executable(icon-extra-libs EXCLUDE_FROM_ALL "${CMAKE_CURRENT_LIST_FILE}")

  # The executable depends on libraries that ICON needs to be linked to:
  target_link_libraries(icon-extra-libs PRIVATE "${ICON_DEP_LIBRARIES}")

  # To be able to extract all linker flags, we create a custom language:
  set_target_properties(icon-extra-libs PROPERTIES LINKER_LANGUAGE
                                                   IconExtraLibs)

  # Define a list of variables that the custom language should inherit from
  # Fortran language to make the <LINK_LIBRARIES> below complete:
  set(fortran_variables CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES
                        CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)

  # Append the RPATH flags if requested:
  set(ICON_ENABLE_RPATHS
      ON
      CACHE BOOL "Enable runtime library search paths")
  if(ICON_ENABLE_RPATHS)
    list(APPEND fortran_variables CMAKE_EXECUTABLE_RUNTIME_Fortran_FLAG)
  endif()

  # Copy the variables from Fortran to the custom language:
  foreach(fortran_variable IN LISTS fortran_variables)
    string(REPLACE Fortran IconExtraLibs icon_flag_language_variable
                   ${fortran_variable})
    set(${icon_flag_language_variable} ${${fortran_variable}})
  endforeach()

  # Set the custom language linker command so that the resulting link.txt file
  # contains only the required linker flags with properly quoted whitespaces:
  set(CMAKE_IconExtraLibs_LINK_EXECUTABLE "<LINK_LIBRARIES>")
endif()

# Allow for 'make test' even if the bundled library does not provide any tests
# or the tests are disabled:
enable_testing()
