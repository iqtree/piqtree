cmake_minimum_required(VERSION 3.21...4.2)

project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

find_package(
  Python 3.13 REQUIRED
  COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT})
find_package(nanobind CONFIG REQUIRED)

# The IQ-TREE library, built by ./build_tools/build_iqtree.sh
set(LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/piqtree/_libiqtree")

if(WIN32)
  add_library(iqtree SHARED IMPORTED)
  set_target_properties(
    iqtree PROPERTIES IMPORTED_IMPLIB "${LIBRARY_DIR}/iqtree.lib"
                      IMPORTED_LOCATION "${LIBRARY_DIR}/iqtree.dll")
  set(IQTREE_FILES "${LIBRARY_DIR}/iqtree.lib" "${LIBRARY_DIR}/iqtree.dll")
else()
  add_library(iqtree STATIC IMPORTED)
  set_target_properties(iqtree PROPERTIES IMPORTED_LOCATION
                                          "${LIBRARY_DIR}/libiqtree.a")
  set(IQTREE_FILES "${LIBRARY_DIR}/libiqtree.a")

  # Link the OpenMP runtime IQ-TREE was built against.
  if(APPLE)
    # Homebrew's OpenMP runtime, from libomp or bundled with llvm.
    foreach(package libomp llvm)
      execute_process(
        COMMAND brew --prefix ${package}
        OUTPUT_VARIABLE brew_prefix
        OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
      list(APPEND omp_hints "${brew_prefix}/lib")
    endforeach()
    find_library(OMP_LIBRARY omp HINTS ${omp_hints} REQUIRED)
    target_link_libraries(iqtree INTERFACE "${OMP_LIBRARY}")
  else()
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    target_link_libraries(iqtree INTERFACE OpenMP::OpenMP_CXX)
  endif()

  find_package(ZLIB REQUIRED)
  target_link_libraries(iqtree INTERFACE ZLIB::ZLIB)
endif()

foreach(file IN LISTS IQTREE_FILES)
  if(NOT EXISTS "${file}")
    message(FATAL_ERROR "${file} not found. "
                        "Build IQ-TREE first with ./build_tools/build_iqtree.sh")
  endif()
endforeach()

if(NB_FREE_THREADED)
  # One abi3t wheel supports free-threaded Python >= 3.15. This needs nanobind's
  # split mode, which loads nanobind from the nanobind-backend package.
  nanobind_add_module(_piqtree STABLE_ABI FREE_THREADED
                      BACKEND_MODULE nanobind_backend
                      src/piqtree/_libiqtree/_piqtree.cpp)
else()
  # STABLE_ABI builds an abi3 extension, so one wheel supports Python >= 3.13.
  nanobind_add_module(_piqtree STABLE_ABI NB_STATIC
                      src/piqtree/_libiqtree/_piqtree.cpp)
endif()
target_link_libraries(_piqtree PRIVATE iqtree)

# Only export the module's init function, rather than every symbol in the
# IQ-TREE library. This also lets the linker discard IQ-TREE code the module
# never uses.
if(APPLE)
  target_link_options(_piqtree PRIVATE
                      "LINKER:-exported_symbol,_PyInit__piqtree")
elseif(NOT WIN32)
  target_link_options(_piqtree PRIVATE "LINKER:--exclude-libs,ALL")
endif()

# Editable installs import piqtree from src/, so install the module (and on
# Windows, IQ-TREE's DLLs) there, where Python and editors can find it.
if(SKBUILD_STATE STREQUAL "editable")
  set(PACKAGE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/piqtree")
else()
  set(PACKAGE_DIR piqtree)
endif()

install(TARGETS _piqtree LIBRARY DESTINATION "${PACKAGE_DIR}")

# Wheels built by cibuildwheel turn this off, as delvewheel vendors the DLLs
# under mangled names instead.
option(PIQTREE_BUNDLE_DLLS "Install IQ-TREE's DLLs into piqtree/_libiqtree" ON)

if(WIN32 AND PIQTREE_BUNDLE_DLLS)
  # Bundle iqtree.dll and the (non-system) DLLs it depends on, searching PATH.
  file(TO_CMAKE_PATH "$ENV{PATH}" dll_search_path)
  install(
    IMPORTED_RUNTIME_ARTIFACTS iqtree RUNTIME_DEPENDENCY_SET iqtree_deps
    RUNTIME DESTINATION "${PACKAGE_DIR}/_libiqtree")
  install(
    RUNTIME_DEPENDENCY_SET iqtree_deps
    PRE_EXCLUDE_REGEXES "api-ms-.*" "ext-ms-.*"
    # Paths may mix separators (e.g. C:\Windows\system32/kernel32.dll).
    POST_EXCLUDE_REGEXES ".*[/\\\\][Ss][Yy][Ss][Tt][Ee][Mm]32[/\\\\].*"
    DIRECTORIES "${LIBRARY_DIR}" ${dll_search_path}
    RUNTIME DESTINATION "${PACKAGE_DIR}/_libiqtree")
endif()
