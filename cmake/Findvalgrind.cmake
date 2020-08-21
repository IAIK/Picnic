# - Try to find valgrind
# Once done this will define
#  VALGRIND_FOUND - System has valgrind
#  VALGRIND_EXECUTABLE - valgrind executable

find_package(PkgConfig)
pkg_check_modules(PC_VALGRIND QUIET valgrind)
set(VALGRIND_VERSION_STRING ${PC_VALGRIND_VERSION})

find_path(VALGRIND_INCLUDE_DIR valgrind.h
          HINTS ${PC_VALGRIND_INCLUDEDIR} ${PC_VALGRIND_INCLUDE_DIRS}
          PATH_SUFFIXES valgrind)
find_program(VALGRIND_EXECUTABLE NAMES valgrind)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VALGRIND_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(valgrind
                                  FOUND_VAR VALGRIND_FOUND
                                  REQUIRED_VARS VALGRIND_INCLUDE_DIR VALGRIND_EXECUTABLE
                                  VERSION_VAR VALGRIND_VERSION_STRING)

mark_as_advanced(VALGRIND_INCLUDE_DIR
                 VALGRIND_VERSION_STRING
                 VALGRIND_EXECUTABLE)
