cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

# cmake setup
# honor visibility properties
cmake_policy(SET CMP0063 NEW)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER "3.8")
  # enforce IPO if enabled
  cmake_policy(SET CMP0069 NEW)
endif()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# project info
project(picnic C)

set(PICNIC_MAJOR_VERSION 2)
set(PICNIC_MINOR_VERSION 1)
set(PICNIC_PATCH_VERSION 1)
set(PICNIC_VERSION ${PICNIC_MAJOR_VERSION}.${PICNIC_MINOR_VERSION}.${PICNIC_PATCH_VERSION})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
endif()

# set required C standard version
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# cmake checks and tools
include(GNUInstallDirs)
include(CheckSIMD)
include(CheckFunctionExists)
include(CheckCCompilerFlag)
include(CheckIncludeFiles)
include(CheckSymbolExists)
include(CheckTypeSize)

function(check_c_compiler_flag_and_add flag result)
  check_c_compiler_flag("${flag}" ${result})
  if(${result})
    add_compile_options("${flag}")
  endif()
endfunction()

# check libraries
find_package(m4ri 20140914)

if(APPLE)
  find_path(SECRUTY_INCLUDE_DIR Security/Security.h)
  find_library(SECURITY_LIBRARY Security)
  mark_as_advanced(SECURITY_INCLUDE_DIR SECURITY_LIBRARY)
  set(HAVE_SECURITY_FRAMEWORK TRUE)
endif()

# check headers
check_include_files(sys/auxv.h HAVE_SYS_AUXV_H)
check_include_files(asm/hwcap.h HAVE_ASM_HWCAP_H)
check_include_files(sys/random.h HAVE_SYS_RANDOM_H)

# check availability of some functions
check_symbol_exists(aligned_alloc stdlib.h HAVE_ALIGNED_ALLOC)
check_symbol_exists(posix_memalign stdlib.h HAVE_POSIX_MEMALIGN)
check_symbol_exists(memalign malloc.h HAVE_MEMALIGN)
check_symbol_exists(getrandom sys/random.h HAVE_GETRANDOM)
check_symbol_exists(getline stdio.h HAVE_GETLINE)

# check supported types
check_type_size(ssize_t SSIZE_T LANGUAGE C)

# check supported compiler flags
check_c_compiler_flag(-march=native CC_SUPPORTS_MARCH_NATIVE)
check_c_compiler_flag(-mtune=native CC_SUPPORTS_MTUNE_NATIVE)
check_c_compiler_flag(-O3 CC_SUPPORTS_03)
check_c_compiler_flag(-fomit-frame-pointer CC_SUPPORTS_FOMIT_FRAME_POINTER)
check_c_compiler_flag_and_add(-Wall CC_SUPPORTS_WALL)
check_c_compiler_flag_and_add(-Wextra CC_SUPPORTS_WEXTRA)
check_c_compiler_flag_and_add(-Wshadow CC_SUPPORTS_WSHADOW)
check_c_compiler_flag_and_add(-Werror=implicit-function-declaration CC_SUPPORTS_WERROR_IMPLICIT_FUNCTION_DECLARATION)

# check SIMD instructions set
check_simd(SSE2 CC_SUPPORTS_SSE2)
check_simd(AVX2 CC_SUPPORTS_AVX2)
check_simd(BMI2 CC_SUPPORTS_BMI2)
check_simd(NEON CC_SUPPORTS_NEON)

# user-settable options
if(APPLE AND ${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  # workaround for broken -march=native support on some versions of GCC on OS X.
  set(DEFAULT_WITH_MARCH_NATIVE OFF)
else()
  set(DEFAULT_WITH_MARCH_NATIVE ON)
endif()

set(WITH_SIMD_OPT ON CACHE BOOL "Enable optimizations via SIMD.")
set(WITH_AVX2 ON CACHE BOOL "Use AVX2 and BMI2 if available.")
set(WITH_SSE2 ON CACHE BOOL "Use SSE2 if available.")
set(WITH_NEON ON CACHE BOOL "Use NEON if available.")
set(WITH_MARCH_NATIVE ${DEFAULT_WITH_MARCH_NATIVE} CACHE BOOL "Build with -march=native -mtune=native (if supported).")
set(WITH_LTO ON CACHE BOOL "Enable link-time optimization (if supported).")
set(SHA3_IMPL "opt64" CACHE STRING "Select SHA3 implementation.")
set(WITH_LOWMC_OPT "OLLE" CACHE STRING "Enable LowMC optimizations.")
set_property(CACHE WITH_LOWMC_OPT PROPERTY STRINGS "OFF" "ORKC" "OLLE")
set(WITH_LOWMC_M1 OFF CACHE BOOL "Enable LowMC instances with 1 Sbox.")
set(WITH_EXTRA_RANDOMNESS OFF CACHE BOOL "Feed extra random bytes to KDF (fault attack counter measure).")
set(WITH_CONFIG_H ON CACHE BOOL "Generate config.h.")
if(MSVC)
  set(USE_STATIC_RUNTIME OFF CACHE BOOL "Use MSVC's static runtime for the static library.")
endif()

# do not build with -rdynamic
string(REGEX REPLACE "-rdynamic" "" CMAKE_EXE_LINKER_FLAGS
  "${CMAKE_EXE_LINKER_FLAGS}")
string(REGEX REPLACE "-rdynamic" "" CMAKE_SHARED_LIBRARY_LINK_C_FLAGS
  "${CMAKE_SHARED_LIBRARY_LINK_C_FLAGS}")

# disable some warnings on MSVC
if(MSVC)
  # function inlining
  list(APPEND WARNINGS 4710 4711)
  # padding of structs
  list(APPEND WARNINGS 4820 4324)
  # .. in include paths
  list(APPEND WARNINGS 4464)
  # undefined macros evaluating to 0 in MSVC headers
  list(APPEND WARNINGS 4668)
  # Spectre mitigation
  list(APPEND WARNINGS 5045)
  # initialization of aggregate types with non-constant values is valid C99
  list(APPEND WARNINGS 4204)
  foreach(warning IN ITEMS ${WARNINGS})
    add_compile_options("/wd${warning}")
  endforeach(warning)
endif()

# enable -march=native -mtune=native if supported
if(WITH_MARCH_NATIVE)
  if (CC_SUPPORTS_MARCH_NATIVE)
    add_compile_options("-march=native")
  endif()
  if (CC_SUPPORTS_MTUNE_NATIVE)
    add_compile_options("-mtune=native")
  endif()
endif()

# enable LTO if supported
if(WITH_LTO)
  if ("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_LESS "3.9")
    check_c_compiler_flag(-flto CC_SUPPORTS_FLTO)
    if(CC_SUPPORTS_FLTO)
      if(CMAKE_COMPILER_IS_GNUCC)
        set(CMAKE_AR "gcc-ar")
        set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
        set(CMAKE_C_ARCHIVE_FINISH true)
      endif()

      add_compile_options(-flto)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto")
    else()
      message(STATUS "IPO is not supported.")
    endif()
  else()
    include(CheckIPOSupported)
    check_ipo_supported(RESULT LTO_SUPPORTED)
    if (LTO_SUPPORTED)
      set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
      message(STATUS "IPO is not supported.")
    endif()
  endif()
endif()

# enable -O3
if(CC_SUPPORTS_03 AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-O3)
endif()
# enable -fomit-frame-pointer
if(CC_SUPPORTS_FOMIT_FRAME_POINTER AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-fomit-frame-pointer)
endif()

if(WITH_CONFIG_H)
  configure_file(config.h.in config.h)
endif()

include_directories(${CMAKE_CURRENT_BINARY_DIR} sha3 "sha3/${SHA3_IMPL}")

# SHAKE implementation
list(APPEND SHA3_SOURCES
     sha3/KeccakHash.c
     sha3/KeccakSpongeWidth1600.c)

if(SHA3_IMPL STREQUAL "avx2")
  list(APPEND SHA3_SOURCES sha3/avx2/KeccakP-1600-AVX2.s)
  set_property(SOURCE sha3/avx2/KeccakP-1600-AVX2.s PROPERTY LANGUAGE C)
  # the times4 variant of Keccack using avx2
  list(APPEND SHA3_SOURCES
       sha3/avx2/KeccakP-1600-times4-SIMD256.c
       sha3/KeccakSpongeWidth1600times4.c
       sha3/KeccakHashtimes4.c)
elseif(SHA3_IMPL STREQUAL "opt64")
  list(APPEND SHA3_SOURCES sha3/opt64/KeccakP-1600-opt64.c)
elseif(SHA3_IMPL STREQUAL "armv8a-neon")
  list(APPEND SHA3_SOURCES sha3/armv8a-neon/KeccakP-1600-armv8a-neon.s)
  set_property(SOURCE sha3/armv8a-neon/KeccakP-1600-armv8a-neon.s PROPERTY LANGUAGE C)
else()
  message(FATAL_ERROR "Unknown SHA3 implementation")
endif()

# Picnic implementation
list(APPEND PICNIC_SOURCES
     aligned_alloc.c
     bitstream.c
     cpu.c
     io.c
     lowmc.c
     lowmc_128_128_20.c
     lowmc_192_192_30.c
     lowmc_256_256_38.c
     mpc_lowmc.c
     mzd_additional.c
     picnic.c
     picnic_impl.c
     picnic2_impl.c
     picnic2_simulate.c
     picnic2_simulate_mul.c
     picnic2_tree.c
     picnic2_types.c
     randomness.c)
if(WITH_LOWMC_M1)
  list(APPEND PICNIC_SOURCES
       lowmc_128_128_182.c
       lowmc_192_192_284.c
       lowmc_256_256_363.c)
endif()
list(APPEND PICNIC_HEADERS picnic.h)

# shared library
add_library(picnic SHARED ${PICNIC_SOURCES} ${SHA3_SOURCES})
set_target_properties(picnic PROPERTIES VERSION ${PICNIC_VERSION} SOVERSION ${PICNIC_MAJOR_VERSION})
if(MSVC)
  set_target_properties(picnic PROPERTIES OUTPUT_NAME libpicnic)
endif()

# static library
add_library(picnic_static STATIC ${PICNIC_SOURCES} ${SHA3_SOURCES})
if(MSVC)
  set_target_properties(picnic_static PROPERTIES OUTPUT_NAME libpicnic_static)
endif()
target_compile_definitions(picnic_static PUBLIC PICNIC_STATIC)

function(apply_base_options lib)
  if(WITH_CONFIG_H)
    target_compile_definitions(${lib} PRIVATE HAVE_CONFIG_H)
  endif()
  if(WITH_LOWMC_OPT STREQUAL "ORKC")
    target_compile_definitions(${lib} PRIVATE REDUCED_ROUND_KEY_COMPUTATION)
  elseif(WITH_LOWMC_OPT STREQUAL "OLLE")
    target_compile_definitions(${lib} PRIVATE REDUCED_ROUND_KEY_COMPUTATION OPTIMIZED_LINEAR_LAYER_EVALUATION)
  endif()
  if(SHA3_IMPL STREQUAL "avx2")
    target_compile_definitions(${lib} PRIVATE WITH_KECCAK_X4)
  endif()
  if(WITH_LOWMC_M1)
    target_compile_definitions(${lib} PRIVATE WITH_LOWMC_M1)
  endif()
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${lib} PRIVATE NDEBUG)
  endif()
endfunction()

function(apply_opt_options lib)
  if(WITH_SIMD_OPT)
    target_compile_definitions(${lib} PRIVATE WITH_OPT)
    if(CC_SUPPORTS_SSE2 AND WITH_SSE2)
      target_compile_definitions(${lib} PRIVATE WITH_SSE2)
      if(CC_SUPPORTS_AVX2 AND CC_SUPPORTS_BMI2 AND WITH_AVX2)
        target_compile_definitions(${lib} PRIVATE WITH_AVX2)
      endif()
    endif()
    if(CC_SUPPORTS_NEON AND WITH_NEON)
      target_compile_definitions(${lib} PRIVATE WITH_NEON)
    endif()
  endif()
endfunction()

function(apply_picnic_options lib)
  apply_base_options(${lib})
  apply_opt_options(${lib})

  set_target_properties(${lib} PROPERTIES C_VISIBILITY_PRESET hidden)
  target_compile_definitions(${lib} PRIVATE
                             WITH_LOWMC_128_128_20
                             WITH_LOWMC_192_192_30
                             WITH_LOWMC_256_256_38)
  if(WITH_LOWMC_M1)
    target_compile_definitions(${lib} PRIVATE
                               WITH_LOWMC_128_128_182
                               WITH_LOWMC_192_192_284
                               WITH_LOWMC_256_256_363)
  endif()
  if(WITH_EXTRA_RANDOMNESS)
    target_compile_definitions(${lib} PRIVATE WITH_EXTRA_RANDOMNESS)
  endif()

  if(WIN32)
    # require new enough Windows for bcrypt to be available
    target_compile_definitions(${lib} PRIVATE "_WIN32_WINNT=0x0601")
    target_link_libraries(${lib} bcrypt)
  endif()

  if(APPLE)
    target_include_directories(${lib} PRIVATE ${SECURTY_INCLUDE_DIR})
    target_link_libraries(${lib} PRIVATE ${SECURITY_LIBRARY})
  endif()
endfunction(apply_picnic_options)

apply_picnic_options(picnic)
apply_picnic_options(picnic_static)

if (WIN32)
  target_compile_definitions(picnic PRIVATE "PICNIC_EXPORT=__declspec(dllexport)")
else()
  target_compile_definitions(picnic PRIVATE "PICNIC_EXPORT=__attribute__((visibility(\"default\")))")
endif()
target_compile_definitions(picnic_static PRIVATE PICNIC_EXPORT=)
if(MSVC AND USE_STATIC_RUNTIME)
  target_compile_options(picnic_static PUBLIC "/MT$<$<CONFIG:Debug>:d>")
endif()
target_include_directories(picnic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(picnic_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# pkg-config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/picnic.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/picnic.pc @ONLY)

# visualize_signature executable
add_executable(visualize_signature tools/visualize_signature.c)
target_link_libraries(visualize_signature picnic_static)
apply_base_options(visualize_signature)

# bench utils library
add_library(bench_utils STATIC tools/bench_utils.c tools/bench_timing.c)
apply_base_options(bench_utils)

# bench exectuable
add_executable(bench tools/bench.c)
target_link_libraries(bench bench_utils picnic)
apply_base_options(bench)

# bench lowmc exectuable
add_executable(bench_lowmc tools/bench_lowmc.c)
target_link_libraries(bench_lowmc bench_utils picnic_static)
apply_base_options(bench_lowmc)
apply_opt_options(bench_lowmc)

# example executable
add_executable(example tools/example.c)
target_link_libraries(example picnic)
apply_base_options(example)

foreach(target IN ITEMS picnic_L1_FS picnic_L1_UR picnic_L3_FS picnic_L3_UR picnic_L5_FS picnic_L5_UR
                        picnic2_L1_FS picnic2_L3_FS picnic2_L5_FS)
  add_library("${target}" STATIC "${target}/sign.c")
  target_include_directories("${target}" PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories("${target}" PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/${target}")
  target_link_libraries("${target}" picnic)
endforeach(target)

# tests
enable_testing()
add_subdirectory(tests)

# install
install(TARGETS picnic DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS picnic_static DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${PICNIC_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picnic.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
