cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# cmake setup
# honor visibility properties
cmake_policy(SET CMP0063 NEW)
# enforce IPO if enabled
cmake_policy(SET CMP0069 NEW)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# project info
set(PICNIC_MAJOR_VERSION 3)
set(PICNIC_MINOR_VERSION 0)
set(PICNIC_PATCH_VERSION 17)
set(PICNIC_VERSION ${PICNIC_MAJOR_VERSION}.${PICNIC_MINOR_VERSION}.${PICNIC_PATCH_VERSION})

project(picnic LANGUAGES C CXX VERSION ${PICNIC_VERSION})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build." FORCE)
endif()

# set required C standard version
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# set required C++ standard version
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
find_package(valgrind)
find_package(Boost 1.65 OPTIONAL_COMPONENTS unit_test_framework program_options)

if(APPLE)
  find_path(SECRUTY_INCLUDE_DIR Security/Security.h)
  find_library(SECURITY_LIBRARY Security)
  mark_as_advanced(SECURITY_INCLUDE_DIR SECURITY_LIBRARY)
  set(HAVE_SECURITY_FRAMEWORK TRUE)
endif()

# check programs
find_program(SED sed)

# check headers
check_include_files(sys/auxv.h HAVE_SYS_AUXV_H)
check_include_files(asm/hwcap.h HAVE_ASM_HWCAP_H)
check_include_files(sys/random.h HAVE_SYS_RANDOM_H)

# check availability of some functions
check_symbol_exists(aligned_alloc stdlib.h HAVE_ALIGNED_ALLOC)
check_symbol_exists(posix_memalign stdlib.h HAVE_POSIX_MEMALIGN)
check_symbol_exists(memalign malloc.h HAVE_MEMALIGN)
check_symbol_exists(getrandom sys/random.h HAVE_GETRANDOM)
check_symbol_exists(arc4random_buf stdlib.h HAVE_ARC4RANDOM_BUF)
check_symbol_exists(getline stdio.h HAVE_GETLINE)
check_symbol_exists(explicit_bzero string.h HAVE_EXPLICIT_BZERO)
check_symbol_exists(consttime_memequal string.h HAVE_CONSTTIME_MEMEQUAL)
check_symbol_exists(timingsafe_bcmp string.h HAVE_TIMINGSAFE_BCMP)

# check supported types
check_type_size(ssize_t SSIZE_T LANGUAGE C)

# check supported compiler flags
check_c_compiler_flag(-march=native CC_SUPPORTS_MARCH_NATIVE)
check_c_compiler_flag(-mtune=native CC_SUPPORTS_MTUNE_NATIVE)
check_c_compiler_flag(-O3 CC_SUPPORTS_O3)
check_c_compiler_flag(-fomit-frame-pointer CC_SUPPORTS_FOMIT_FRAME_POINTER)
check_c_compiler_flag_and_add(-Wall CC_SUPPORTS_WALL)
check_c_compiler_flag_and_add(-Wextra CC_SUPPORTS_WEXTRA)
check_c_compiler_flag_and_add(-Wshadow CC_SUPPORTS_WSHADOW)
check_c_compiler_flag_and_add(-Werror=vla CC_SUPPORTS_WERROR_VLA)

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

set(DEFAULT_SHA3_IMPL "opt64")
if(NOT DEFINED CMAKE_SIZEOF_VOID_P)
  message(WARNING "Unable to set default SHA3 implementation based on the platform since CMAKE_SIZEOF_VOID_P is not defined.")
else()
  if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(DEFAULT_SHA3_IMPL "plain32")
  endif()
endif()

set(WITH_ZKBPP ON CACHE BOOL "Enable Picnic instances based on ZKB++.")
set(WITH_KKW ON CACHE BOOL "Enable Picnic instances based on KKW.")
set(WITH_UNRUH ON CACHE BOOL "Enable Picnic instances based on ZKB++ and Unruh.")
if(NOT WITH_ZKBPP AND NOT WITH_KKW)
  message(FATAL_ERROR "At least one of WITH_ZKBPP and WITH_KKW is required.")
endif()
if(NOT WITH_ZKBPP AND WITH_UNRUH)
  message(STATUS "WITH_ZKBPP is required for WITH_UNRUH, so assuming WITH_UNRUH=OFF.")
endif()
set(WITH_LOWMC_128_128_20 ON CACHE BOOL "Enable LowMC instance 128-128-20.")
set(WITH_LOWMC_192_192_30 ON CACHE BOOL "Enable LowMC instance 192-192-30.")
set(WITH_LOWMC_256_256_38 ON CACHE BOOL "Enable LowMC instance 256-256-38.")
set(WITH_LOWMC_129_129_4 ON CACHE BOOL "Enable LowMC instance 129-129-4.")
set(WITH_LOWMC_192_192_4 ON CACHE BOOL "Enable LowMC instance 192-192-4.")
set(WITH_LOWMC_255_255_4 ON CACHE BOOL "Enable LowMC instance 255-255-4.")
if(WITH_KKW AND NOT WITH_LOWMC_129_129_4 AND NOT WITH_LOWMC_192_192_4 AND NOT WITH_LOWMC_255_255_4)
  message(FATAL_ERROR "At least one of WITH_LOWMC_129_129_4, WITH_LOWMC_192_192_4, and WITH_LOWMC_255_255_4 is required for WITH_KKW.")
endif()
if(NOT WITH_LOWMC_128_128_20 AND NOT WITH_LOWMC_192_192_30 AND NOT WITH_LOWMC_256_256_38 AND NOT WITH_LOWMC_129_129_4 AND NOT WITH_LOWMC_192_192_4 AND NOT WITH_LOWMC_255_255_4)
  message(FATAL_ERROR "At least one of the LowMC parameters is required.")
endif()

set(WITH_SIMD_OPT ON CACHE BOOL "Enable optimizations via SIMD.")
set(WITH_AVX2 ON CACHE BOOL "Use AVX2 and BMI2 if available.")
set(WITH_SSE2 ON CACHE BOOL "Use SSE2 if available.")
set(WITH_NEON ON CACHE BOOL "Use NEON if available.")
set(WITH_MARCH_NATIVE ${DEFAULT_WITH_MARCH_NATIVE} CACHE BOOL "Build with -march=native -mtune=native (if supported).")
set(WITH_LTO ON CACHE BOOL "Enable link-time optimization (if supported).")
set(WITH_SHA3_IMPL ${DEFAULT_SHA3_IMPL} CACHE STRING "Select SHA3 implementation.")
set_property(CACHE WITH_SHA3_IMPL PROPERTY STRINGS "opt64" "plain32" "avx2" "armv8a-neon" "s390-cpacf")
set(WITH_EXTRA_RANDOMNESS OFF CACHE BOOL "Feed extra random bytes to KDF (fault attack counter measure).")
set(WITH_CONFIG_H ON CACHE BOOL "Generate config.h. Disabling this option is discouraged. It is only available to test builds produced for SUPERCOP.")
if(MSVC)
  set(USE_STATIC_RUNTIME OFF CACHE BOOL "Use MSVC's static runtime for the static library.")
endif()
set(WITH_EXTENDED_TESTS ON CACHE BOOL "Run extended tests.")
set(WITH_KATS_TESTS ON CACHE BOOL "Run known-answer tests.")
set(WITH_VALGRIND_CT_TESTS OFF CACHE BOOL "Build crypto_api tests with constant-time checks. Enabling this option is discouraged. It is only available to test builds for SUPERCOP/TIMECOP.")

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
  # no function prototype in MSVC headers
  list(APPEND WARNINGS 4255)
  # Spectre mitigation
  list(APPEND WARNINGS 5045)
  # initialization of aggregate types with non-constant values is valid C99
  list(APPEND WARNINGS 4204)
  # nameless structs/unions are part of C11
  list(APPEND WARNINGS 4201)
  # conditional expression is constant
  list(APPEND WARNINGS 4127)
  # unary minus of unsigned types
  list(APPEND WARNINGS 4146)
  foreach(warning IN ITEMS ${WARNINGS})
    add_compile_options("/wd${warning}")
  endforeach(warning)
  # "deprecation warnings" of libc functions
  add_compile_options(-D_CRT_SECURE_NO_WARNINGS)
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

# check if LTO supported
if(WITH_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT LTO_SUPPORTED)
  if (LTO_SUPPORTED)
    # enable -fomit-frame-pointer
    if(CC_SUPPORTS_FOMIT_FRAME_POINTER AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
      add_compile_options(-fomit-frame-pointer)
    endif()
  else()
    message(STATUS "IPO is not supported.")
  endif()
else()
set(LTO_SUPPORTED FALSE)
endif()

message(STATUS "Picnic configuration:")
message(STATUS "ZKB++ enabled: ${WITH_ZKBPP}")
message(STATUS "KKW enabled: ${WITH_KKW}")
message(STATUS "SHA3 implementation: ${WITH_SHA3_IMPL}")
message(STATUS "LTO enabled: ${LTO_SUPPORTED}")
message(STATUS "SIMD enabled: ${WITH_SIMD_OPT}")
if(WITH_SIMD_OPT)
  message(STATUS "SSE2 supported: ${CC_SUPPORTS_SSE2} enabled: ${WITH_SSE2}")
  message(STATUS "AVX2 supported: ${CC_SUPPORTS_AVX2} enabled: ${WITH_AVX2}")
  message(STATUS "NEON supported: ${CC_SUPPORTS_NEON} enabled: ${WITH_NEON}")
endif()

# enable -O3
if(CC_SUPPORTS_O3 AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-O3)
endif()

if(WITH_CONFIG_H)
  configure_file(config.h.in config.h)
endif()

# generate sign.c for NIST / SUPERCOP
if(SED)
  list(APPEND ALL_PARAMETERS Picnic_L1_FS Picnic_L1_UR Picnic_L1_full Picnic3_L1)
  list(APPEND ALL_PARAMETERS Picnic_L3_FS Picnic_L3_UR Picnic_L3_full Picnic3_L3)
  list(APPEND ALL_PARAMETERS Picnic_L5_FS Picnic_L5_UR Picnic_L5_full Picnic3_L5)
  list(APPEND ALL_DIRECTORIES picnic_L1_FS picnic_L1_UR picnic_L1_full picnic3_L1)
  list(APPEND ALL_DIRECTORIES picnic_L3_FS picnic_L3_UR picnic_L3_full picnic3_L3)
  list(APPEND ALL_DIRECTORIES picnic_L5_FS picnic_L5_UR picnic_L5_full picnic3_L5)

  foreach(idx RANGE 11)
    list(GET ALL_PARAMETERS ${idx} param)
    list(GET ALL_DIRECTORIES ${idx} dir)
    add_custom_command(
      OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/crypto_sign.c"
      COMMAND ${SED}
      -n
      "s/PICNIC_INSTANCE/${param}/g;w ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/crypto_sign.c"
      "${CMAKE_CURRENT_SOURCE_DIR}/crypto_sign.c.in"
      DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/crypto_sign.c.in"
      VERBATIM)
  endforeach()
endif()

# SHAKE implementation
if(NOT WITH_SHA3_IMPL STREQUAL "s390-cpacf")
  list(APPEND SHA3_SOURCES
       sha3/KeccakHash.c
       sha3/KeccakSponge.c
       sha3/KeccakSpongetimes4.c
       sha3/KeccakHashtimes4.c)

  if(WITH_SHA3_IMPL STREQUAL "avx2")
    list(APPEND SHA3_SOURCES
         sha3/avx2/KeccakP-1600-AVX2.s
         sha3/avx2/KeccakP-1600-times4-SIMD256.c)
    set_property(SOURCE sha3/avx2/KeccakP-1600-AVX2.s PROPERTY LANGUAGE C)
  elseif(WITH_SHA3_IMPL STREQUAL "opt64")
    list(APPEND SHA3_SOURCES
         sha3/opt64/KeccakP-1600-opt64.c
         sha3/opt64/KeccakP-1600-times4-on1.c)
  elseif(WITH_SHA3_IMPL STREQUAL "plain32")
    list(APPEND SHA3_SOURCES
         sha3/plain32/KeccakP-1600-inplace32BI.c
         sha3/plain32/KeccakP-1600-times4-on1.c)
  elseif(WITH_SHA3_IMPL STREQUAL "armv8a-neon")
    list(APPEND SHA3_SOURCES
         sha3/armv8a-neon/KeccakP-1600-armv8a-neon.s
         sha3/opt64/KeccakP-1600-times4-on1.c)
    set_property(SOURCE sha3/armv8a-neon/KeccakP-1600-armv8a-neon.s PROPERTY LANGUAGE C)
  else()
    message(FATAL_ERROR "Unknown SHA3 implementation")
  endif()

  # 32 bit ARM is in general unhappy with unaligned access, hence disable it in this case
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "^arm" AND (WITH_SHA3_IMPL STREQUAL "opt64" OR WITH_SHA3_IMPL STREQUAL "plain32"))
    set(DISABLE_MISALIGNED_ACCESS TRUE)
  endif()
endif()

# Picnic implementation
list(APPEND PICNIC_SOURCES
     compat.c
     cpu.c
     io.c
     lowmc.c
     $<$<BOOL:${WITH_LOWMC_129_129_4}>:lowmc_129_129_4.c>
     $<$<BOOL:${WITH_LOWMC_192_192_4}>:lowmc_192_192_4.c>
     $<$<BOOL:${WITH_LOWMC_255_255_4}>:lowmc_255_255_4.c>
     mzd_additional.c
     picnic.c
     picnic_instances.c
     randomness.c)
if(WITH_ZKBPP)
  list(APPEND PICNIC_SOURCES
       $<$<BOOL:${WITH_LOWMC_128_128_20}>:lowmc_128_128_20.c>
       $<$<BOOL:${WITH_LOWMC_192_192_30}>:lowmc_192_192_30.c>
       $<$<BOOL:${WITH_LOWMC_256_256_38}>:lowmc_256_256_38.c>
       mpc_lowmc.c
       picnic_impl.c)
endif()
if(WITH_KKW)
  list(APPEND PICNIC_SOURCES
       picnic3_impl.c
       picnic3_simulate.c
       picnic3_tree.c
       picnic3_types.c)
endif()
list(APPEND PICNIC_HEADERS picnic.h)

# shared library
add_library(picnic SHARED ${PICNIC_SOURCES} ${SHA3_SOURCES})
set_target_properties(picnic PROPERTIES
                      VERSION ${PICNIC_VERSION}
                      SOVERSION ${PICNIC_MAJOR_VERSION}
                      INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED})
if(MSVC)
  set_target_properties(picnic PROPERTIES OUTPUT_NAME libpicnic)
endif()

# static library
add_library(picnic_static STATIC ${PICNIC_SOURCES} ${SHA3_SOURCES})
set_target_properties(picnic_static PROPERTIES
                      OUTPUT_NAME picnic
                      ARCHIVE_OUTPUT_DIRECTORY static
                      INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED})
target_compile_definitions(picnic_static PUBLIC PICNIC_STATIC)
# build static library with valgrind if enabled for tests
if(VALGRIND_FOUND AND WITH_VALGRIND_CT_TESTS)
  target_compile_definitions(picnic_static PUBLIC WITH_VALGRIND)
endif()

function(apply_base_options lib)
  target_compile_definitions(${lib} PRIVATE
                             $<$<BOOL:${WITH_ZKBPP}>:WITH_ZKBPP>
                             $<$<AND:$<BOOL:${WITH_UNRUH}>,$<BOOL:${WITH_ZKBPP}>>:WITH_UNRUH>
                             $<$<BOOL:${WITH_KKW}>:WITH_KKW>
                             $<$<BOOL:${WITH_CONFIG_H}>:HAVE_CONFIG_H>
                             $<$<BOOL:${DISABLE_MISALIGNED_ACCESS}>:NO_MISALIGNED_ACCESSES>)
  target_include_directories(${lib} PRIVATE
                             "${CMAKE_BINARY_DIR}"
                             "${CMAKE_SOURCE_DIR}"
                             "${CMAKE_SOURCE_DIR}/sha3"
                             "${CMAKE_SOURCE_DIR}/sha3/${WITH_SHA3_IMPL}")
  if(WITH_ZKBPP)
    target_compile_definitions(${lib} PRIVATE
                               $<$<BOOL:${WITH_LOWMC_128_128_20}>:WITH_LOWMC_128_128_20>
                               $<$<BOOL:${WITH_LOWMC_192_192_30}>:WITH_LOWMC_192_192_30>
                               $<$<BOOL:${WITH_LOWMC_256_256_38}>:WITH_LOWMC_256_256_38>)
  endif()
  target_compile_definitions(${lib} PRIVATE
                             $<$<BOOL:${WITH_LOWMC_129_129_4}>:WITH_LOWMC_129_129_4>
                             $<$<BOOL:${WITH_LOWMC_192_192_4}>:WITH_LOWMC_192_192_4>
                             $<$<BOOL:${WITH_LOWMC_255_255_4}>:WITH_LOWMC_255_255_4>)
  if(WITH_SHA3_IMPL STREQUAL "s390-cpacf")
    target_compile_definitions(${lib} PRIVATE WITH_SHAKE_S390_CPACF)
  else()
    target_compile_definitions(${lib} PRIVATE WITH_KECCAK_X4)
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
  target_compile_definitions(${lib} PRIVATE $<$<BOOL:${WITH_EXTRA_RANDOMNESS}>:WITH_EXTRA_RANDOMNESS>)

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

# pkg-config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/picnic.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/picnic.pc @ONLY)

add_library(utils STATIC tools/utils.cpp)
target_link_libraries(utils PUBLIC picnic)
target_include_directories(utils PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/tools")

add_library(utils_static STATIC tools/utils.cpp)
target_link_libraries(utils_static PUBLIC picnic_static)
target_include_directories(utils_static PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/tools")

if(Boost_PROGRAM_OPTIONS_FOUND)
  # bench exectuable
  add_executable(bench tools/bench.cpp)
  set_target_properties(bench PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}})
  target_link_libraries(bench utils Boost::program_options)
  apply_base_options(bench)

  # bench lowmc exectuable
  add_executable(bench_lowmc tools/bench_lowmc.cpp)
  set_target_properties(bench_lowmc PROPERTIES INTERPROCEDURAL_OPTIMIZATION ${LTO_SUPPORTED}})
  target_link_libraries(bench_lowmc utils_static Boost::program_options)
  apply_base_options(bench_lowmc)
  apply_opt_options(bench_lowmc)
endif()

# example executable
add_executable(example tools/example.c)
target_link_libraries(example picnic)
apply_base_options(example)

if(WITH_ZKBPP)
  if(WITH_LOWMC_128_128_20)
    list(APPEND api_targets picnic_L1_FS)
    if(WITH_UNRUH)
      list(APPEND api_targets picnic_L1_UR)
    endif()
  endif()
  if(WITH_LOWMC_192_192_30)
    list(APPEND api_targets picnic_L3_FS)
    if(WITH_UNRUH)
      list(APPEND api_targets picnic_L3_UR)
    endif()
  endif()
  if(WITH_LOWMC_256_256_38)
    list(APPEND api_targets picnic_L5_FS)
    if(WITH_UNRUH)
      list(APPEND api_targets picnic_L5_UR)
    endif()
  endif()
  if(WITH_LOWMC_129_129_4)
    list(APPEND api_targets picnic_L1_full)
  endif()
  if(WITH_LOWMC_192_192_4)
    list(APPEND api_targets picnic_L3_full)
  endif()
  if(WITH_LOWMC_255_255_4)
    list(APPEND api_targets picnic_L5_full)
  endif()
endif()
if(WITH_KKW)
  if(WITH_LOWMC_129_129_4)
    list(APPEND api_targets picnic3_L1)
  endif()
  if(WITH_LOWMC_192_192_4)
    list(APPEND api_targets picnic3_L3)
  endif()
  if(WITH_LOWMC_255_255_4)
    list(APPEND api_targets picnic3_L5)
  endif()
endif()

foreach(target IN ITEMS ${api_targets})
  add_library("${target}" STATIC "${target}/crypto_sign.c")
  target_include_directories("${target}" PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories("${target}" PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/${target}")
  if (VALGRIND_FOUND AND WITH_VALGRIND_CT_TESTS)
    target_link_libraries("${target}" picnic_static)
  else()
    target_link_libraries("${target}" picnic)
  endif()
endforeach(target)

# tests
enable_testing()
add_subdirectory(tests)

# install
install(TARGETS picnic picnic_static
        EXPORT picnic-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${PICNIC_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/picnic.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
install(EXPORT picnic-targets
        NAMESPACE picnic::
        FILE picnic-config.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/picnic)
