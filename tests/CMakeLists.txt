# Generate tests for all parameters
function(add_picnic_tests name target m1)
  if(WITH_ZKBPP)
    list(APPEND ALL_PARAMS
         picnic_L1_FS
         picnic_L1_UR
         picnic_L3_FS
         picnic_L3_UR
         picnic_L5_FS
         picnic_L5_UR
         picnic_L1_full
         picnic_L3_full
         picnic_L5_full)
  endif()
  if(WITH_KKW)
    list(APPEND ALL_PARAMS
         picnic3_L1_FS
         picnic3_L3_FS
         picnic3_L5_FS)
  endif()

  foreach(param IN LISTS ALL_PARAMS)
    add_test(NAME "${name}_${param}" COMMAND ${target} ${param})
  endforeach(param)
endfunction()

if(M4RI_FOUND)
  foreach(target IN ITEMS mzd)
    add_executable("${target}_test" "${target}_test.c")
    target_link_libraries("${target}_test" picnic_static M4RI::m4ri m)
    apply_base_options("${target}_test")
    apply_opt_options("${target}_test")

    add_test(NAME ${target} COMMAND "${target}_test")
  endforeach(target)
endif()

list(APPEND test_static_targets shift kdf_shake256 shake bitstream lowmc extended_picnic)
list(APPEND test_targets)
if(WITH_KKW)
  list(APPEND test_static_targets transpose)
endif()
if(NOT WIN32)
  list(APPEND test_targets picnic)
else()
  list(APPEND test_static_targets picnic)
endif()
foreach(target IN ITEMS ${test_static_targets})
  add_executable("${target}_test" "${target}_test.c")
  target_link_libraries("${target}_test" picnic_static)
  if(${target} STREQUAL lowmc)
    apply_picnic_options("${target}_test")
  else()
    apply_base_options("${target}_test")
    apply_opt_options("${target}_test")
  endif()

  if(${target} STREQUAL picnic OR ${target} STREQUAL extended_picnic)
    add_picnic_tests(${target} "${target}_test" ON)
  else()
    add_test(NAME ${target} COMMAND "${target}_test")
  endif()
endforeach(target)

foreach(target IN ITEMS ${test_targets})
  add_executable("${target}_test" "${target}_test.c")
  target_link_libraries("${target}_test" picnic)
  apply_base_options("${target}_test")

  if(${target} STREQUAL picnic OR ${target} STREQUAL extended_picnic)
    add_picnic_tests(${target} "${target}_test" ON)
  else()
    add_test(NAME ${target} COMMAND "${target}_test")
  endif()
endforeach(target)

if(NOT WITH_EXTRA_RANDOMNESS AND WITH_CONFIG_H)
  add_executable(kats_test kats_test.c)
  if (NOT WIN32)
    target_link_libraries(kats_test picnic)
  else()
    target_link_libraries(kats_test picnic_static)
  endif()
  apply_base_options(kats_test)
  target_compile_definitions(kats_test PRIVATE "-DKATDIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")

  add_picnic_tests(kats kats_test OFF)
endif()

if(NOT WIN32)
  foreach(target IN ITEMS ${api_targets})
    add_executable("api_${target}_test" api_test.c)
    apply_base_options("api_${target}_test")
    target_link_libraries("api_${target}_test" "${target}")
    add_test(NAME "api_${target}" COMMAND "api_${target}_test")
  endforeach(target)
endif()
