include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

# Picnic tests and extended tests
if(Boost_UNIT_TEST_FRAMEWORK_FOUND)
  add_executable(picnic_test
    picnic_test.cpp
    $<$<BOOL:${WITH_EXTENDED_TESTS}>:extended_picnic_test.cpp>)

  if(NOT WIN32 AND NOT CYGWIN)
    target_compile_definitions(picnic_test PRIVATE BOOST_TEST_DYN_LINK=1)
    target_link_libraries(picnic_test utils Boost::unit_test_framework)
  else()
    target_link_libraries(picnic_test utils_static Boost::unit_test_framework)
  endif()

  add_test(NAME picnic COMMAND picnic_test)

  foreach(target IN ITEMS ${api_targets})
    STRING(TOLOWER ${target} param)
    add_executable("${target}_test" picnic_param_test.cpp)
    target_compile_definitions("${target}_test" PRIVATE PARAM=${param} "HEADER=<${param}.h>")
    target_include_directories("${target}_test" PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../${target}")

    if(NOT WIN32 AND NOT CYGWIN)
      target_compile_definitions(picnic_test PRIVATE BOOST_TEST_DYN_LINK=1)
      target_link_libraries("${target}_test" utils Boost::unit_test_framework)
    else()
      target_link_libraries("${target}_test" utils_static Boost::unit_test_framework)
    endif()

    add_test(NAME "${target}" COMMAND "${target}_test")
  endforeach(target)
endif()

# Generate tests for all parameters
function(add_picnic_tests name target)
  foreach(param IN LISTS api_targets)
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

list(APPEND test_static_targets shift kdf_shake256 shake lowmc bitstream)
list(APPEND test_targets)

foreach(target IN ITEMS ${test_static_targets})
  add_executable("${target}_test" "${target}_test.c")
  add_test(NAME ${target} COMMAND "${target}_test")
  target_link_libraries("${target}_test" picnic_static)
  apply_base_options("${target}_test")
  apply_opt_options("${target}_test")
endforeach(target)

foreach(target IN ITEMS ${test_targets})
  add_executable("${target}_test" "${target}_test.c")
  add_test(NAME ${target} COMMAND "${target}_test")
  target_link_libraries("${target}_test" picnic)
  apply_base_options("${target}_test")
endforeach(target)

if(WITH_KATS_TESTS AND NOT WITH_EXTRA_RANDOMNESS)
  add_executable(kats_test kats_test.cpp)

  if(NOT WIN32 AND NOT CYGWIN)
    target_link_libraries(kats_test utils)
  else()
    target_link_libraries(kats_test utils_static)
  endif()

  apply_base_options(kats_test)
  target_compile_definitions(kats_test PRIVATE "-DKATDIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")

  add_picnic_tests(kats kats_test)
endif()

if(NOT WIN32 AND NOT CYGWIN)
  foreach(target IN ITEMS ${api_targets})
    add_executable("api_${target}_test" api_test.c)
    apply_base_options("api_${target}_test")
    target_link_libraries("api_${target}_test" "${target}")

    if(VALGRIND_FOUND AND WITH_VALGRIND_CT_TESTS)
      # if valgrind CT tests are enabled, execute tests with valgrind
      add_test(NAME "api_${target}" COMMAND "${VALGRIND_EXECUTABLE}" -q --error-exitcode=1 --track-origins=yes "$<TARGET_FILE:api_${target}_test>")
    else()
      add_test(NAME "api_${target}" COMMAND "api_${target}_test")
    endif()
  endforeach(target)
endif()
