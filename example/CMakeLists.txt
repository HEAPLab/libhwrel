add_executable (libhwrel_bsc_test example.cpp)
target_link_libraries(libhwrel_bsc_test hwrel_bsc)

target_include_directories (hwrel_bsc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories (hwrel_bsc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../bsc)
