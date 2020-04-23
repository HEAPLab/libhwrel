add_executable (libhwrel_bsc_test example.cpp)
target_link_libraries(libhwrel_bsc_test libhwrel_bsc)

target_include_directories (libhwrel_bsc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories (libhwrel_bsc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../bsc)
