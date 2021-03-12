add_library (hwrel_bsc SHARED bschwrel.cpp)

target_include_directories (hwrel_bsc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)

install(TARGETS hwrel_bsc
    ARCHIVE  DESTINATION lib
    LIBRARY  DESTINATION lib)
install(FILES bschwrel.h ../libhwrel.h DESTINATION include/libhwrel)

