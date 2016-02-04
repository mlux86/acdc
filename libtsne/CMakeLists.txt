add_library(tsne
	sptree.cpp
	tsne.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/src)

set_target_properties(tsne PROPERTIES COMPILE_FLAGS "-O2 -Wno-sign-compare -Wno-unused-variable -Wno-unused-result -Wno-unused-parameter -Wno-maybe-uninitialized")