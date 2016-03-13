add_library(fastcluster
	fastclusterAcdc.cpp
)

set_target_properties(fastcluster PROPERTIES COMPILE_FLAGS "-O2")