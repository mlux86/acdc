add_library(dip
	dip.c
	DipStatistic.cpp
)

set_target_properties(dip PROPERTIES COMPILE_FLAGS "-O2")