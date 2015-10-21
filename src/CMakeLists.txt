add_library(mainprog
	Util.cpp
	BarnesHutSNEBridge.cpp
	SequenceVectorizer.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/tsne)

add_executable(main main.cpp)

target_link_libraries(main mainprog tsne ${Boost_LIBRARIES})
