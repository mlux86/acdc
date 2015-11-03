add_library(mainprog
	Opts.cpp
	Util.cpp
	BarnesHutSNEAdapter.cpp
	SequenceVectorizer.cpp
	Clustering.cpp
	TarjansAlgorithm.cpp
	Kmeans.cpp
	ClusterAnalysis.cpp
	Logger.cpp
)

include_directories(${CMAKE_SOURCE_DIR}/libtsne)
include_directories(${CMAKE_SOURCE_DIR}/libdip)

add_executable(acdc main.cpp)

target_link_libraries(acdc mainprog tsne dip ${Boost_LIBRARIES} ${SEQAN_LIBRARIES})
