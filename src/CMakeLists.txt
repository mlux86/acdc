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
	VisualizationServer.cpp
	KrakenAdapter.cpp
)

add_executable(acdc main.cpp)

target_link_libraries(acdc mainprog webserver tsne dip ${Boost_LIBRARIES} ${SEQAN_LIBRARIES} ${JSONCPP_LIBRARY})
