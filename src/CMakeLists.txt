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
	ResultIO.cpp
	KrakenAdapter.cpp
)

add_executable(acdc main.cpp)

target_link_libraries(acdc mainprog tsne jsoncpp dip pthread m ${Boost_LIBRARIES} ${JSONCPP_LIBRARY})
