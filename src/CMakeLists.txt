add_library(mainprog
	Opts.cpp
	SequenceUtil.cpp
	MatrixUtil.cpp
	IOUtil.cpp
	MLUtil.cpp
	BarnesHutSNEAdapter.cpp
	SequenceVectorizer.cpp
	Clustering.cpp
	HierarchicalClustering.cpp
	TarjansAlgorithm.cpp
	Kmeans.cpp
	ClusterAnalysis.cpp
	Logger.cpp
	KrakenAdapter.cpp
	ResultIO.cpp
)

add_executable(acdc main.cpp)

target_link_libraries(acdc mainprog tsne jsoncpp dip pthread fastcluster m ${Boost_LIBRARIES} ${JSONCPP_LIBRARY})
