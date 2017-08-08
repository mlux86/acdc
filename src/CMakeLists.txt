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
	ClusteringUtil.cpp
	TarjansAlgorithm.cpp
	ContaminationDetection.cpp
	Decontamination.cpp
	Logger.cpp
	KrakenAdapter.cpp
	ResultIO.cpp
	RnammerAdapter.cpp
)

add_executable(acdc main.cpp)

target_link_libraries(acdc mainprog tsne jsoncpp dip pthread fastcluster m ${Boost_LIBRARIES} ${JSONCPP_LIBRARY})
