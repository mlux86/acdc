add_library(mainprog
	Opts.cpp
	SequenceUtil.cpp
	MatrixUtil.cpp
	IOUtil.cpp
	MLUtil.cpp
	BarnesHutSNEAdapter.cpp
	SequenceVectorizer.cpp
	ClusteringAlgorithm.cpp
	HierarchicalClustering.cpp
	ClusterEstimator.cpp
	ConnectedComponentsEstimator.cpp
	ExternalValidityEstimator.cpp
	ClusterPostProcessing.cpp
	ClusteringUtil.cpp
	TarjansAlgorithm.cpp
	ContaminationDetection.cpp
	Decontamination.cpp
	Logger.cpp
	KrakenAdapter.cpp
	ResultIO.cpp
	RnammerAdapter.cpp
	TaxonomyAnnotation.cpp
)

add_executable(acdc main.cpp)

target_link_libraries(acdc mainprog yaml-cpp tsne dip pthread fastcluster m ${Boost_LIBRARIES})
