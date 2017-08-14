message(STATUS ${CMAKE_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(acdc-test 
	testmain.cpp
	SequenceVectorizer.cpp
	DipStatistic.cpp
	TarjansAlgorithm.cpp
)

target_link_libraries(acdc-test mainprog yaml-cpp tsne dip pthread m ${Boost_LIBRARIES})

add_test(all acdc-test)