message(STATUS ${CMAKE_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(acdc-test 
	testmain.cpp
	SequenceVectorizer.cpp
	DipStatistic.cpp
)

target_link_libraries(acdc-test mainprog tsne dip jsoncpp pthread m ${Boost_LIBRARIES} ${JSONCPP_LIBRARY})

add_test(all acdc-test)