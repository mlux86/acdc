message(STATUS ${CMAKE_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/tsne)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(testmain 
	testmain.cpp
	SequenceVectorizer.cpp
)

target_link_libraries(testmain mainprog tsne ${Boost_LIBRARIES})

add_test(all testmain)