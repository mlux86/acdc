add_library(webserver
	Controller.cpp
	WebServer.cpp
)

target_link_libraries(webserver ${MHD_LIBRARIES})