#pragma once

#include <condition_variable>
#include <mutex>
#include <vector>
#include <microhttpd.h>
#include "Controller.h"

class WebServer
{

private:
    unsigned port;
    struct MHD_Daemon * daemon;
    std::vector<Controller *> controllers;

    bool stopped = false;
    std::mutex stoppedMutex;
    std::condition_variable stoppedCv;

    static int requestHandler( void * cls, struct MHD_Connection * connection,
                                const char * url, const char * method, const char * version,
                                const char * upload_data, size_t * upload_data_size, void ** ptr);
public:
    WebServer(unsigned p);
    ~WebServer();

    void addController(Controller * controller);
    void start();
    void stop();

};
