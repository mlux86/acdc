#include "WebServer.h"
#include "Controller.h"

#include <stdexcept>

int WebServer::requestHandler(
    void * cls, struct MHD_Connection * connection,
    const char * url, const char * method, const char * version,
    const char * upload_data, size_t * upload_data_size, void ** ptr) 
{
    WebServer * server = (WebServer *)cls;

    Controller * controller = 0;
    for (auto c : server->controllers)
    {
        if (c->validPath(url, method))
        {
            controller = c;
            break;
        }
    }

    if (!controller)
    {
        struct MHD_Response * response = MHD_create_response_from_buffer(0, 0, MHD_RESPMEM_PERSISTENT);
        return MHD_queue_response (connection, MHD_HTTP_NOT_FOUND, response);
    }

    return controller->handleRequest(connection, url, method, upload_data, upload_data_size);
}

WebServer::WebServer(unsigned p)
{
    port = p;
    daemon = nullptr;
}

WebServer::~WebServer()
{
}

void WebServer::addController(Controller * controller)
{
    controllers.push_back(controller);
}

void WebServer::start()
{
    daemon = MHD_start_daemon(MHD_USE_SELECT_INTERNALLY | MHD_USE_THREAD_PER_CONNECTION, port, NULL, NULL, &requestHandler, this, MHD_OPTION_END);

    if(!daemon)
        throw std::runtime_error("Failed to start webserver!");

    stopped = false;
    std::unique_lock<std::mutex> lock(stoppedMutex);
    while (!stopped) 
    {
        stoppedCv.wait(lock);
    }    
}

void WebServer::stop()
{
    MHD_stop_daemon(daemon);
    stopped = true;
    stoppedCv.notify_one();
}
