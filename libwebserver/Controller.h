#pragma once

#include <microhttpd.h>
#include <sstream>

class Controller
{

public:
    virtual bool validPath(const char * path, const char * method) = 0;

    virtual int handleRequest(  struct MHD_Connection * connection,
                                const char * url, const char * method, 
                                const char * upload_data, size_t * upload_data_size) = 0;

};

class DynamicController : public Controller
{

public:
    virtual bool validPath(const char * path, const char * method) = 0;

    virtual void createResponse(struct MHD_Connection * connection,
                                const char * url, const char * method, 
                                const char * upload_data, size_t * upload_data_size, 
                                std::stringstream & response) = 0;

    virtual int handleRequest(  struct MHD_Connection * connection,
                                const char * url, const char * method, 
                                const char * upload_data, size_t * upload_data_size);
};

class SimpleGetController : public DynamicController
{
private:
    std::string path;

public:
    SimpleGetController(const std::string path_);
    ~SimpleGetController();

    bool validPath(const char * path, const char * method);

    void createResponse(struct MHD_Connection * connection,
                                const char * url, const char * method, 
                                const char * upload_data, size_t * upload_data_size, 
                                std::stringstream & response);

    virtual void respond(std::stringstream & response) = 0;
};

