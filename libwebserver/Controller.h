#pragma once

#include <microhttpd.h>
#include <sstream>
#include <map>
#include <string>

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
    virtual ~SimpleGetController();

    virtual bool validPath(const char * path, const char * method);

    virtual void createResponse(struct MHD_Connection * connection,
                                const char * url, const char * method, 
                                const char * upload_data, size_t * upload_data_size, 
                                std::stringstream & response);

    virtual void respond(std::stringstream & response, const std::map<std::string, std::string> params) = 0;

    static int MHDCollectParams(void * cls, enum MHD_ValueKind kind, const char * key, const char * value);
};

