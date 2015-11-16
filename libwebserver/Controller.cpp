#include "Controller.h"

#include <cstring>
#include <map>
#include <string>
#include <iostream>

int DynamicController::handleRequest(  struct MHD_Connection * connection,
                            const char * url, const char * method, 
                            const char * upload_data, size_t * upload_data_size)
{
    std::stringstream response_string;
    createResponse(connection, url, method, upload_data, upload_data_size, response_string);

    //Send response.
    struct MHD_Response * response = MHD_create_response_from_buffer(
        strlen(response_string.str().c_str()),
        (void *)response_string.str().c_str(), 
        MHD_RESPMEM_MUST_COPY);

    int ret = MHD_queue_response(connection, MHD_HTTP_OK, response);
    
    MHD_destroy_response(response);
    
    return ret;
}

SimpleGetController::SimpleGetController(const std::string path_) : path(path_)
{
}

SimpleGetController::~SimpleGetController()
{
}

bool SimpleGetController::validPath(const char * path_, const char * method_)
{
	return strcmp(path_, this->path.c_str()) == 0 && strcmp("GET", method_) == 0;
}

void SimpleGetController::createResponse(	struct MHD_Connection * connection,
                                			const char * url, const char * method, 
                                			const char * upload_data, size_t * upload_data_size, 
                                			std::stringstream & response)
{
    std::map<std::string, std::string> params;
    MHD_get_connection_values(connection, MHD_GET_ARGUMENT_KIND, MHDCollectParams, &params);
	respond(response, params);
}

int SimpleGetController::MHDCollectParams(void * cls, enum MHD_ValueKind kind, const char * key, const char * value)
{
    std::map<std::string, std::string> * params = static_cast< std::map<std::string, std::string> * >(cls);
    std::string k(key);
    std::string v(value);
    (*params)[k] = v;
    return MHD_YES;
}