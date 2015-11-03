#ifndef __Logger__
#define __Logger__

#define ELOG Logger::getInstance().error()
#define ILOG Logger::getInstance().info()
#define VLOG Logger::getInstance().verbose()
#define DLOG Logger::getInstance().debug()

#include <iostream>
#include <fstream>

enum LogLevel
{
	Off,
	Error,
	Info,
	Verbose,
	Debug
};

class Logger
{

private:
	Logger();
	~Logger();
	Logger(Logger const &) = delete;
    void operator=(Logger const &) = delete;

    LogLevel level = Off;

    std::ofstream nullOut; 

public:
	static Logger & getInstance();

	void setLevel(LogLevel newLevel);

	std::ostream & log(LogLevel lvl);
	std::ostream & error();
	std::ostream & info();
	std::ostream & verbose();
	std::ostream & debug();

};

#endif 
