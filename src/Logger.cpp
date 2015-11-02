#include "Logger.h"

#include <thread>

Logger::Logger()
{
	nullOut.open("/dev/null");
}

Logger::~Logger()
{
	nullOut.close();
}

Logger & Logger::getInstance()
{
	static Logger instance;
	return instance;
}

void Logger::setLevel(LogLevel newLevel)
{
	this->level = newLevel;
}

std::ostream & Logger::log(LogLevel lvl)
{	
	std::ostream * out;

	switch (this->level)
	{
		case Off:
			out = &nullOut; break;
		case Error:
			out = (lvl == Error) ? &std::cerr : &nullOut; break;
		case Info:
			out = (lvl == Error || lvl == Info) ? &std::cout : &nullOut; break;
		case Verbose:
			out = (lvl == Error || lvl == Info || lvl == Verbose) ? &std::cout : &nullOut; break;
		case Debug:
			out = (lvl == Error || lvl == Info || lvl == Verbose || lvl == Debug) ? &std::cout : &nullOut; break;
		default:
			out = &nullOut; break;
	}

	std::hash<std::thread::id> hasher;
	char buff[100];
	snprintf(buff, sizeof(buff), "%lx", hasher(std::this_thread::get_id()));
	std::string buffAsStdStr = buff;

	(*out) << "[thread " << buffAsStdStr << "] ";
	return *out;
}

std::ostream & Logger::error()
{
	return log(Error);
}

std::ostream & Logger::info()
{
	return log(Info);
}

std::ostream & Logger::verbose()
{
	return log(Verbose);
}

std::ostream & Logger::debug()
{
	return log(Debug);
}