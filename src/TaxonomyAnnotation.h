#pragma once

#include <string>
#include <yaml-cpp/yaml.h>

// forward declaration
struct ResultContainer;

enum Type
{
	Annotated, Estimated
};

enum ContaminationState
{
	NA, Clean, Contamination
};

struct Taxonomy
{
	Type type;
	std::string identifier;
	double confidence;
	ContaminationState state;
};

YAML::Emitter & operator << (YAML::Emitter & out, const Taxonomy & t);

class TaxonomyAnnotation
{
private:
	TaxonomyAnnotation();
	~TaxonomyAnnotation();

	static void updateContaminationConfidences(ResultContainer & res);
	static void updateContaminationState(ResultContainer & res);

public:
	static void initialize(ResultContainer & res);
	static void annotateFromFile(ResultContainer & res, const std::string & filename);
	static void annotateUnknown(ResultContainer & res);

	
};