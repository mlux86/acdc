#pragma once

#include <string>

// forward declaration
struct ResultContainer;

enum ContaminationState
{
	NA, Clean, Contaminated
};

struct Taxonomy
{
	std::string identifier;
	std::string confidence;
	ContaminationState state;
};

class TaxonomyAnnotation
{
private:
	TaxonomyAnnotation();
	~TaxonomyAnnotation();

	static void updateContaminationState(ResultContainer & res);

public:
	static void annotateFromFile(ResultContainer & res, const std::string & filename);
	static void annotateContig(ResultContainer & res, const std::string & contigId, const std::string & species);
	static void annotateUnknown(ResultContainer & res);

	
};