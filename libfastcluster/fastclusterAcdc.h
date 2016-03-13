#pragma once

// Same methods as in fastcluster.cpp
enum ClusterMethods 
{
	FC_METHOD_METR_SINGLE           = 0,
	FC_METHOD_METR_COMPLETE         = 1,
	FC_METHOD_METR_AVERAGE          = 2,
	FC_METHOD_METR_WEIGHTED         = 3,
	FC_METHOD_METR_WARD             = 4,
	FC_METHOD_METR_WARD_D           = FC_METHOD_METR_WARD,
	FC_METHOD_METR_CENTROID         = 5,
	FC_METHOD_METR_MEDIAN           = 6,
	FC_METHOD_METR_WARD_D2          = 7,
};

void fc_linkage(double * const D_, double * const Z_, long int N_, unsigned char method);