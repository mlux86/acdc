#pragma once

extern "C" {
	void diptest(const double x[], const int n, double *dip, int *lo_hi, int *ifault, int *gcm, int *lcm, int *mn, int *mj);
}
