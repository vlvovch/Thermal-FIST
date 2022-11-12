#include "HRGRealGas/ExcludedVolumeModels.h"
#include <cstdio>

namespace thermalfist {

	double ExcludedVolumeModelBase::etasolBinarySearch(double etatil) const
	{
		double eps = 1.e-10;
		double left = 0., right = 1.;
		if (EtaMax() >= 0.0) 
			right = EtaMax() * (1. - eps);
		else {
			while (right / f(right) < etatil) right *= 2.;
		}
		double center = (left + right) / 2.;
		double valleft = left / f(left) - etatil;
		double valcenter = 0.;
		while ((right - left) / center > eps) {
			valcenter = center / f(center) - etatil;
			if (valleft * valcenter < 0.) {
				right = center;
			}
			else {
				left = center;
				valleft = valcenter;
			}
			center = (left + right) / 2.;
		}

		//printf("\n");
		//printf("%15lf %15lf %15E\n", etatil, center / etatil, center / etatil - f(center));

		return etatil * f(center);
	}

	double ExcludedVolumeModelCS::df(int n, double eta) const {
		if (n == 0)
			return f(eta);
		if (n == 1)
			return d1f(eta);
		if (n == 2)
			return d2f(eta);
		if (n == 3)
			return d3f(eta);
		if (n == 4)
			return d4f(eta);

		printf("**ERROR** ExcludedVolumeModelCS::df(n,eta): n = %lf not supported!", n);
		exit(-1);
		return 0.;
	}

	double ExcludedVolumeModelTVM::df(int n, double eta) const {
		if (n == 0)
			return f(eta);
		if (n == 1)
			return d1f(eta);
		if (n == 2)
			return d2f(eta);
		if (n == 3)
			return d3f(eta);
		if (n == 4)
			return d4f(eta);

		printf("**ERROR** ExcludedVolumeModelTVM::df(n,eta): n = %lf not supported!", n);
		exit(-1);
		return 0.;
	}

} // namespace thermalfist