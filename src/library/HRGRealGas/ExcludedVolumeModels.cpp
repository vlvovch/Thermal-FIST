#include "HRGRealGas/ExcludedVolumeModels.h"
#include <cstdio>
#include <stdexcept>
#include <string>

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
		switch (n) {
			case 0:
				return f(eta);
			case 1:
				return d1f(eta);
			case 2:
				return d2f(eta);
			case 3:
				return d3f(eta);
			case 4:
				return d4f(eta);
			default:
				throw std::invalid_argument("ExcludedVolumeModelCS::df(n,eta): n = " + std::to_string(n) + " not supported!");
		}
	}

	double ExcludedVolumeModelTVM::df(int n, double eta) const {
		switch (n) {
			case 0:
				return f(eta);
			case 1:
				return d1f(eta);
			case 2:
				return d2f(eta);
			case 3:
				return d3f(eta);
			case 4:
				return d4f(eta);
			default:
				throw std::invalid_argument("ExcludedVolumeModelTVM::df(n,eta): n = " + std::to_string(n) + " not supported!");
		}
	}


} // namespace thermalfist