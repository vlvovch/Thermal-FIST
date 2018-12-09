#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdio>

#include "HRGBase.h"
#include "HRGEV.h"
#include "HRGFit.h"
#include "HRGVDW/ThermalModelVDWFull.h"
#include "HRGEV/ExcludedVolumeHelper.h"

#include "ThermalFISTConfig.h"

using namespace std;

using namespace thermalfist;

// A macro for generating input files with excluded-volume/van der Waals parameters
// Modes:
// 0 - Constant
// 1 - Meson/Baryon(CI)
// 2 - Meson/Baryon(CII)
// 3 - Meson/Baryon(CIII)
// 4 - Bag model
// 5 - 2b
// 6 - 4b
// 7 - s-inv

void Diagonal(std::string filename, ThermalParticleSystem *TPS, int mode = 0, double r1 = 0., double r2 = 0., double r3 = 0., double r4 = 0.) {
	ThermalModelEVDiagonal model(TPS);
	vector<double> radii(TPS->Particles().size(), 0.);
	for (int i = 0; i < TPS->Particles().size(); ++i) {
		ThermalParticle &tpart = TPS->Particle(i);
		if (mode == 0) {
			radii[i] = r1;
		}
		if (mode >= 1 && mode <= 3) {
			if (tpart.BaryonCharge() == 0)
				radii[i] = r2;
			else
				radii[i] = r1 * pow(abs(tpart.BaryonCharge()), 1. / 3.);
		}
		if (mode == 4) {
			radii[i] = r1 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
		}
		if (mode == 5) {
			if (tpart.BaryonCharge() == 0)
				radii[i] = r2 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
			else
				radii[i] = r1 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
		}
		if (mode == 6) {
			if (tpart.BaryonCharge() == 0) {
				if (tpart.Strangeness() == 0)
					radii[i] = r2 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
				else 
					radii[i] = r4 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
			}
			else {
				if (tpart.Strangeness() == 0)
					radii[i] = r1 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
				else
					radii[i] = r3 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
			}
		}
		if (mode == 7) {
			if (tpart.Strangeness() == 0)
				radii[i] = r1 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
			else
				radii[i] = r2 * pow(tpart.Mass() / TPS->ParticleByPDG(3122).Mass() , -1. / 3.);
		}
	}

	model.FillVirial(radii);
	model.WriteInteractionParameters(filename);
}


void Crossterms(std::string filename, ThermalParticleSystem *TPS, int mode = 0, double r1 = 0., double r2 = 0., double r3 = 0., double r4 = 0.) {
	ThermalModelEVCrossterms model(TPS);
	vector<double> radii(TPS->Particles().size(), 0.);
	for (int i = 0; i < TPS->Particles().size(); ++i) {
		ThermalParticle &tpart = TPS->Particle(i);
		if (mode == 0) {
			radii[i] = r1;
		}
		if (mode >= 1 && mode <= 3) {
			//if (tpart.BaryonCharge() == 0)
			//	radii[i] = r2;
			//else
			//	radii[i] = r1 * pow(abs(tpart.BaryonCharge()), 1. / 3.);
			int N = TPS->Particles().size();
			for (int j = 0; j < N; ++j) {
				ThermalParticle &tpart2 = TPS->Particle(j);
				double tr1 = r1, tr2 = r2;
				//bool mes1 = (tpart.BaryonCharge() == 0), mes2 = (tpart2.BaryonCharge() == 0);

				if (tpart.BaryonCharge() == 0)
					tr1 = r2;
				else
					tr1 = r1 * pow(abs(tpart.BaryonCharge()), 1. / 3.);

				if (tpart2.BaryonCharge() == 0)
					tr2 = r2;
				else
					tr2 = r1 * pow(abs(tpart2.BaryonCharge()), 1. / 3.);
				double tb = CuteHRGHelper::brr(tr1, tr2);

				if (mode == 2 || mode == 3) {
					if (tpart.BaryonCharge() * tpart2.BaryonCharge() < 0)
						tb = 0.;
				}

				if (mode == 3) {
					if ((tpart.BaryonCharge() == 0 && tpart2.BaryonCharge() !=0) 
						|| (tpart2.BaryonCharge() == 0 && tpart.BaryonCharge() != 0))
						tb = 0.;
				}

				model.SetVirial(i, j, tb);
			}
		}
		if (mode == 4) {
			radii[i] = r1 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
		}
		if (mode == 5) {
			if (tpart.BaryonCharge() == 0)
				radii[i] = r2 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
			else
				radii[i] = r1 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
		}
		if (mode == 6) {
			if (tpart.BaryonCharge() == 0) {
				if (tpart.Strangeness() == 0)
					radii[i] = r2 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
				else
					radii[i] = r4 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
			}
			else {
				if (tpart.Strangeness() == 0)
					radii[i] = r1 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
				else
					radii[i] = r3 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
			}
		}
		if (mode == 7) {
			if (tpart.Strangeness() == 0)
				radii[i] = r1 * pow(tpart.Mass() / xMath::mnucleon(), 1. / 3.);
			else
				radii[i] = r2 * pow(tpart.Mass() / TPS->ParticleByPDG(3122).Mass(), -1. / 3.);
		}
	}

	if (!(mode >= 1 && mode <= 3))
		model.FillVirial(radii);
	model.WriteInteractionParameters(filename);
}


void QvdWHRG(std::string filename, ThermalParticleSystem *TPS, double a = 0.329, double b = 3.42) {
	ThermalModelVDWFull model(TPS);

	for (int i1 = 0; i1 < model.TPS()->Particles().size(); ++i1) {
		for (int i2 = 0; i2 < model.TPS()->Particles().size(); ++i2) {
			const ThermalParticle &part1 = model.TPS()->Particles()[i1];
			const ThermalParticle &part2 = model.TPS()->Particles()[i2];
			int B1 = part1.BaryonCharge();
			int B2 = part2.BaryonCharge();

			// Meson-meson, meson-baryon, baryon-antibaryon non-interacting
			if (!(B1 * B2 > 0)) {
				model.SetVirial(i1, i2, 0.);
				model.SetAttraction(i1, i2, 0.);
				continue;
			}

			// First excluded-volume, proportional to baryon charge if light nuclei included
			if (B1 == B2) {
				double tb = b * abs(B1);
				model.SetVirial(i1, i2, tb);
			}
			else {
				// Following prescription in nucl-th/9906068, Eqs. (53) and (54)
				double b11 = b * abs(B1);
				double b22 = b * abs(B2);
				double r1 = CuteHRGHelper::rv(b11);
				double r2 = CuteHRGHelper::rv(b22);
				double b12sym = CuteHRGHelper::brr(r1, r2);
				double b12 = 2. * b11 * b12sym / (b11 + b22);
				double b21 = 2. * b22 * b12sym / (b11 + b22);
				model.SetVirial(i1, i2, b12);
			}

			// QvdW attraction for baryon-baryon pairs only
			if ((B1 == 1 && B2 == 1) || (B1 == -1 && B2 == -1)) {
				model.SetAttraction(i1, i2, a);
			}
			else {
				model.SetAttraction(i1, i2, 0.);
			}

		}
	}

	model.WriteInteractionParameters(filename);
}

// Temperature dependence of EV-HRG thermodynamics at zero chemical potential
int main(int argc, char *argv[])
{
	string directory = string(INPUT_FOLDER) + "/list/PDG2014/";
	ThermalParticleSystem TPS(directory + "list-withnuclei.dat");
	
	// arXiv:1512.08046
	Diagonal(directory + "1512.08046.DiagonalEV.Rm00Rb03.dat", &TPS, 1, 0.3);
	Diagonal(directory + "1512.08046.DiagonalEV.Bag05.dat", &TPS, 4, 0.5);

	// arXiv:1606.06542
	Diagonal(directory + "1606.06542.DiagonalEV.sinv.dat", &TPS, 7, 0.49, 0.42);

	// arXiv:1609.03975
	QvdWHRG(directory + "1609.03975.QvdWHRG.dat", &TPS, 0.329, 3.42);

	// arXiv:1708.02852
	Crossterms(directory + "1708.02852.ImMu.dat", &TPS, 3, CuteHRGHelper::rv(1.), 0.);



	return 0;
}
