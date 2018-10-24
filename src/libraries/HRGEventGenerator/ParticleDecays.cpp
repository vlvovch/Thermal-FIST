#include "HRGEventGenerator/ParticleDecays.h"

#include "HRGBase/xMath.h"
#include "HRGEventGenerator/RandomGenerators.h"


namespace ParticleDecays {

	// Square of probability density for m12 in 3-body decay
	double ThreeBodym12F2(double m12, double M, double m1, double m2, double m3) {
		double m12k = m12 * m12;
		return (m12k - (m1+m2)*(m1+m2)) * (m12k - (m1-m2)*(m1-m2)) * (M*M - (m12+m3)*(m12+m3)) * (M*M - (m12-m3)*(m12-m3)) / m12k;
	}

	// Ternary search for maximum of probability density for m12
	double TernaryThreeBodym12Maximum(double M, double m1_, double m2_, double m3_) {
		double l = m1_+m2_, r = M-m3_;
		double m1 = 0., m2 = 0.;
		m1 = l + (r-l) / 3.;
		m2 = r - (r-l) / 3.;
		while ((m2-m1)>1e-8) {
			double f1 = ThreeBodym12F2(m1, M, m1_, m2_, m3_);
			double f2 = ThreeBodym12F2(m2, M, m1_, m2_, m3_);
			if (f1<f2) {
				l = m1;
			}
			else {
				r = m2;
			}
			m1 = l + (r-l) / 3.;
			m2 = r - (r-l) / 3.;
		}
		return sqrt(ThreeBodym12F2((m1+m2)/2., M, m1_, m2_, m3_));
	}

	int threebodysucc = 0, threebodytot = 0;
	
	// Random sample for m12 in a 3-body decay
	double GetRandomThreeBodym12(double M, double m1_, double m2_, double m3_, double fm12max) {
		while (true) {
			threebodytot++;
			double x0 = m1_ + m2_ + (M - m1_ - m2_ - m3_) * RandomGenerators::randgenMT.randDblExc();
			double y0 = fm12max * RandomGenerators::randgenMT.randDblExc();
			if (y0*y0<ThreeBodym12F2(x0, M, m1_, m2_, m3_)) { 
				threebodysucc++; 
				return x0; 
			}
		}
		return 0.;
	}




	SimpleParticle LorentzBoost(const SimpleParticle &part, double vx, double vy, double vz) {
		SimpleParticle ret = part;
		double v2 = vx*vx + vy*vy + vz*vz;
		double gamma = 1./sqrt(1.-v2);
		ret.p0 = gamma * part.p0 - gamma * (vx * part.px + vy * part.py + vz * part.pz);
		ret.px = -gamma * vx * part.p0 + (1. + (gamma-1.)*vx*vx/v2) * part.px +
			(gamma-1.)*vx*vy/v2 * part.py + (gamma-1.)*vx*vz/v2 * part.pz;
		ret.py = -gamma * vy * part.p0 + (1. + (gamma-1.)*vy*vy/v2) * part.py +
			(gamma-1.)*vy*vx/v2 * part.px + (gamma-1.)*vy*vz/v2 * part.pz;
		ret.pz = -gamma * vz * part.p0 + (1. + (gamma-1.)*vz*vz/v2) * part.pz +
			(gamma-1.)*vz*vx/v2 * part.px + (gamma-1.)*vz*vy/v2 * part.py;
		return ret;
	}

	std::vector<SimpleParticle> TwoBodyDecay(const SimpleParticle & Mother, double m1, int pdg1, double m2, int pdg2) {
		std::vector<SimpleParticle> ret(0);
		ret.push_back(Mother);
		ret.push_back(Mother);
		ret[0].PDGID = pdg1;
		ret[0].m = m1;
		ret[1].PDGID = pdg2;
		ret[1].m = m2;

		double vx = Mother.px / Mother.p0;
		double vy = Mother.py / Mother.p0;
		double vz = Mother.pz / Mother.p0;
		SimpleParticle Mo = LorentzBoost(Mother, vx, vy, vz);
		double ten1 = (Mo.m*Mo.m - m2*m2 + m1*m1) / 2./ Mo.m;
		double tp   = sqrt(ten1*ten1 - m1*m1);
		double tphi = 2. * xMath::Pi() * RandomGenerators::randgenMT.rand();
		double cthe = 2. * RandomGenerators::randgenMT.rand() - 1.;
		double sthe = sqrt(1.- cthe*cthe);
		ret[0].px = tp * cos(tphi) * sthe;
		ret[0].py = tp * sin(tphi) * sthe;
		ret[0].pz = tp * cthe;
		ret[0].p0 = ten1;
		ret[1].px = -ret[0].px;
		ret[1].py = -ret[0].py;
		ret[1].pz = -ret[0].pz;
		ret[1].p0 = sqrt(m2*m2 + ret[1].px*ret[1].px + ret[1].py*ret[1].py + ret[1].pz*ret[1].pz);

		ret[0] = LorentzBoost(ret[0], -vx, -vy, -vz);
		ret[1] = LorentzBoost(ret[1], -vx, -vy, -vz);

		ret[0].MotherPDGID = Mother.PDGID;
		ret[1].MotherPDGID = Mother.PDGID;

#ifdef DEBUGDECAYS
		for(int i=0;i<ret.size();++i)
			if (ret[i].px!=ret[i].px) {
				std::cout << "2BodyBug!\n";
			}
#endif

		return ret;
	}

	std::vector<SimpleParticle> ManyBodyDecay(const SimpleParticle & Mother, std::vector<double> masses, std::vector<int> pdgs) {
		std::vector<SimpleParticle> ret(0);
		if (masses.size()<1) return ret;

		// If only one daughter listed, assume a radiative decay A -> B + gamma
		if (masses.size() == 1)
		{
			masses.push_back(0.);
			pdgs.push_back(22);
		}
		SimpleParticle Mother2 = Mother;
		// Mass validation
		double tmasssum = 0.;
		for(int i=0;i<masses.size();++i) tmasssum += masses[i];
		if (Mother2.m < tmasssum) Mother2.m = tmasssum + 1e-7;

		if (masses.size()==2) return TwoBodyDecay(Mother2, masses[0], pdgs[0], masses[1], pdgs[1]);
		double tmin = 0.;
		for(int i=0;i<masses.size()-1;++i) tmin += masses[i];
		double tmax = Mother2.m - masses[masses.size()-1];
		double mijk = 0.;
		if (masses.size()==3) {
			mijk = GetRandomThreeBodym12(Mother2.m, masses[0], masses[1], masses[2], 1.01*TernaryThreeBodym12Maximum(Mother2.m, masses[0], masses[1], masses[2]));
		}
		else // More than 3 body decay kinematics are only approximate!
		{
			mijk = tmin + (tmax-tmin) * RandomGenerators::randgenMT.rand();
		}
		std::vector<SimpleParticle> ret1 = TwoBodyDecay(Mother2, mijk, 11111, masses[masses.size()-1], pdgs[pdgs.size()-1]);
		ret.push_back(ret1[1]);
		masses.resize(masses.size()-1);
		pdgs.resize(pdgs.size()-1);
		ret1 = ManyBodyDecay(ret1[0], masses, pdgs);
		for(int i=0;i<ret1.size();++i)
			ret.push_back(ret1[i]);

		for (int i = 0; i < ret.size(); ++i)
			ret[i].MotherPDGID = Mother.PDGID;

#ifdef DEBUGDECAYS
		for(int i=0;i<ret.size();++i)
			if (ret[i].px!=ret[i].px) std::cout << "**WARNING** NaN in NBodyDecay procedure output!\n";
#endif

		return ret;
	}
	
}
