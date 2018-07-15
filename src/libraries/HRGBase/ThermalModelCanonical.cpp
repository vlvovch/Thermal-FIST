#include "HRGBase/ThermalModelCanonical.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "HRGBase/xMath.h"
#include "HRGBase/NumericalIntegration.h"

using namespace std;


ThermalModelCanonical::ThermalModelCanonical(ThermalParticleSystem *TPS_, const ThermalModelParameters& params) :
	ThermalModelBase(TPS_, params)
{

	m_TAG = "ThermalModelCanonical";

	m_Ensemble = CE;
	m_InteractionModel = Ideal;
}


ThermalModelCanonical::~ThermalModelCanonical(void)
{
}

void ThermalModelCanonical::ChangeTPS(ThermalParticleSystem *TPS_) {
	ThermalModelBase::ChangeTPS(TPS_);
}

void ThermalModelCanonical::CalculateQuantumNumbersRange(bool doubleRange)
{
	m_BMAX = 0;
	m_QMAX = 0;
	m_SMAX = 0;
	m_CMAX = 0;

	m_OnlyBose = true;

	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		ThermalParticle &part = m_TPS->Particle(i);
			if (part.Statistics()!=0 && part.CalculationType() == IdealGasFunctions::ClusterExpansion) {
				m_BMAX = max(m_BMAX, abs(part.BaryonCharge() * part.ClusterExpansionOrder()));
				m_QMAX = max(m_QMAX, abs(part.ElectricCharge() * part.ClusterExpansionOrder()));
				m_SMAX = max(m_SMAX, abs(part.Strangeness() * part.ClusterExpansionOrder()));
				m_CMAX = max(m_CMAX, abs(part.Charm() * part.ClusterExpansionOrder()));
				if (abs(part.BaryonCharge()) > 1 || (abs(part.BaryonCharge()) & 1))
					m_OnlyBose = false;
			}
			else {
				m_BMAX = max(m_BMAX, abs(part.BaryonCharge()));
				m_QMAX = max(m_QMAX, abs(part.ElectricCharge()));
				m_SMAX = max(m_SMAX, abs(part.Strangeness()));
				m_CMAX = max(m_CMAX, abs(part.Charm()));
			}
	}

	if (doubleRange) {
		m_BMAX *= 2;
		m_QMAX *= 2;
		m_SMAX *= 2;
		m_CMAX *= 2;
	}

	printf("BMAX = %d\tQMAX = %d\tSMAX = %d\tCMAX = %d\n", m_BMAX, m_QMAX, m_SMAX, m_CMAX);

	m_QNMap.clear();
	m_QNvec.resize(0);

	m_Corr.resize(0);
	m_PartialZ.resize(0);

	int ind = 0;
	for (int iB = -m_BMAX; iB <= m_BMAX; ++iB)
		for (int iQ = -m_QMAX; iQ <= m_QMAX; ++iQ)
			for (int iS = -m_SMAX; iS <= m_SMAX; ++iS) 
				for (int iC = -m_CMAX; iC <= m_CMAX; ++iC) {
				
				QuantumNumbers qn(iB, iQ, iS, iC);
				m_QNMap[qn] = ind;
				m_QNvec.push_back(qn);

				m_PartialZ.push_back(0.);
				m_Corr.push_back(1.);
				

				ind++;
			}

}

void ThermalModelCanonical::SetParameters(double T, double gammaS, double V, int B, int Q, int S, int C) {
	m_Parameters.T = T;
	m_Parameters.gammaS = gammaS;
	m_Parameters.V = V;
	m_Parameters.B = B;
	m_Parameters.Q = Q;
	m_Parameters.S = S;
	m_Parameters.C = C;
	m_Volume = V;
	m_Calculated = false;
}

void ThermalModelCanonical::SetParameters(const ThermalModelParameters& params) {
	m_Parameters = params;
	m_Volume = m_Parameters.V;
	m_Calculated = false;
}

void ThermalModelCanonical::SetStatistics(bool stats) {
	m_QuantumStats = stats;
	for (unsigned int i = 0; i < m_TPS->Particles().size(); ++i) {
			m_TPS->Particle(i).UseStatistics(stats);
	}
	CalculateQuantumNumbersRange();
}


void ThermalModelCanonical::CalculateDensities() {
	m_FluctuationsCalculated = false;

	if (m_PartialZ.size() == 0)
		CalculateQuantumNumbersRange();

	if (m_BMAX<=1)
		CalculatePartitionFunctionsBoseOnly();
	else
		CalculatePartitionFunctions();

	for (int i = 0; i < m_densities.size(); ++i) {
		ThermalParticle &tpart = m_TPS->Particle(i);
		m_densities[i] = 0.;

		if (tpart.Statistics() == 0
			|| tpart.CalculationType() != IdealGasFunctions::ClusterExpansion
			|| (tpart.BaryonCharge() != 0 && m_OnlyBose)) 
		{
			int ind = m_QNMap[QuantumNumbers(tpart.BaryonCharge(), tpart.ElectricCharge(), tpart.Strangeness(), tpart.Charm())];

			if (ind < m_Corr.size())
				m_densities[i] = m_Corr[ ind ] * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
		}
		else {
			for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
				int ind = m_QNMap[QuantumNumbers(n*tpart.BaryonCharge(), n*tpart.ElectricCharge(), n*tpart.Strangeness(), n*tpart.Charm())];
				if (ind < m_Corr.size())
					m_densities[i] += m_Corr[ ind ] * tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
			}
		}
	}

	CalculateFeeddown();

	m_Calculated = true;
	m_LastCalculationSuccessFlag = true;
}

void ThermalModelCanonical::CalculatePartitionFunctions()
{
	vector<double> Ns(m_PartialZ.size(), 0.);

	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		ThermalParticle &tpart = m_TPS->Particle(i);

		if (tpart.Statistics() == 0
			|| tpart.CalculationType() != IdealGasFunctions::ClusterExpansion) {
			int ind = m_QNMap[QuantumNumbers(tpart.BaryonCharge(), tpart.ElectricCharge(), tpart.Strangeness(), tpart.Charm())];
			if (ind < Ns.size())
				Ns[ind] += tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
		}
		else {
			for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
				int ind = m_QNMap[QuantumNumbers(n*tpart.BaryonCharge(), n*tpart.ElectricCharge(), n*tpart.Strangeness(), n*tpart.Charm())];
				if (ind < Ns.size())
					Ns[ind] += tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.) / static_cast<double>(n); // TODO: Check
			}
		}
	}

	for (int i = 0; i < Ns.size(); ++i) {
		Ns[i] *= m_Parameters.V;
		m_PartialZ[i] = 0.;
	}

	int nmax = max(3, (int)sqrt(m_Parameters.B*m_Parameters.B + m_Parameters.Q*m_Parameters.Q + m_Parameters.S*m_Parameters.S + m_Parameters.C*m_Parameters.C));
	if (m_Parameters.B == 0 && m_Parameters.Q == 0 && m_Parameters.S == 0 && m_Parameters.C == 0)
		nmax = 6;

	double dphi = xMath::Pi() / nmax;

	int nmaxB = nmax;
	double dphiB = xMath::Pi() / nmaxB;

	for (int iB = 0; iB < 2 * nmaxB; ++iB) {
		double aB = iB * dphiB;
		if (iB >= nmaxB) aB = xMath::Pi() - (iB + 1) * dphiB;
		double bB = aB + dphiB;
		vector<double> xlegB, wlegB;
		NumericalIntegration::GetCoefsIntegrateLegendre10(aB, bB, &xlegB, &wlegB);
		for (int iS = 0; iS < 2 * nmax; ++iS) {
			double aS = iS * dphi;
			if (iS >= nmax) aS = xMath::Pi() - (iS + 1) * dphi;
			double bS = aS + dphi;
			vector<double> xlegS, wlegS;
			NumericalIntegration::GetCoefsIntegrateLegendre10(aS, bS, &xlegS, &wlegS);
			for (int iQ = 0; iQ < nmax; ++iQ) {
				double aQ = iQ * dphi;
				double bQ = aQ + dphi;
				vector<double> xlegQ, wlegQ;
				NumericalIntegration::GetCoefsIntegrateLegendre10(aQ, bQ, &xlegQ, &wlegQ);

				int maxC = 2 * nmax;
				if (m_CMAX == 0)
					maxC = 1;

				for (int iC = 0; iC < maxC; ++iC) {
					vector<double> xlegC, wlegC;
					if (m_CMAX != 0) {
						double aC = iC * dphi;
						if (iC >= nmax) aC = xMath::Pi() - (iC + 1) * dphi;
						double bC = aC + dphi;
						NumericalIntegration::GetCoefsIntegrateLegendre10(aC, bC, &xlegC, &wlegC);
					}
					else {
						xlegC.resize(1);
						xlegC[0] = 0.;
						wlegC.resize(1);
						wlegC[0] = 1.;
					}

					for (int iBt = 0; iBt < xlegB.size(); ++iBt) {
						for (int iSt = 0; iSt < xlegS.size(); ++iSt) {
							for (int iQt = 0; iQt < xlegQ.size(); ++iQt) {
								for (int iCt = 0; iCt < xlegC.size(); ++iCt) {
									vector<double> cosph(m_PartialZ.size(), 0.), sinph(m_PartialZ.size(), 0.);
									double wx = 0., wy = 0., mx = 0.;
									for (int i = 0; i < m_PartialZ.size(); ++i) {
										int tB = m_QNvec[i].B;// GetB(m_BSQinv[i]);
										int tQ = m_QNvec[i].Q;// GetQ(m_BSQinv[i]);
										int tS = m_QNvec[i].S;// GetS(m_BSQinv[i]);
										int tC = m_QNvec[i].C;// GetS(m_BSQinv[i]);
										cosph[i] = cos(tB*xlegB[iBt] + tS*xlegS[iSt] + tQ*xlegQ[iQt] + tC*xlegC[iCt]);

										mx += Ns[i] * cosph[i];
									}
									for (int iN = 0; iN < m_PartialZ.size(); ++iN) {
										int tBg = m_Parameters.B - m_QNvec[iN].B;
										int tQg = m_Parameters.Q - m_QNvec[iN].Q;
										int tSg = m_Parameters.S - m_QNvec[iN].S;
										int tCg = m_Parameters.C - m_QNvec[iN].C;
										m_PartialZ[iN] += wlegB[iBt] * wlegS[iSt] * wlegQ[iQt] * wlegC[iCt] * cos(tBg*xlegB[iBt] + tSg*xlegS[iSt] + tQg*xlegQ[iQt] + tCg*xlegC[iCt]) * exp(mx);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	for (int iN = 0; iN < m_PartialZ.size(); ++iN) {
		m_PartialZ[iN] /= pow(2. * xMath::Pi(), 3);
		m_PartialZ[iN] *= 2.; // Because phiQ was integrated over [0,\pi] only 
		if (m_CMAX != 0)
			m_PartialZ[iN] /= 2. * xMath::Pi();
	}


	m_Corr.resize(m_PartialZ.size());
	for (int iN = 0; iN < m_PartialZ.size(); ++iN) {
		m_Corr[iN] = m_PartialZ[iN] / m_PartialZ[m_QNMap[QuantumNumbers(0, 0, 0, 0)]];
	}
}

void ThermalModelCanonical::CalculatePartitionFunctionsBoseOnly()
{
	vector<double> Ns(m_PartialZ.size(), 0.);

	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		ThermalParticle &tpart = m_TPS->Particle(i);
		if (abs(tpart.BaryonCharge()) > 1) continue;

		if (tpart.Statistics() == 0
			|| tpart.CalculationType() != IdealGasFunctions::ClusterExpansion
			|| tpart.BaryonCharge() != 0) {
			int ind = m_QNMap[QuantumNumbers(tpart.BaryonCharge(), tpart.ElectricCharge(), tpart.Strangeness(), tpart.Charm())];
			if (ind < Ns.size())
				Ns[ind] += tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
		}
		else {
			for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
				int ind = m_QNMap[QuantumNumbers(n*tpart.BaryonCharge(), n*tpart.ElectricCharge(), n*tpart.Strangeness(), n*tpart.Charm())];
				if (ind < Ns.size())
					Ns[ind] += tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.) / static_cast<double>(n); // TODO: Check
			}
		}
	}

	for (int i = 0; i < Ns.size(); ++i) {
		Ns[i] *= m_Parameters.V;
		m_PartialZ[i] = 0.;
	}

	int nmax = max(3, (int)sqrt(m_Parameters.B*m_Parameters.B + m_Parameters.Q*m_Parameters.Q + m_Parameters.S*m_Parameters.S + m_Parameters.C*m_Parameters.C));
	if (m_Parameters.B == 0 && m_Parameters.Q == 0 && m_Parameters.S == 0 && m_Parameters.C == 0)
		nmax = 6;

	double dphi = xMath::Pi() / nmax;
	for (int iS = 0; iS < 2 * nmax; ++iS) {
		double aS = iS * dphi;
		if (iS >= nmax) aS = xMath::Pi() - (iS + 1) * dphi;
		double bS = aS + dphi;
		vector<double> xlegS, wlegS;
		NumericalIntegration::GetCoefsIntegrateLegendre10(aS, bS, &xlegS, &wlegS);
		for (int iQ = 0; iQ < nmax; ++iQ) {
			double aQ = iQ * dphi;
			double bQ = aQ + dphi;
			vector<double> xlegQ, wlegQ;
			NumericalIntegration::GetCoefsIntegrateLegendre10(aQ, bQ, &xlegQ, &wlegQ);

			int maxC = 2 * nmax;
			if (m_CMAX == 0)
				maxC = 1;

			for (int iC = 0; iC < maxC; ++iC) {
				vector<double> xlegC, wlegC;
				if (m_CMAX != 0) {
					double aC = iC * dphi;
					if (iC >= nmax) aC = xMath::Pi() - (iC + 1) * dphi;
					double bC = aC + dphi;
					NumericalIntegration::GetCoefsIntegrateLegendre10(aC, bC, &xlegC, &wlegC);
				}
				else {
					xlegC.resize(1);
					xlegC[0] = 0.;
					wlegC.resize(1);
					wlegC[0] = 1.;
				}


				for (int iSt = 0; iSt < xlegS.size(); ++iSt) {
					for (int iQt = 0; iQt < xlegQ.size(); ++iQt) {
						for (int iCt = 0; iCt < xlegC.size(); ++iCt) {
							vector<double> cosph(m_PartialZ.size(), 0.), sinph(m_PartialZ.size(), 0.);
							double wx = 0., wy = 0., mx = 0.;
							for (int i = 0; i < m_PartialZ.size(); ++i) {
								int tQ = m_QNvec[i].Q;// GetQ(m_BSQinv[i]);
								int tS = m_QNvec[i].S;// GetS(m_BSQinv[i]);
								int tC = m_QNvec[i].C;// GetS(m_BSQinv[i]);
								cosph[i] = cos(tS*xlegS[iSt] + tQ*xlegQ[iQt] + tC*xlegC[iCt]);
								sinph[i] = sin(tS*xlegS[iSt] + tQ*xlegQ[iQt] + tC*xlegC[iCt]);

								if (m_QNvec[i].B == 1) {
									wx += Ns[i] * cosph[i];
									wy += Ns[i] * sinph[i];
								}
								else if (m_QNvec[i].B == 0) {
									mx += Ns[i] * cosph[i];
								}
							}
							double wmod = sqrt(wx*wx + wy*wy);
							double warg = atan2(wy, wx);
							for (int iN = 0; iN < m_PartialZ.size(); ++iN) {
								int tBg = m_Parameters.B - m_QNvec[iN].B;
								int tQg = m_Parameters.Q - m_QNvec[iN].Q;
								int tSg = m_Parameters.S - m_QNvec[iN].S;
								int tCg = m_Parameters.C - m_QNvec[iN].C;
								m_PartialZ[iN] += wlegS[iSt] * wlegQ[iQt] * wlegC[iCt] * cos(tSg*xlegS[iSt] + tQg*xlegQ[iQt] + tCg*xlegC[iCt] - tBg * warg) * exp(mx) * xMath::BesselI(tBg, 2. * wmod);
							}
						}
					}
				}
			}
		}
	}

	for (int iN = 0; iN < m_PartialZ.size(); ++iN) {
		m_PartialZ[iN] /= pow(2. * xMath::Pi(), 2);// *2.;
		m_PartialZ[iN] *= 2.; // Because phiQ was integrated over [0,\pi] only 
		if (m_CMAX != 0)
			m_PartialZ[iN] /= 2. * xMath::Pi();
	}


	m_Corr.resize(m_PartialZ.size());
	for (int iN = 0; iN < m_PartialZ.size(); ++iN) {
		m_Corr[iN] = m_PartialZ[iN] / m_PartialZ[m_QNMap[QuantumNumbers(0, 0, 0, 0)]];
	}
}


double ThermalModelCanonical::CalculateParticleScaledVariance(int part)
{
	ThermalParticle &tpart = m_TPS->Particle(part);
	double ret1 = 0., ret2 = 0., ret3 = 0.;
	

	if (tpart.Statistics() == 0
		|| tpart.CalculationType() != IdealGasFunctions::ClusterExpansion
		|| (tpart.BaryonCharge() != 0 && m_OnlyBose))
	{
		int ind  = m_QNMap[QuantumNumbers(tpart.BaryonCharge(), tpart.ElectricCharge(), tpart.Strangeness(), tpart.Charm())];
		int ind2 = m_QNMap[QuantumNumbers(2 * tpart.BaryonCharge(), 2 * tpart.ElectricCharge(), 2 * tpart.Strangeness(), 2 * tpart.Charm())];

		ret1 = 1.;
		if (ind < m_Corr.size() && ind2 < m_Corr.size())
			ret2 = m_Corr[ind2] / m_Corr[ind] * m_Parameters.V * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
			
		if (ind < m_Corr.size())
			ret3 = -m_Corr[ind] * m_Parameters.V * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
	}
	else {
		double ret1num = 0., ret1zn = 0.;
		for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
			int ind = m_QNMap[QuantumNumbers(n*tpart.BaryonCharge(), n*tpart.ElectricCharge(), n*tpart.Strangeness(), n*tpart.Charm())];
			
			double densityClusterN = tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);

			if (ind < m_Corr.size()) {
				ret1num += m_Corr[ind] * n * densityClusterN;
				ret1zn  += m_Corr[ind] * densityClusterN;
			}

			for (int n2 = 1; n2 <= tpart.ClusterExpansionOrder(); ++n2) {
				if (m_QNMap.count(QuantumNumbers((n + n2)*tpart.BaryonCharge(), (n + n2)*tpart.ElectricCharge(), (n + n2)*tpart.Strangeness(), (n + n2)*tpart.Charm())) != 0) {
					int ind2 = m_QNMap[QuantumNumbers((n + n2)*tpart.BaryonCharge(), (n + n2)*tpart.ElectricCharge(), (n + n2)*tpart.Strangeness(), (n + n2)*tpart.Charm())];
					if (ind < m_Corr.size() && ind2 < m_Corr.size())
						ret2 += densityClusterN * m_Corr[ind2] * m_Parameters.V * tpart.DensityCluster(n2, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
				}
			}
		}
		ret1 = ret1num / ret1zn;
		ret2 = ret2 / ret1zn;
		ret3 = -ret1zn * m_Parameters.V;
	}
	return ret1 + ret2 + ret3;
}

void ThermalModelCanonical::CalculateTwoParticleCorrelations() {
	int NN = m_densities.size();

	vector<double> yld(NN, 0);
	vector<double> ret1num(NN, 0);
	vector< vector<double> > ret2num(NN, vector<double>(NN, 0.));

	for (int i = 0; i < NN; ++i)
		yld[i] = m_densities[i] * m_Parameters.V;

	for (int i = 0; i < NN; ++i) {
		ThermalParticle &tpart = m_TPS->Particle(i);
		if (tpart.Statistics() == 0
			|| tpart.CalculationType() != IdealGasFunctions::ClusterExpansion
			|| (tpart.BaryonCharge() != 0 && m_OnlyBose))
		{
			ret1num[i] = yld[i];
		}
		else {
			for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
				int ind = m_QNMap[QuantumNumbers(n*tpart.BaryonCharge(), n*tpart.ElectricCharge(), n*tpart.Strangeness(), n*tpart.Charm())];

				double densityClusterN = tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);

				if (ind < m_Corr.size())
					ret1num[i] += m_Corr[ind] * n * densityClusterN * m_Parameters.V;
			}
		}
	}

	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			ThermalParticle &tpart1 = m_TPS->Particle(i);
			ThermalParticle &tpart2 = m_TPS->Particle(j);

			int n1max = tpart1.ClusterExpansionOrder();
			int n2max = tpart2.ClusterExpansionOrder();

			if (tpart1.Statistics() == 0 || tpart1.CalculationType() != IdealGasFunctions::ClusterExpansion
				|| (tpart1.BaryonCharge() != 0 && m_OnlyBose))
				n1max = 1;
			if (tpart2.Statistics() == 0 || tpart2.CalculationType() != IdealGasFunctions::ClusterExpansion
				|| (tpart2.BaryonCharge() != 0 && m_OnlyBose))
				n2max = 1;


			for (int n1 = 1; n1 <= n1max; ++n1) {
				for (int n2 = 1; n2 <= n2max; ++n2) {
					int ind = m_QNMap[QuantumNumbers(
						n1*tpart1.BaryonCharge() + n2*tpart2.BaryonCharge(),
						n1*tpart1.ElectricCharge() + n2*tpart2.ElectricCharge(),
						n1*tpart1.Strangeness() + n2*tpart2.Strangeness(),
						n1*tpart1.Charm() + n2*tpart2.Charm())];

					double densityClusterN1 = tpart1.DensityCluster(n1, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
					double densityClusterN2 = tpart2.DensityCluster(n2, m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);

					if (ind < m_Corr.size())
						ret2num[i][j] += m_Corr[ind] * densityClusterN1 * densityClusterN2 * m_Parameters.V * m_Parameters.V;
				}
			}
		}
	}


	m_PrimCorrel.resize(NN);
	for (int i = 0; i < NN; ++i) 
		m_PrimCorrel[i].resize(NN);
	m_TotalCorrel = m_PrimCorrel;

	for (int i = 0; i < NN; ++i) {
		for (int j = 0; j < NN; ++j) {
			m_PrimCorrel[i][j] = 0.;
			if (i == j)
				m_PrimCorrel[i][j] += ret1num[i] / yld[i];
			m_PrimCorrel[i][j] += ret2num[i][j] / yld[i];
			m_PrimCorrel[i][j] += -yld[j];
			m_PrimCorrel[i][j] *= yld[i] / m_Parameters.V / m_Parameters.T;
			if (yld[i] == 0.0)
				m_PrimCorrel[i][j] = 0.;
		}
	}

	for (int i = 0; i < NN; ++i) {
		m_wprim[i] = m_PrimCorrel[i][i];
		if (m_densities[i] > 0.) m_wprim[i] *= m_Parameters.T / m_densities[i];
		else m_wprim[i] = 1.;
	}

	CalculateSusceptibilityMatrix();
	CalculateTwoParticleFluctuationsDecays();

}

void ThermalModelCanonical::CalculateFluctuations() {
	CalculateTwoParticleCorrelations();

	m_FluctuationsCalculated = true;

	for (int i = 0; i < m_wprim.size(); ++i) {
		m_skewprim[i] = 1.;
		m_kurtprim[i] = 1.;
		m_skewtot[i]  = 1.;
		m_kurttot[i]  = 1.;
	}
}

double ThermalModelCanonical::CalculateEnergyDensity() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;

	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		ThermalParticle &tpart = m_TPS->Particle(i);
		{
			if (tpart.Statistics() == 0
				|| tpart.CalculationType() != IdealGasFunctions::ClusterExpansion
				|| (tpart.BaryonCharge() != 0 && m_OnlyBose)) {
				int ind = m_QNMap[QuantumNumbers(tpart.BaryonCharge(), tpart.ElectricCharge(), tpart.Strangeness(), tpart.Charm())];

				if (ind < m_Corr.size())
					ret += m_Corr[ ind ] * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, 0.);
			}
			else {
				for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
					int ind = m_QNMap[QuantumNumbers(n*tpart.BaryonCharge(), n*tpart.ElectricCharge(), n*tpart.Strangeness(), n*tpart.Charm())];
					if (ind < m_Corr.size())
						ret += m_Corr[ ind ] * tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::EnergyDensity, m_UseWidth, 0.);
				}
			}
		}
	}

	return ret;
}

double ThermalModelCanonical::CalculatePressure() {
	if (!m_Calculated) CalculateDensities();
	double ret = 0.;

	for (int i = 0; i < m_TPS->Particles().size(); ++i) {
		ThermalParticle &tpart = m_TPS->Particle(i);
		{
			if (tpart.Statistics() == 0
				|| tpart.CalculationType() != IdealGasFunctions::ClusterExpansion
				|| (tpart.BaryonCharge() != 0 && m_OnlyBose)) {
				int ind = m_QNMap[QuantumNumbers(tpart.BaryonCharge(), tpart.ElectricCharge(), tpart.Strangeness(), tpart.Charm())];

				if (ind < m_Corr.size())
					ret += m_Corr[ind] * tpart.DensityCluster(1, m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, 0.);
			}
			else {
				for (int n = 1; n <= tpart.ClusterExpansionOrder(); ++n) {
					int ind = m_QNMap[QuantumNumbers(n*tpart.BaryonCharge(), n*tpart.ElectricCharge(), n*tpart.Strangeness(), n*tpart.Charm())];

					if (ind < m_Corr.size())
						ret += m_Corr[ind] * tpart.DensityCluster(n, m_Parameters, IdealGasFunctions::Pressure, m_UseWidth, 0.);
				}
			}
		}
	}

	return ret;
}

double ThermalModelCanonical::CalculateEntropyDensity()
{
	return (CalculateEnergyDensity() / m_Parameters.T) + log(m_PartialZ[m_QNMap[QuantumNumbers(0, 0, 0, 0)]]) / m_Parameters.V;
}

double ThermalModelCanonical::GetGCEDensity(int i) const
{
	return m_TPS->Particles()[i].Density(m_Parameters, IdealGasFunctions::ParticleDensity, m_UseWidth, 0.);
}

