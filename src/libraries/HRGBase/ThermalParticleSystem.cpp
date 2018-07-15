#include "HRGBase/ThermalParticleSystem.h"

#include <fstream>
#include <algorithm>
#include <queue>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>
#include <cmath>
#include <cstdlib>

using namespace std;

namespace /*ThermalParticleSystemNameSpace*/ {
	bool cmpParticleMass(const ThermalParticle &a, const ThermalParticle &b) {
		return (a.Mass() < b.Mass());
	}
}

ThermalParticleSystem::ThermalParticleSystem(std::string InputFile, bool GenAntiP, double mcut)
{
	m_NumberOfParticles = 0;
	m_Particles.resize(0);
	m_PDGtoID.clear();
	m_IDtoPDG.resize(0);

	LoadTable(InputFile, GenAntiP, mcut);
}


ThermalParticleSystem::~ThermalParticleSystem(void)
{
}

std::vector<ParticleDecay> ThermalParticleSystem::GetDecaysFromAntiParticle(const std::vector<ParticleDecay> &Decays) {
	std::vector<ParticleDecay> ret = Decays;
	for (unsigned int i = 0; i < ret.size(); ++i) {
		for (unsigned int j = 0; j < ret[i].mDaughters.size(); ++j) {
			if (m_PDGtoID.count(-ret[i].mDaughters[j]) > 0) ret[i].mDaughters[j] = -ret[i].mDaughters[j];
		}
	}
	return ret;
}

void ThermalParticleSystem::FillDecayProperties()
{
	for (int i = 0; i < m_Particles.size(); ++i) {
		if (m_Particles[i].Decays().size() != 0) {
			double tsumb = 0.;

			for (int j = 0; j < m_Particles[i].Decays().size(); ++j) {

				m_Particles[i].Decays()[j].mPole = m_Particles[i].Mass();

				double M0 = 0.;
				double tS = 0.;
				for (int k = 0; k < m_Particles[i].Decays()[j].mDaughters.size(); ++k) {
					int tid = PdgToId(m_Particles[i].Decays()[j].mDaughters[k]);
					if (tid != -1) {
						M0 += m_Particles[tid].Mass();
						tS += max(0., (m_Particles[tid].Degeneracy() - 1.) / 2.);
					}
				}
				m_Particles[i].Decays()[j].mM0 = M0;
				m_Particles[i].Decays()[j].mL = abs(max(0., (m_Particles[i].Degeneracy() - 1.) / 2.) - tS);

				//printf("%15s %d: %lf\n", m_Particles[i].Name().c_str(), j, m_Particles[i].Decays()[j].mL);

				//tsumb += m_Particles[i].Decays()[j].mBratioAtPole;
				tsumb += m_Particles[i].Decays()[j].mBratio;
				m_Particles[i].Decays()[j].mBratioAverage = m_Particles[i].Decays()[j].mBratio;
			}
		}
		m_Particles[i].CalculateAndSetDynamicalThreshold();
		m_Particles[i].FillCoefficientsDynamical();
	}
}

void ThermalParticleSystem::FillDecayThresholds()
{
	for (int i = 0; i < m_Particles.size(); ++i) {
		if (m_Particles[i].Decays().size() != 0) {
			for (int j = 0; j < m_Particles[i].Decays().size(); ++j) {
				double M0 = 0.;
				for (int k = 0; k < m_Particles[i].Decays()[j].mDaughters.size(); ++k) {
					if (PdgToId(m_Particles[i].Decays()[j].mDaughters[k]) != -1)
						M0 += m_Particles[PdgToId(m_Particles[i].Decays()[j].mDaughters[k])].Mass();
				}
				m_Particles[i].Decays()[j].mM0 = M0;
			}
			m_Particles[i].FillCoefficients();
		}
	}
}

void ThermalParticleSystem::FillResonanceDecays() {
	for (int i = 0; i < m_Particles.size(); ++i) {
		m_Particles[i].DecayContributions().resize(0);
		m_Particles[i].DecayContributionsSigmas().resize(0);
		m_Particles[i].DecayProbabilities().resize(0);
		m_Particles[i].DecayCumulants().resize(0);
	}
	for (int i = m_Particles.size() - 1; i >= 0; i--)
		if (!m_Particles[i].IsStable()) {
			GoResonance(i, i, 1.);
		}

	for (int i = 0; i < m_Particles.size(); ++i) {
		for (int j = 0; j < m_Particles[i].DecayContributions().size(); ++j) {
			vector<double> tmp = GoResonanceDecayProbs(m_Particles[i].DecayContributions()[j].second, i, true);
			if (tmp.size() > 1) m_Particles[i].DecayProbabilities().push_back(make_pair(tmp, m_Particles[i].DecayContributions()[j].second));
		}
		for (int j = 0; j < m_Particles[i].DecayProbabilities().size(); ++j) {
			double tmp = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0.;
			for (int jj = 0; jj < m_Particles[i].DecayProbabilities()[j].first.size(); ++jj) {
				tmp += m_Particles[i].DecayProbabilities()[j].first[jj] * jj;
				tmp2 += m_Particles[i].DecayProbabilities()[j].first[jj] * jj * jj;
				tmp3 += m_Particles[i].DecayProbabilities()[j].first[jj] * jj * jj * jj;
				tmp4 += m_Particles[i].DecayProbabilities()[j].first[jj] * jj * jj * jj * jj;
			}
			double n2 = 0., n3 = 0., n4 = 0.;
			n2 = tmp2 - tmp*tmp;
			n3 = tmp3 - 3. * tmp2 * tmp + 2. * tmp * tmp * tmp;
			n4 = tmp4 - 4. * tmp3 * tmp + 6. * tmp2 * tmp * tmp - 3. * tmp * tmp * tmp * tmp - 3. * n2 * n2;
			m_Particles[i].DecayContributionsSigmas().push_back(make_pair(tmp2 - tmp*tmp, m_Particles[i].DecayProbabilities()[j].second));
			vector<double> moments(0);
			moments.push_back(tmp);
			moments.push_back(n2);
			moments.push_back(n3);
			moments.push_back(n4);
			m_Particles[i].DecayCumulants().push_back(make_pair(moments, m_Particles[i].DecayProbabilities()[j].second));
		}
	}

	for (int i = 0; i < m_Particles.size(); ++i) {
		vector<int> nchtyp(0);
		nchtyp.push_back(0);
		nchtyp.push_back(1);
		nchtyp.push_back(-1);

		m_Particles[i].Nch().resize(0);
		m_Particles[i].DeltaNch().resize(0);

		for (int nti = 0; nti < 3; nti++) {
			vector<double> prob = GoResonanceDecayProbsCharge(i, nchtyp[nti], true);
			double tmp = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0.;
			for (int jj = 0; jj < prob.size(); ++jj) {
				tmp  += prob[jj] * jj;
				tmp2 += prob[jj] * jj * jj;
				tmp3 += prob[jj] * jj * jj * jj;
				tmp4 += prob[jj] * jj * jj * jj * jj;
			}
			double n2 = 0., n3 = 0., n4 = 0.;
			n2 = tmp2 - tmp*tmp;
			n3 = tmp3 - 3. * tmp2 * tmp + 2. * tmp * tmp * tmp;
			n4 = tmp4 - 4. * tmp3 * tmp + 6. * tmp2 * tmp * tmp - 3. * tmp * tmp * tmp * tmp - 3. * n2 * n2;
			m_Particles[i].Nch().push_back(tmp);
			m_Particles[i].DeltaNch().push_back(n2);
		}
		
	}
}

void ThermalParticleSystem::GoResonance(int ind, int startind, double BR) {
	if (ind != startind && m_Particles[ind].DecayContributions().size() > 0 && m_Particles[ind].DecayContributions()[m_Particles[ind].DecayContributions().size() - 1].second == startind)
	{
		m_Particles[ind].DecayContributions()[m_Particles[ind].DecayContributions().size() - 1].first += BR;
	}
	else if (ind != startind) m_Particles[ind].DecayContributions().push_back(make_pair(BR, startind));
	if (!m_Particles[ind].IsStable()) {
		for (int i = 0; i < m_Particles[ind].Decays().size(); ++i) {
			double tbr = m_Particles[ind].Decays()[i].mBratio;

			// TODO: Fix(?) for canonical ensemble
			if (m_ResonanceWidthIntegrationType == ThermalParticle::eBW && ind == startind)
				tbr = m_Particles[ind].Decays()[i].mBratioAverage;

			for (int j = 0; j < m_Particles[ind].Decays()[i].mDaughters.size(); ++j) {
				if (m_PDGtoID.count(m_Particles[ind].Decays()[i].mDaughters[j]) != 0) 
					GoResonance(m_PDGtoID[m_Particles[ind].Decays()[i].mDaughters[j]], startind, BR*tbr);
			}
		}
	}
}

std::vector<double> ThermalParticleSystem::GoResonanceDecayProbs(int ind, int goalind, bool firstdecay) {
	std::vector<double> ret(1, 0.);
	if (m_Particles[ind].IsStable()) {
		if (ind == goalind) ret.push_back(1.);
		else ret[0] = 1.;
		return ret;
	}
	else if (ind == goalind) {
		ret.push_back(1.);
		return ret;
	}
	else {
		ret[0] = 0.;
		vector<double> tret;
		for (int i = 0; i < m_Particles[ind].Decays().size(); ++i) {
			double tbr = m_Particles[ind].Decays()[i].mBratio;
			if (m_ResonanceWidthIntegrationType == ThermalParticle::eBW && firstdecay)
				tbr = m_Particles[ind].Decays()[i].mBratioAverage;

			tret.resize(1);
			tret[0] = 1.;
			for (int j = 0; j < m_Particles[ind].Decays()[i].mDaughters.size(); ++j) {
				if (m_PDGtoID.count(m_Particles[ind].Decays()[i].mDaughters[j]) != 0) {
					vector<double> tmp = GoResonanceDecayProbs(m_PDGtoID[m_Particles[ind].Decays()[i].mDaughters[j]], goalind);
					vector<double> tmp2(tret.size() + tmp.size() - 1, 0.);
					for (int i1 = 0; i1 < tret.size(); ++i1)
						for (int i2 = 0; i2 < tmp.size(); ++i2)
							tmp2[i1 + i2] += tret[i1] * tmp[i2];
					tret = tmp2;
				}
			}
			if (ret.size() < tret.size()) ret.resize(tret.size(), 0.);
			for (int j = 0; j < tret.size(); ++j)
				//ret[j] += m_Particles[ind].Decays()[i].mBratio * tret[j];
				ret[j] += tbr * tret[j];
		}
		double totprob = 0.;
		for (int i = 0; i < ret.size(); ++i)
			totprob += ret[i];
		if (totprob > 1.) {
			for (int i = 0; i < ret.size(); ++i)
				ret[i] *= 1. / totprob;
		}
		else {
			ret[0] += 1. - totprob;
		}
		return ret;
	}
	return ret;
}

std::vector<double> ThermalParticleSystem::GoResonanceDecayProbsCharge(int ind, int nch, bool firstdecay)
{
	bool fl = false;
	int tQ = m_Particles[ind].ElectricCharge();
	if (nch == 0 && tQ != 0)
		fl = true;
	if (nch == 1 && tQ > 0)
		fl = true;
	if (nch == -1 && tQ < 0)
		fl = true;

	std::vector<double> ret(1, 0.);
	if (m_Particles[ind].IsStable()) {
		if (fl) ret.push_back(1.);
		else ret[0] = 1.;
		return ret;
	}
	else {
		ret[0] = 0.;
		vector<double> tret;
		for (int i = 0; i < m_Particles[ind].Decays().size(); ++i) {
			double tbr = m_Particles[ind].Decays()[i].mBratio;
			if (m_ResonanceWidthIntegrationType == ThermalParticle::eBW && firstdecay)
				tbr = m_Particles[ind].Decays()[i].mBratioAverage;

			tret.resize(1);
			tret[0] = 1.;
			for (int j = 0; j < m_Particles[ind].Decays()[i].mDaughters.size(); ++j) {
				if (m_PDGtoID.count(m_Particles[ind].Decays()[i].mDaughters[j]) != 0) {
					vector<double> tmp = GoResonanceDecayProbsCharge(m_PDGtoID[m_Particles[ind].Decays()[i].mDaughters[j]], nch);
					vector<double> tmp2(tret.size() + tmp.size() - 1, 0.);
					for (int i1 = 0; i1 < tret.size(); ++i1)
						for (int i2 = 0; i2 < tmp.size(); ++i2)
							tmp2[i1 + i2] += tret[i1] * tmp[i2];
					tret = tmp2;
				}
			}
			if (ret.size() < tret.size()) 
				ret.resize(tret.size(), 0.);
			for (int j = 0; j < tret.size(); ++j)
				//ret[j] += m_Particles[ind].Decays()[i].mBratio * tret[j];
				ret[j] += tbr * tret[j];
		}
		double totprob = 0.;
		for (int i = 0; i < ret.size(); ++i)
			totprob += ret[i];
		if (totprob > 1.) {
			for (int i = 0; i < ret.size(); ++i)
				ret[i] *= 1. / totprob;
		}
		else {
			ret[0] += 1. - totprob;
		}
		return ret;
	}
	return ret;
}


namespace CuteHRGHelper {
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim)) {
			elems.push_back(item);
		}
		return elems;
	}

	std::vector<std::string> split(const std::string &s, char delim) {
			std::vector<std::string> elems;
			split(s, delim, elems);
			return elems;
		}
}

void ThermalParticleSystem::LoadTable(std::string InputFile, bool GenAntiP, double mcut) {
	
	m_NumberOfParticles = 0;
	m_Particles.resize(0);
	m_PDGtoID.clear();
	m_IDtoPDG.resize(0);

	m_NumBaryons = m_NumCharged = m_NumStrange = m_NumCharmed = 0;

	ifstream fin;
	fin.open(InputFile.c_str());
	if (fin.is_open()) {

		char tmpc[2000];
		fin.getline(tmpc, 2000);
		string tmp = string(tmpc);
		vector<string> elems = CuteHRGHelper::split(tmp, '#');

		int flnew = 0;
		if (elems.size() < 1 || CuteHRGHelper::split(elems[0], ';').size() < 4)
			flnew = 1;
		else
			flnew = 0;

		fin.clear();
		fin.seekg(0, ios::beg);

		if (flnew == 1)
			LoadTable_NewFormat(fin);
		else
			LoadTable_OldFormat(fin);

		fin.close();

	}

	string decayprefix = "";

	for (int i = InputFile.size() - 1; i >= 0; --i)
		if (InputFile[i] == '\\' || InputFile[i] == '/')
		{
			decayprefix = InputFile.substr(0, i + 1);
			break;
		}

	LoadDecays((decayprefix + "decays.dat").c_str(), GenAntiP);

}

void ThermalParticleSystem::LoadTable_OldFormat(std::ifstream & fin, bool GenerateAntiParticles, double mcut)
{
	if (fin.is_open()) {
		string tmp;
		char tmpc[500];
		fin.getline(tmpc, 500);
		tmp = string(tmpc);
		while (1) {
			vector<string> fields = CuteHRGHelper::split(tmp, ';');
			if (fields.size() < 14) break;
			int stable, pdgid, spin, stat, str, bary, chg, charm;
			double mass, width, threshold, abss, absc, radius = 0.5;
			string name, decayname = "";
			stable = atoi(fields[0].c_str());
			name = fields[1];
			pdgid = atoi(fields[2].c_str());
			spin = atoi(fields[3].c_str());
			stat = atoi(fields[4].c_str());
			mass = atof(fields[5].c_str());
			str = atoi(fields[6].c_str());
			bary = atoi(fields[7].c_str());
			chg = atoi(fields[8].c_str());
			charm = atoi(fields[9].c_str());
			abss = atof(fields[10].c_str());
			absc = atof(fields[11].c_str());
			width = atof(fields[12].c_str());
			threshold = atof(fields[13].c_str());
			if (fields.size() >= 15) radius = atof(fields[14].c_str());
			if (fields.size() == 16) decayname = fields[15];

			if (mass > mcut) {
				fin.getline(tmpc, 500);
				tmp = string(tmpc);
				continue;
			}

			if (bary != 0)  m_NumBaryons++;
			if (chg != 0)   m_NumCharged++;
			if (str != 0)   m_NumStrange++;
			if (charm != 0) m_NumCharmed++;

			m_Particles.push_back(ThermalParticle((bool)stable, name, pdgid, spin, stat, mass, str, bary, chg, abss, width, threshold, charm, absc, radius));
			m_NumberOfParticles++;
			m_IDtoPDG.push_back(pdgid);
			m_PDGtoID[pdgid] = m_IDtoPDG.size() - 1;

			if (GenerateAntiParticles && !(bary == 0 && chg == 0 && str == 0 && charm == 0)) {

				if (bary != 0)  m_NumBaryons++;
				if (chg != 0)   m_NumCharged++;
				if (str != 0)   m_NumStrange++;
				if (charm != 0) m_NumCharmed++;

				if (bary == 0 && name[name.size() - 1] == '+')
					name[name.size() - 1] = '-';
				else if (bary == 0 && name[name.size() - 1] == '-')
					name[name.size() - 1] = '+';
				else
					name = "anti-" + name; 
				m_Particles.push_back(ThermalParticle((bool)stable, name, -pdgid, spin, stat, mass, -str, -bary, -chg, abss, width, threshold, -charm, absc, radius));
				m_Particles[m_Particles.size() - 1].SetAntiParticle(true);
				m_IDtoPDG.push_back(-pdgid);
				m_PDGtoID[-pdgid] = m_IDtoPDG.size() - 1;
			}

			fin.getline(tmpc, 500);
			tmp = string(tmpc);
		}

		sort(m_Particles.begin(), m_Particles.end(), ::cmpParticleMass);
		m_PDGtoID.clear();
		m_IDtoPDG.resize(0);
		for (int i = 0; i < m_Particles.size(); ++i) {
			m_IDtoPDG.push_back(m_Particles[i].PdgId());
			m_PDGtoID[m_Particles[i].PdgId()] = m_IDtoPDG.size() - 1;
		}
	}
}

void ThermalParticleSystem::LoadTable_NewFormat(std::ifstream & fin, bool GenerateAntiParticles, double mcut)
{
	if (fin.is_open()) {
		char cc[2000];
		while (!fin.eof()) {
			fin.getline(cc, 2000);
			string tmp = string(cc);
			vector<string> elems = CuteHRGHelper::split(tmp, '#');
			if (elems.size() < 1)
				continue;

			istringstream iss(elems[0]);

			int stable, pdgid, stat, str, bary, chg, charm;
			double mass, degeneracy, width, threshold, abss, absc;
			string name;

			if (iss >> pdgid 
				>> name
				>> stable
				>> mass
				>> degeneracy
				>> stat
				>> bary
				>> chg
				>> str
				>> charm
				>> abss
				>> absc
				>> width
				>> threshold) {

				if (mass > mcut) 
					continue;

				if (bary != 0)  m_NumBaryons++;
				if (chg != 0)   m_NumCharged++;
				if (str != 0)   m_NumStrange++;
				if (charm != 0) m_NumCharmed++;

				m_Particles.push_back(ThermalParticle((bool)stable, name, pdgid, degeneracy, stat, mass, str, bary, chg, abss, width, threshold, charm, absc));
				m_NumberOfParticles++;
				m_IDtoPDG.push_back(pdgid);
				m_PDGtoID[pdgid] = m_IDtoPDG.size() - 1;

				if (GenerateAntiParticles && !(bary == 0 && chg == 0 && str == 0 && charm == 0)) {

					if (bary != 0)  m_NumBaryons++;
					if (chg != 0)   m_NumCharged++;
					if (str != 0)   m_NumStrange++;
					if (charm != 0) m_NumCharmed++;

					if (bary == 0 && name[name.size() - 1] == '+')
						name[name.size() - 1] = '-';
					else if (bary == 0 && name[name.size() - 1] == '-')
						name[name.size() - 1] = '+';
					else
						name = "anti-" + name; 
					m_Particles.push_back(ThermalParticle((bool)stable, name, -pdgid, degeneracy, stat, mass, -str, -bary, -chg, abss, width, threshold, -charm, absc));
					m_Particles[m_Particles.size() - 1].SetAntiParticle(true);
					m_IDtoPDG.push_back(-pdgid);
					m_PDGtoID[-pdgid] = m_IDtoPDG.size() - 1;
				}
			}
		}

		sort(m_Particles.begin(), m_Particles.end(), ::cmpParticleMass);
		m_PDGtoID.clear();
		m_IDtoPDG.resize(0);
		for (int i = 0; i < m_Particles.size(); ++i) {
			m_IDtoPDG.push_back(m_Particles[i].PdgId());
			m_PDGtoID[m_Particles[i].PdgId()] = m_IDtoPDG.size() - 1;
		}
	}
}

void ThermalParticleSystem::WriteTableToFile(std::string OutputFile, bool WriteAntiParticles)
{
	std::ofstream fout(OutputFile.c_str());
	if (fout.is_open()) {
		fout << "#"
			<< std::setw(14) << "pdgid"
			<< std::setw(15) << "name"
			<< std::setw(15) << "stable"
			<< std::setw(15) << "mass[GeV]"
			<< std::setw(15) << "degeneracy"
			<< std::setw(15) << "statistics"
			<< std::setw(15) << "B"
			<< std::setw(15) << "Q"
			<< std::setw(15) << "S"
			<< std::setw(15) << "C"
			<< std::setw(15) << "|S|"
			<< std::setw(15) << "|C|"
			<< std::setw(15) << "width[GeV]"
			<< std::setw(15) << "threshold[GeV]"
			<< std::endl;

		for (int i = 0; i < m_Particles.size(); ++i) {
			const ThermalParticle& part = m_Particles[i];
			if (part.PdgId() < 0 && !WriteAntiParticles)
				continue;

			fout << std::setw(15) << part.PdgId()
				<< std::setw(15) << part.Name()
				<< std::setw(15) << static_cast<int>(part.IsStable())
				<< std::setw(15) << part.Mass()
				<< std::setw(15) << part.Degeneracy()
				<< std::setw(15) << part.Statistics()
				<< std::setw(15) << part.BaryonCharge()
				<< std::setw(15) << part.ElectricCharge()
				<< std::setw(15) << part.Strangeness()
				<< std::setw(15) << part.Charm()
				<< std::setw(15) << part.AbsoluteStrangeness()
				<< std::setw(15) << part.AbsoluteCharm()
				<< std::setw(15) << part.ResonanceWidth()
				<< std::setw(15) << part.DecayThresholdMass()
				<< std::endl;
		}
		fout.close();
	}
}

void ThermalParticleSystem::LoadDecays(std::string DecaysFile, bool GenerateAntiParticles)
{
	for (int i = 0; i < m_Particles.size(); ++i)
		m_Particles[i].ClearDecays();

	vector< vector<ParticleDecay> > decays(0);
	vector<int> pdgids(0);
	map<int, int> decaymap;
	decaymap.clear();
	
	ifstream fin(DecaysFile.c_str());
	if (fin.is_open()) {
		int decaypartnumber = 0;
		fin >> decaypartnumber;
		decays.reserve(decaypartnumber);

		for (unsigned int i = 0; i < decaypartnumber; ++i) {
			int pdgid, decaysnumber, tmpid, daughters;
			double bratio;
			fin >> pdgid >> decaysnumber;
			decaymap[pdgid] = i;
			decays.push_back(vector<ParticleDecay>(0));
			pdgids.push_back(pdgid);
			for (unsigned int j = 0; j < decaysnumber; ++j) {
				ParticleDecay decay;
				fin >> bratio;
				decay.mBratio = bratio / 100.;
				fin >> daughters;
				decay.mDaughters.reserve(daughters);
				for (unsigned int k = 0; k < daughters; ++k) {
					fin >> tmpid;
					decay.mDaughters.push_back(tmpid);
				}
				decays[i].push_back(decay);
			}
		}
		fin.close();
	}

	for (int i = 0; i < m_Particles.size(); ++i) {
		if (decaymap.count(m_Particles[i].PdgId()) != 0)
			m_Particles[i].SetDecays(decays[decaymap[m_Particles[i].PdgId()]]);
	}

	if (GenerateAntiParticles) {
		for (int i = 0; i < m_Particles.size(); ++i) {
			if (m_Particles[i].PdgId() < 0)
				m_Particles[i].SetDecays(GetDecaysFromAntiParticle(m_Particles[m_PDGtoID[-m_Particles[i].PdgId()]].Decays()));
		}
	}

	for (int i = 0; i < m_Particles.size(); ++i)
		m_Particles[i].SetDecaysOriginal(m_Particles[i].Decays());
	fin.close();

	FillDecayProperties();
	FillDecayThresholds();
	ProcessDecays();
}

std::string ThermalParticleSystem::GetNameFromPDG(int pdgid) {
	if (pdgid == 1) return string("Npart");
	if (pdgid == 33340) return string("Omega+Omegabar");
	if (pdgid == 22120) return string("p+n");
	if (m_PDGtoID.count(pdgid) == 0) return string("???");
	else return m_Particles[m_PDGtoID[pdgid]].Name();
}

void ThermalParticleSystem::NormalizeBranchingRatios() {
	for (int i = 0; i < m_Particles.size(); ++i) m_Particles[i].NormalizeBranchingRatios();
	ProcessDecays();
}


void ThermalParticleSystem::RestoreBranchingRatios() {
	for (int i = 0; i < m_Particles.size(); ++i) m_Particles[i].RestoreBranchingRatios();
	ProcessDecays();
}

void ThermalParticleSystem::SetCalculationType(IdealGasFunctions::QStatsCalculationType type)
{
	for (int i = 0; i < m_Particles.size(); ++i)
		m_Particles[i].SetCalculationType(type);
}

void ThermalParticleSystem::SetClusterExpansionOrder(int order)
{
	for (int i = 0; i < m_Particles.size(); ++i)
		m_Particles[i].SetClusterExpansionOrder(order);
}

void ThermalParticleSystem::SetResonanceWidthShape(ThermalParticle::ResonanceWidthShape shape)
{
	for (int i = 0; i < m_Particles.size(); ++i)
		m_Particles[i].SetResonanceWidthShape(shape);
}

void ThermalParticleSystem::SetResonanceWidthIntegrationType(ThermalParticle::ResonanceWidthIntegration type)
{
	bool dodecays = (type != m_ResonanceWidthIntegrationType);

	m_ResonanceWidthIntegrationType = type; 
	
	for (int i = 0; i < m_Particles.size(); ++i)
		m_Particles[i].SetResonanceWidthIntegrationType(type);

	if (dodecays)
		ProcessDecays();
}

void ThermalParticleSystem::SeparateDecaysIntoWeakAndStrong() {
	set<int> weakPDG;
	weakPDG.insert(310);
	weakPDG.insert(3122); weakPDG.insert(-3122);
	weakPDG.insert(3222); weakPDG.insert(-3222);
	weakPDG.insert(3112); weakPDG.insert(-3112);
	weakPDG.insert(3322); weakPDG.insert(-3322);
	weakPDG.insert(3312); weakPDG.insert(-3312);
	weakPDG.insert(3334); weakPDG.insert(-3334);

	for (int i = 0; i < m_Particles.size(); ++i) {
		if (weakPDG.count(m_Particles[i].PdgId()) > 0) {
			m_Particles[i].SetDecayType(2);
		}
		else if (!m_Particles[i].IsStable()) m_Particles[i].SetDecayType(1);
		else m_Particles[i].SetDecayType(0);
	}
}

const ThermalParticle & ThermalParticleSystem::Particle(int id) const
{
	if (id < 0 || id >= m_Particles.size()) {
		printf("**ERROR** ThermalParticleSystem::Particle(int id): id is out of bounds!");
		exit(1);
	}
	return m_Particles[id];
}

ThermalParticle & ThermalParticleSystem::Particle(int id)
{
	if (id < 0 || id >= m_Particles.size()) {
		printf("**ERROR** ThermalParticleSystem::Particle(int id): id is out of bounds!");
		exit(1);
	}
	return m_Particles[id];
}

ThermalParticle & ThermalParticleSystem::ParticleByPDG(int pdgid)
{
	if (m_PDGtoID.count(pdgid)==0) {
		printf("**ERROR** ThermalParticleSystem::ParticleByPDG(int pdgid): pdgid is unknown!");
		exit(1);
	}
	return m_Particles[m_PDGtoID[pdgid]];
}

void ThermalParticleSystem::FillResonanceWeakDecays() {
	for (int i = 0; i < m_Particles.size(); ++i) {
		m_Particles[i].WeakDecayContributions().resize(0);
	}
	for (int i = m_Particles.size() - 1; i >= 0; i--)
		if (m_Particles[i].DecayType() != 0) {
			GoResonanceWeak(i, i, 1.);
		}
}

void ThermalParticleSystem::GoResonanceWeak(int ind, int startind, double BR) {
	if (ind != startind && m_Particles[ind].WeakDecayContributions().size() > 0 && m_Particles[ind].WeakDecayContributions()[m_Particles[ind].WeakDecayContributions().size() - 1].second == startind)
	{
		m_Particles[ind].WeakDecayContributions()[m_Particles[ind].WeakDecayContributions().size() - 1].first += BR;
	}
	else if (ind != startind) m_Particles[ind].WeakDecayContributions().push_back(make_pair(BR, startind));
	if (m_Particles[ind].DecayType() != 0) {
		for (int i = 0; i < m_Particles[ind].Decays().size(); ++i) {
			double tbr = m_Particles[ind].Decays()[i].mBratio;

			// TODO: Fix(?) for canonical ensemble
			if (m_ResonanceWidthIntegrationType == ThermalParticle::eBW && ind == startind)
				tbr = m_Particles[ind].Decays()[i].mBratioAverage;

			for (int j = 0; j < m_Particles[ind].Decays()[i].mDaughters.size(); ++j) {
				if (m_PDGtoID.count(m_Particles[ind].Decays()[i].mDaughters[j]) != 0) 
					GoResonanceWeak(m_PDGtoID[m_Particles[ind].Decays()[i].mDaughters[j]], startind, BR*tbr);
			}
		}
	}
}
