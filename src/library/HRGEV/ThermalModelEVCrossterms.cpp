/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2015-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEV/ThermalModelEVCrossterms.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>

#include "HRGBase/xMath.h"
#include "HRGEV/ExcludedVolumeHelper.h"

#include <Eigen/Dense>

using namespace std;

namespace thermalfist {

  void ThermalModelEVCrossterms::ReadInteractionParameters(const std::string & filename)
  {
    m_Virial = std::vector< std::vector<double> >(m_TPS->Particles().size(), std::vector<double>(m_TPS->Particles().size(), 0.));

    ifstream fin(filename.c_str());
    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1)
        continue;
      istringstream iss(elems[0]);
      long long pdgid1, pdgid2;
      double b;
      if (iss >> pdgid1 >> pdgid2 >> b) {
        int ind1 = m_TPS->PdgToId(pdgid1);
        int ind2 = m_TPS->PdgToId(pdgid2);
        if (ind1 != -1 && ind2 != -1)
          m_Virial[ind1][ind2] = b;
      }
    }
    fin.close();
  }

  void ThermalModelEVCrossterms::WriteInteractionParameters(const std::string & filename)
  {
    ofstream fout(filename.c_str());
    fout << "# List of crossterms parameters to be used in the Crossterms excluded-volume HRG model"
      << std::endl;
    fout << "# Only particle pairs with a non-zero eigenvolume parameter are listed here"
      << std::endl;
    /*fout << "#" << std::setw(14) << "pdg_i"
      << std::setw(15) << "pdg_j"
      << std::setw(15) << "b_{ij}[fm^3]"
      << std::endl;*/
    fout << "#" << " " << "pdg_i"
      << " " << "pdg_j"
      << " " << "b_{ij}[fm^3]"
      << std::endl;
    for (int i = 0; i < m_TPS->ComponentsNumber(); ++i) {
      for (int j = 0; j < m_TPS->ComponentsNumber(); ++j) {
        if (m_Virial[i][j] != 0.) {
          fout << " " << m_TPS->Particle(i).PdgId();
          fout << " " << m_TPS->Particle(j).PdgId();
          fout << " " << m_Virial[i][j];
          fout << std::endl;
        }
      }
    }
    fout.close();
  }

  void ThermalModelEVCrossterms::SetAttraction(int i, int j, double a)
  {
    if (a != 0.0) {
      printf("**WARNING** ThermalModelEVCrossterms::SetAttraction(): Trying to include attraction into ThermalModelEVCrossterms()!\n");
    }
    //m_Attr[i][j] = 0.0;
  }

  void ThermalModelEVCrossterms::SetMultipleSolutionsMode(bool search)
  {
    if (search) {
      printf("**WARNING** ThermalModelEVCrossterms::SetMultipleSolutionsMode(): Trying to search for multiple solutions in ThermalModelEVCrossterms()! There is no need.\n");
    }
    m_SearchMultipleSolutions = false;
  }

} // namespace thermalfist