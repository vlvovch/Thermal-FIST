/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2016-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEV/ExcludedVolumeHelper.h"


double thermalfist::CuteHRGHelper::btilrr(double r1, double r2)
{
  double bij = brr(r1, r2);
  if (bij == 0.0) return 0.;
  double bii = brr(r1, r1);
  if (bii == 0.0) return 0.;
  double bjj = brr(r2, r2);
  return 2. * bij * bii / (bii + bjj);
}

std::vector< std::vector<double> > thermalfist::CuteHRGHelper::bijMatrix(const ThermalModelBase * model)
{
  std::vector< std::vector<double> > ret(model->ComponentsNumber());
  for (int i = 0; i < model->ComponentsNumber(); ++i) {
    ret[i].resize(model->ComponentsNumber());
    for (int j = 0; j < model->ComponentsNumber(); ++j) {
      ret[i][j] = model->RepulsionCoefficient(i, j);
    }
  }
  return ret;
}

std::vector< std::vector<double> > thermalfist::CuteHRGHelper::aijMatrix(const ThermalModelBase * model)
{
  std::vector< std::vector<double> > ret(model->ComponentsNumber());
  for (int i = 0; i < model->ComponentsNumber(); ++i) {
    ret[i].resize(model->ComponentsNumber());
    for (int j = 0; j < model->ComponentsNumber(); ++j) {
      ret[i][j] = model->AttractionCoefficient(i, j);
    }
  }
  return ret;
}

std::vector<std::vector<double>> thermalfist::GetBaryonBaryonInteractionMatrix(const ThermalParticleSystem *TPS, double param) {
  std::vector<std::vector<double>> ret(TPS->Particles().size(), std::vector<double>(TPS->Particles().size(), 0.));
  for (int i1 = 0; i1 < TPS->Particles().size(); ++i1) {
    for (int i2 = 0; i2 < TPS->Particles().size(); ++i2) {
      // Baryon charge of first species
      int B1 = TPS->Particles()[i1].BaryonCharge();
      // Baryon charge of second species
      int B2 = TPS->Particles()[i2].BaryonCharge();
      if ((B1 > 0 && B2 > 0) || (B1 < 0 && B2 < 0)) {
        ret[i1][i2] = param;
      } else {
        ret[i1][i2] = 0.;
      }
    }
  }
  return ret;
}