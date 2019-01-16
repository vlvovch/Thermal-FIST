/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2016-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGEV/ExcludedVolumeHelper.h"

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
