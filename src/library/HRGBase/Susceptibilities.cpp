/*
 * Thermal-FIST package
 *
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "HRGBase/Susceptibilities.h"

#include <iostream>
#include <sstream>
#include <cassert>

#include "HRGBase/xMath.h"


// For time keeping
// Windows
#ifdef _WIN32
#include <Windows.h>
#else
#include <time.h>
#include <sys/time.h>
#endif

using namespace std;

namespace thermalfist {

  std::map<std::vector<int>, double> ComputeBQSSusceptibilities(ThermalModelBase *model, int order, double dmuTnum) {
    assert(order >= 0 && order <= 4);

    std::map<std::vector<int>, double> ret;

    double T = model->Parameters().T;
    double muB = model->Parameters().muB;
    double muQ = model->Parameters().muQ;
    double muS = model->Parameters().muS;

    model->CalculatePrimordialDensities();

    // First the scaled pressure p/T^4
    ret[std::vector<int>{0, 0, 0}] = model->Pressure() / xMath::GeVtoifm3() / pow(T, 4);

    if (order == 0)
      return ret;

    // Next the scaled densities n/T^3
    ret[std::vector<int>{1, 0, 0}] = model->BaryonDensity() / xMath::GeVtoifm3() / pow(T, 3);
    ret[std::vector<int>{0, 1, 0}] = model->ElectricChargeDensity() / xMath::GeVtoifm3() / pow(T, 3);
    ret[std::vector<int>{0, 0, 1}] = model->StrangenessDensity() / xMath::GeVtoifm3() / pow(T, 3);

    if (order == 1)
      return ret;

    model->CalculateFluctuations();

    // Next the scaled susceptibilities
    ret[std::vector<int>{2, 0, 0}] = model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge);
    ret[std::vector<int>{1, 1, 0}] = model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::ElectricCharge);
    ret[std::vector<int>{1, 0, 1}] = model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge);
    ret[std::vector<int>{0, 2, 0}] = model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
    ret[std::vector<int>{0, 1, 1}] = model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge);
    ret[std::vector<int>{0, 0, 2}] = model->Susc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge);

    if (order == 2)
      return ret;

    // Third order and above, need to use central finite difference method

    // dmuB
    model->SetBaryonChemicalPotential(muB + dmuTnum * T);
    auto SuscmuBp = ComputeBQSSusceptibilities(model, 2);
    model->SetBaryonChemicalPotential(muB - dmuTnum * T);
    auto SuscmuBm = ComputeBQSSusceptibilities(model, 2);
    model->SetBaryonChemicalPotential(muB);

    // dmuQ
    model->SetElectricChemicalPotential(muQ + dmuTnum * T);
    auto SuscmuQp = ComputeBQSSusceptibilities(model, 2);
    model->SetElectricChemicalPotential(muQ - dmuTnum * T);
    auto SuscmuQm = ComputeBQSSusceptibilities(model, 2);
    model->SetElectricChemicalPotential(muQ);

    // dmuS
    model->SetStrangenessChemicalPotential(muS + dmuTnum * T);
    auto SuscmuSp = ComputeBQSSusceptibilities(model, 2);
    model->SetStrangenessChemicalPotential(muS - dmuTnum * T);
    auto SuscmuSm = ComputeBQSSusceptibilities(model, 2);
    model->SetStrangenessChemicalPotential(muS);

    // Iterate over all the keys and calculate the third order susceptibilities
    for(int imuB = 0; imuB <= 3; imuB++) {
      for(int imuQ = 0; imuQ <= 3; imuQ++) {
        for(int imuS = 0; imuS <= 3; imuS++) {
          if(imuB + imuQ + imuS == 3) {
            std::vector<int> key{imuB, imuQ, imuS};
            if (imuB > 0) {
              std::vector<int> key2 = std::vector<int>{imuB - 1, imuQ, imuS};
              ret[key] = (SuscmuBp[key2] - SuscmuBm[key2]) / (2 * dmuTnum);
            } else if (imuS > 0) {
              std::vector<int> key2 = std::vector<int>{imuB, imuQ, imuS - 1};
              ret[key] = (SuscmuSp[key2] - SuscmuSm[key2]) / (2 * dmuTnum);
            } else {
              assert(imuQ > 0);
              std::vector<int> key2 = std::vector<int>{imuB, imuQ - 1, imuS};
              ret[key] = (SuscmuQp[key2] - SuscmuQm[key2]) / (2 * dmuTnum);
            }
          }
        }
      }
    }

    if (order == 3)
      return ret;

    // Iterate over all the keys and calculate the fourth order susceptibilities
    for(int imuB = 0; imuB <= 4; imuB++) {
      for(int imuQ = 0; imuQ <= 4; imuQ++) {
        for(int imuS = 0; imuS <= 4; imuS++) {
          if(imuB + imuQ + imuS == 4) {
            std::vector<int> key{imuB, imuQ, imuS};
            if (imuB > 1) {
              std::vector<int> key2 = std::vector<int>{imuB - 2, imuQ, imuS};
              ret[key] = (SuscmuBp[key2] - 2 * ret[key2] + SuscmuBm[key2]) / (dmuTnum * dmuTnum);
            } else if (imuS > 1) {
              std::vector<int> key2 = std::vector<int>{imuB, imuQ, imuS - 2};
              ret[key] = (SuscmuSp[key2] - 2 * ret[key2] + SuscmuSm[key2]) / (dmuTnum * dmuTnum);
            } else {
              assert(imuQ > 1);
              std::vector<int> key2 = std::vector<int>{imuB, imuQ - 2, imuS};
              ret[key] = (SuscmuQp[key2] - 2 * ret[key2] + SuscmuQm[key2]) / (dmuTnum * dmuTnum);
            }
          }
        }
      }
    }

    return ret;
  }

  std::map<std::vector<int>, double> ComputeBQSSusceptibilitiesDerivativeT(ThermalModelBase *model, int order, double dmuTnum) {
    assert(order >= 0 && order <= 4);

    std::map<std::vector<int>, double> ret;

    double T = model->Parameters().T;
    double muB = model->Parameters().muB;
    double muQ = model->Parameters().muQ;
    double muS = model->Parameters().muS;

    model->CalculatePrimordialDensities();

    // First the T derivative of the scaled pressure p/T^4
    // The derivative of the scaled pressure p/T^4 with respect to temperature T
    // is given by the expression: d(p/T^4)/dT = (s/T^4) - (4*p/T^5)
    // where s is the entropy density and p is the pressure.
    ret[std::vector<int>{0, 0, 0}] = (model->EntropyDensity() - 4. * model->Pressure() / T) / xMath::GeVtoifm3() / pow(T, 4);

    if (order == 0)
      return ret;

    model->CalculateTemperatureDerivatives();

    // Next the T-derivatives of the scaled densities n/T^3
    // d(n/T^3)/dT = (dn/dT)/T^3 - 3*n/T^4
    ret[std::vector<int>{1, 0, 0}] = (model->ConservedChargeDensitydT(ConservedCharge::BaryonCharge) - 3. * model->BaryonDensity() / T) / xMath::GeVtoifm3() / pow(T, 3);
    ret[std::vector<int>{0, 1, 0}] = (model->ConservedChargeDensitydT(ConservedCharge::ElectricCharge) - 3. * model->ElectricChargeDensity() / T) / xMath::GeVtoifm3() / pow(T, 3);
    ret[std::vector<int>{0, 0, 1}] = (model->ConservedChargeDensitydT(ConservedCharge::StrangenessCharge) - 3. * model->StrangenessDensity() / T) / xMath::GeVtoifm3() / pow(T, 3);

    if (order == 1)
      return ret;

    model->CalculateFluctuations();

    // Next the scaled susceptibilities
    ret[std::vector<int>{2, 0, 0}] = model->dSuscdT(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge);
    ret[std::vector<int>{1, 1, 0}] = model->dSuscdT(ConservedCharge::BaryonCharge, ConservedCharge::ElectricCharge);
    ret[std::vector<int>{1, 0, 1}] = model->dSuscdT(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge);
    ret[std::vector<int>{0, 2, 0}] = model->dSuscdT(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge);
    ret[std::vector<int>{0, 1, 1}] = model->dSuscdT(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge);
    ret[std::vector<int>{0, 0, 2}] = model->dSuscdT(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge);

    if (order == 2)
      return ret;

    // Third order and above, need to use central finite difference method

    auto Susc = ComputeBQSSusceptibilities(model, order);

    // dmuB
    model->SetBaryonChemicalPotential(muB + dmuTnum * T);
    auto dSuscmuBp = ComputeBQSSusceptibilitiesDerivativeT(model, 2);
    model->SetBaryonChemicalPotential(muB - dmuTnum * T);
    auto dSuscmuBm = ComputeBQSSusceptibilitiesDerivativeT(model, 2);
    model->SetBaryonChemicalPotential(muB);

    // dmuQ
    model->SetElectricChemicalPotential(muQ + dmuTnum * T);
    auto dSuscmuQp = ComputeBQSSusceptibilitiesDerivativeT(model, 2);
    model->SetElectricChemicalPotential(muQ - dmuTnum * T);
    auto dSuscmuQm = ComputeBQSSusceptibilitiesDerivativeT(model, 2);
    model->SetElectricChemicalPotential(muQ);

    // dmuS
    model->SetStrangenessChemicalPotential(muS + dmuTnum * T);
    auto dSuscmuSp = ComputeBQSSusceptibilitiesDerivativeT(model, 2);
    model->SetStrangenessChemicalPotential(muS - dmuTnum * T);
    auto dSuscmuSm = ComputeBQSSusceptibilitiesDerivativeT(model, 2);
    model->SetStrangenessChemicalPotential(muS);

    // Iterate over all the keys and calculate the third order susceptibilities
    for(int imuB = 0; imuB <= 3; imuB++) {
      for(int imuQ = 0; imuQ <= 3; imuQ++) {
        for(int imuS = 0; imuS <= 3; imuS++) {
          if(imuB + imuQ + imuS == 3) {
            std::vector<int> key{imuB, imuQ, imuS};
            ret[key] = Susc[key] / T;
            if (imuB > 0) {
              std::vector<int> key2 = std::vector<int>{imuB - 1, imuQ, imuS};
              ret[key] += (dSuscmuBp[key2] - dSuscmuBm[key2]) / (2 * dmuTnum);
            } else if (imuS > 0) {
              std::vector<int> key2 = std::vector<int>{imuB, imuQ, imuS - 1};
              ret[key] += (dSuscmuSp[key2] - dSuscmuSm[key2]) / (2 * dmuTnum);
            } else {
              assert(imuQ > 0);
              std::vector<int> key2 = std::vector<int>{imuB, imuQ - 1, imuS};
              ret[key] += (dSuscmuQp[key2] - dSuscmuQm[key2]) / (2 * dmuTnum);
            }
          }
        }
      }
    }

    if (order == 3)
      return ret;

    // Iterate over all the keys and calculate the fourth order susceptibilities
    for(int imuB = 0; imuB <= 4; imuB++) {
      for(int imuQ = 0; imuQ <= 4; imuQ++) {
        for(int imuS = 0; imuS <= 4; imuS++) {
          if(imuB + imuQ + imuS == 4) {
            std::vector<int> key{imuB, imuQ, imuS};
            ret[key] = 2. * Susc[key] / T;
            if (imuB > 1) {
              std::vector<int> key2 = std::vector<int>{imuB - 2, imuQ, imuS};
              ret[key] += (dSuscmuBp[key2] - 2 * ret[key2] + dSuscmuBm[key2]) / (dmuTnum * dmuTnum);
            } else if (imuS > 1) {
              std::vector<int> key2 = std::vector<int>{imuB, imuQ, imuS - 2};
              ret[key] += (dSuscmuSp[key2] - 2 * ret[key2] + dSuscmuSm[key2]) / (dmuTnum * dmuTnum);
            } else {
              assert(imuQ > 1);
              std::vector<int> key2 = std::vector<int>{imuB, imuQ - 2, imuS};
              ret[key] += (dSuscmuQp[key2] - 2 * ret[key2] + dSuscmuQm[key2]) / (dmuTnum * dmuTnum);
            }
          }
        }
      }
    }

    return ret;
  }

} // namespace thermalfist
