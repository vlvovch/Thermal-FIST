/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "resultdialog.h"

#include <QLayout>
#include <QHeaderView>
#include <QTimer>
#include <QLabel>
#include <QApplication>
#include <QDebug>
#include <iostream>

#include "HRGBase/ThermalModelBase.h"
#include "HRGBase/xMath.h"

using namespace thermalfist;

ResultDialog::ResultDialog(QWidget *parent, ThermalModelBase *mod, ChargesFluctuations *flucts_in) :
  QDialog(parent), model(mod), flucts(flucts_in)
{
  QVBoxLayout *layout = new QVBoxLayout();

  QFont font = QApplication::font();
  font.setPointSize(font.pointSize() + 2);
  QLabel *labParam = new QLabel(tr("Parameters:"));
  labParam->setFont(font);

  QFont font2("Monospace");
  font2.setStyleHint(QFont::TypeWriter);
  parameters = new QTextEdit();
  parameters->setReadOnly(true);
  parameters->setFont(font2);
  //data->setAcceptRichText(true);
  parameters->setPlainText(GetParameters());
  parameters->setMinimumHeight(350);
  parameters->setMinimumWidth(400);


  layout->addWidget(labParam);
  layout->addWidget(parameters);

  QLabel *labResults = new QLabel(tr("Calculation results:"));
  labResults->setFont(font);

  results = new QTextEdit();
  results->setReadOnly(true);
  results->setFont(font2);
  //data->setAcceptRichText(true);
  results->setPlainText(GetResults());
  results->setMinimumHeight(350);
  results->setMinimumWidth(550);

  QVBoxLayout *layoutResults = new QVBoxLayout();
  layoutResults->setAlignment(Qt::AlignTop);
  layoutResults->addWidget(labResults);
  layoutResults->addWidget(results);

  QHBoxLayout *mainLayout = new QHBoxLayout();

  mainLayout->addLayout(layout);
  mainLayout->addLayout(layoutResults);

  setLayout(mainLayout);

  this->setWindowTitle(tr("Equation of state properties"));

  QTimer::singleShot(0, this, SLOT(checkFixTableSize()));
}

QString ResultDialog::GetParameters() {
  QString ret = "";

  char cc[500];

  sprintf(cc, "%-20s = ", "T");
  //ret += "T\t\t= ";
  ret += QString(cc);
  ret += QString::number(model->Parameters().T*1.e3) + " MeV\r\n";

  if (model->Ensemble() != ThermalModelBase::CE) {
    sprintf(cc, "%-20s = ", "\\mu_B");
    ret += QString(cc);
    ret += QString::number(model->Parameters().muB*1.e3) + " MeV\r\n";

    if (model->Ensemble() != ThermalModelBase::SCE) {
      sprintf(cc, "%-20s = ", "\\mu_S");
      ret += QString(cc);
      ret += QString::number(model->Parameters().muS*1.e3) + " MeV\r\n";
    }

    sprintf(cc, "%-20s = ", "\\mu_Q");
    ret += QString(cc);
    ret += QString::number(model->Parameters().muQ*1.e3) + " MeV\r\n";

    if (model->TPS()->hasCharmed()) {
      sprintf(cc, "%-20s = ", "\\mu_C");
      ret += QString(cc);
      ret += QString::number(model->Parameters().muC*1.e3) + " MeV\r\n";
    }
  }
  else {
    sprintf(cc, "%-20s = ", "B");
    ret += QString(cc);
    ret += QString::number(model->Parameters().B) + "\r\n";

    sprintf(cc, "%-20s = ", "S");
    ret += QString(cc);
    ret += QString::number(model->Parameters().S) + "\r\n";

    sprintf(cc, "%-20s = ", "Q");
    ret += QString(cc);
    ret += QString::number(model->Parameters().Q) + "\r\n";

    if (model->TPS()->hasCharmed()) {
      sprintf(cc, "%-20s = ", "C");
      ret += QString(cc);
      ret += QString::number(model->Parameters().C) + "\r\n";
    }
  }

  sprintf(cc, "%-20s = ", "\\gamma_q");
  ret += QString(cc);
  ret += QString::number(model->Parameters().gammaq) + "\r\n";

  sprintf(cc, "%-20s = ", "\\gamma_S");
  ret += QString(cc);
  ret += QString::number(model->Parameters().gammaS) + "\r\n";

  if (model->TPS()->hasCharmed()) {
    sprintf(cc, "%-20s = ", "\\gamma_C");
    ret += QString(cc);
    ret += QString::number(model->Parameters().gammaC) + "\r\n";
  }

  sprintf(cc, "%-20s = ", "Volume");
  ret += QString(cc);
  ret += QString::number(model->Parameters().V) + " fm^3\r\n";

  sprintf(cc, "%-20s = ", "Finite widths");
  ret += QString(cc);
  if (model->UseWidth()) ret += "Yes\r\n";
  else ret += "No\r\n";

  sprintf(cc, "%-20s = ", "Quantum statistics");
  ret += QString(cc);
  if (model->QuantumStatistics()) ret += "Yes\r\n";
  else ret += "No\r\n";

  return ret;
}

QString ResultDialog::GetResults() {
  QString ret = "";

  char cc[500];

  sprintf(cc, "%-25s = ", "Total hadron density");
  ret += QString(cc);
  ret += QString::number(model->CalculateHadronDensity()) + " fm^-3\r\n";

  sprintf(cc, "%-25s = ", "Net baryon density");
  ret += QString(cc);
  ret += QString::number(model->CalculateBaryonDensity()) + " fm^-3\r\n";

  sprintf(cc, "%-25s = ", "Electric charge density");
  ret += QString(cc);
  ret += QString::number(model->CalculateChargeDensity()) + " fm^-3\r\n";

  sprintf(cc, "%-25s = ", "Net strangeness density");
  ret += QString(cc);
  ret += QString::number(model->CalculateStrangenessDensity()) + " fm^-3\r\n";

  sprintf(cc, "%-25s = ", "Net charm density");
  ret += QString(cc);
  ret += QString::number(model->CalculateCharmDensity()) + " fm^-3\r\n";

  ret += "\r\n";

  sprintf(cc, "%-25s = ", "Net baryon number");
  ret += QString(cc);
  ret += QString::number(model->CalculateBaryonDensity()*model->Volume()) + "\r\n";

  sprintf(cc, "%-25s = ", "Net electric charge");
  ret += QString(cc);
  ret += QString::number(model->CalculateChargeDensity()*model->Volume()) + "\r\n";

  sprintf(cc, "%-25s = ", "Net strangeness");
  ret += QString(cc);
  ret += QString::number(model->CalculateStrangenessDensity()*model->Volume()) + "\r\n";

  sprintf(cc, "%-25s = ", "Net charm");
  ret += QString(cc);
  ret += QString::number(model->CalculateCharmDensity()*model->Volume()) + "\r\n";

  ret += "\r\n";

  sprintf(cc, "%-25s = ", "Energy density");
  ret += QString(cc);
  ret += QString::number(model->CalculateEnergyDensity()*1.e3) + " MeV/fm^3\r\n";

  sprintf(cc, "%-25s = ", "Pressure");
  ret += QString(cc);
  ret += QString::number(model->CalculatePressure()*1.e3, 'g', 15) + " MeV/fm^3\r\n";
  //ret += QString::number(model->CalculatePressure() * 1.e3) + " MeV/fm^3\r\n";

  sprintf(cc, "%-25s = ", "Entropy density");
  ret += QString(cc);
  ret += QString::number(model->CalculateEntropyDensity()) + " fm^-3\r\n";

  //if (model->TAG()!="ThermalModel") {
  if (model->Ensemble() == ThermalModelBase::GCE && model->InteractionModel() == ThermalModelBase::Ideal) {
    sprintf(cc, "%-25s = ", "Baryon entropy fraction");
    ret += QString(cc);
    ret += QString::number(model->CalculateBaryonMatterEntropyDensity() / model->CalculateEntropyDensity()) + "\r\n";

    sprintf(cc, "%-25s = ", "Meson entropy fraction");
    ret += QString(cc);
    ret += QString::number(model->CalculateMesonMatterEntropyDensity() / model->CalculateEntropyDensity()) + "\r\n";
  }

  ret += "\r\n";

  sprintf(cc, "%-25s = ", "p/T^4");
  ret += QString(cc);
  ret += QString::number((model->CalculatePressure()) / model->Parameters().T / model->Parameters().T / model->Parameters().T / model->Parameters().T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm()) + "\r\n";

  sprintf(cc, "%-25s = ", "(e-3p)/T^4");
  ret += QString(cc);
  ret += QString::number((model->CalculateEnergyDensity() - 3.*model->CalculatePressure()) / model->Parameters().T / model->Parameters().T / model->Parameters().T / model->Parameters().T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm()) + "\r\n";

  sprintf(cc, "%-25s = ", "e/T^4");
  ret += QString(cc);
  ret += QString::number((model->CalculateEnergyDensity()) / model->Parameters().T / model->Parameters().T / model->Parameters().T / model->Parameters().T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm()) + "\r\n";

  sprintf(cc, "%-25s = ", "s/T^3");
  ret += QString(cc);
  ret += QString::number((model->CalculateEntropyDensity()) / model->Parameters().T / model->Parameters().T / model->Parameters().T / xMath::GeVtoifm() / xMath::GeVtoifm() / xMath::GeVtoifm()) + "\r\n";

  if (model->IsFluctuationsCalculated()) {
    ret += "\r\n";

    sprintf(cc, "%-25s = ", "\\chi2B");
    ret += QString(cc);
    ret += QString::number(model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi2Q");
    ret += QString(cc);
    ret += QString::number(model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi2S");
    ret += QString(cc);
    ret += QString::number(model->Susc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge)) + "\r\n";

    if (model->TPS()->hasCharmed()) {
      sprintf(cc, "%-25s = ", "\\chi2C");
      ret += QString(cc);
      ret += QString::number(model->Susc(ConservedCharge::CharmCharge, ConservedCharge::CharmCharge)) + "\r\n";
    }

    sprintf(cc, "%-25s = ", "\\chi11BQ");
    ret += QString(cc);
    ret += QString::number(model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::ElectricCharge)) + "\r\n";

    
    sprintf(cc, "%-25s = ", "\\chi11QS");
    ret += QString(cc);
    ret += QString::number(model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi11BS");
    ret += QString(cc);
    ret += QString::number(model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge)) + "\r\n";
      
    if (model->TPS()->hasStrange()) {
      sprintf(cc, "%-25s = ", "CBS");
      ret += QString(cc);
      ret += QString::number(-3. * model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge) / model->Susc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge)) + "\r\n";
    }

    sprintf(cc, "%-25s = ", "\\chi11BS/\\chi2S");
    ret += QString(cc);
    ret += QString::number(model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge) / model->Susc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi11QS/\\chi2S");
    ret += QString(cc);
    ret += QString::number(model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge) / model->Susc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge)) + "\r\n";
    
    sprintf(cc, "%-25s = ", "\\chi11QB/\\chi2B");
    ret += QString(cc);
    ret += QString::number(model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::BaryonCharge) / model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge)) + "\r\n";

    if (model->TPS()->hasCharmed()) {
      sprintf(cc, "%-25s = ", "\\chi11BC");
      ret += QString(cc);
      ret += QString::number(model->Susc(ConservedCharge::BaryonCharge, ConservedCharge::CharmCharge)) + "\r\n";

      sprintf(cc, "%-25s = ", "\\chi11QC");
      ret += QString(cc);
      ret += QString::number(model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::CharmCharge)) + "\r\n";

      sprintf(cc, "%-25s = ", "\\chi11SC");
      ret += QString(cc);
      ret += QString::number(model->Susc(ConservedCharge::StrangenessCharge, ConservedCharge::CharmCharge)) + "\r\n";
    }

    ret += "\r\n";

    sprintf(cc, "%-25s = ", "\\chi2prot");
    ret += QString(cc);
    ret += QString::number(model->ProxySusc(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi2Q");
    ret += QString(cc);
    ret += QString::number(model->ProxySusc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi2kaon");
    ret += QString(cc);
    ret += QString::number(model->ProxySusc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge)) + "\r\n";


    sprintf(cc, "%-25s = ", "\\chi11Q,p");
    ret += QString(cc);
    ret += QString::number(model->ProxySusc(ConservedCharge::BaryonCharge, ConservedCharge::ElectricCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi11Q,k");
    ret += QString(cc);
    ret += QString::number(model->ProxySusc(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi11p,k");
    ret += QString(cc);
    ret += QString::number(model->ProxySusc(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi11p,k/\\chi2k");
    ret += QString(cc);
    ret += QString::number(model->ProxySusc(ConservedCharge::BaryonCharge, ConservedCharge::StrangenessCharge) / model->ProxySusc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi11Q,k/\\chi2k");
    ret += QString(cc);
    ret += QString::number(model->ProxySusc(ConservedCharge::ElectricCharge, ConservedCharge::StrangenessCharge) / model->ProxySusc(ConservedCharge::StrangenessCharge, ConservedCharge::StrangenessCharge)) + "\r\n";

    sprintf(cc, "%-25s = ", "\\chi11Q,p/\\chi2p");
    ret += QString(cc);
    ret += QString::number(model->ProxySusc(ConservedCharge::ElectricCharge, ConservedCharge::BaryonCharge) / model->ProxySusc(ConservedCharge::BaryonCharge, ConservedCharge::BaryonCharge)) + "\r\n";


    ret += "\r\n";

    sprintf(cc, "%-25s = ", "Primordial Nch");
    ret += QString(cc);
    ret += QString::number(model->ChargedMultiplicity(0)) + "\r\n";

    sprintf(cc, "%-25s = ", "Primordial N+");
    ret += QString(cc);
    ret += QString::number(model->ChargedMultiplicity(1)) + "\r\n";

    sprintf(cc, "%-25s = ", "Primordial N-");
    ret += QString(cc);
    ret += QString::number(model->ChargedMultiplicity(-1)) + "\r\n";


    sprintf(cc, "%-25s = ", "Primordial w[Nch]");
    ret += QString(cc);
    ret += QString::number(model->ChargedScaledVariance(0)) + "\r\n";

    sprintf(cc, "%-25s = ", "Primordial w[N+]");
    ret += QString(cc);
    ret += QString::number(model->ChargedScaledVariance(1)) + "\r\n";

    sprintf(cc, "%-25s = ", "Primordial w[N-]");
    ret += QString(cc);
    ret += QString::number(model->ChargedScaledVariance(-1)) + "\r\n";

    ret += "\r\n";

    {
      double Np = model->ChargedMultiplicity(1);
      double Nm = model->ChargedMultiplicity(-1);
      double wplus = model->ChargedScaledVariance(1);
      double wminus = model->ChargedScaledVariance(-1);
      double covNpm =
              (Np * wplus + Nm * wminus - (model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge)
                                           * model->Volume() * pow(model->Parameters().T * xMath::GeVtoifm(), 3))) / 2.;

      double nudyn = wplus / Np + wminus / Nm - 2. * covNpm / Np / Nm - (Np + Nm) / Np / Nm;

      sprintf(cc, "%-25s = ", "Primordial <Nch>*nu_dyn");
      ret += QString(cc);
      ret += QString::number((Np + Nm) * nudyn) + "\r\n";

      double Dprim = 4. * model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge)
                     * model->Volume() * pow(model->Parameters().T * xMath::GeVtoifm(), 3) /
                     model->ChargedMultiplicity(0);

      sprintf(cc, "%-25s = ", "Primordial D");
      ret += QString(cc);
      ret += QString::number(Dprim) + "\r\n";

      ret += "\r\n";
    }


    sprintf(cc, "%-25s = ", "Final Nch");
    ret += QString(cc);
    ret += QString::number(model->ChargedMultiplicityFinal(0)) + "\r\n";

    sprintf(cc, "%-25s = ", "Final N+");
    ret += QString(cc);
    ret += QString::number(model->ChargedMultiplicityFinal(1)) + "\r\n";

    sprintf(cc, "%-25s = ", "Final N-");
    ret += QString(cc);
    ret += QString::number(model->ChargedMultiplicityFinal(-1)) + "\r\n";


    sprintf(cc, "%-25s = ", "Final w[Nch]");
    ret += QString(cc);
    ret += QString::number(model->ChargedScaledVarianceFinal(0)) + "\r\n";

    sprintf(cc, "%-25s = ", "Final w[N+]");
    ret += QString(cc);
    ret += QString::number(model->ChargedScaledVarianceFinal(1)) + "\r\n";

    sprintf(cc, "%-25s = ", "Final w[N-]");
    ret += QString(cc);
    ret += QString::number(model->ChargedScaledVarianceFinal(-1)) + "\r\n";

    ret += "\r\n";

    {
      double Np = model->ChargedMultiplicityFinal(1);
      double Nm = model->ChargedMultiplicityFinal(-1);
      double wplus = model->ChargedScaledVarianceFinal(1);
      double wminus = model->ChargedScaledVarianceFinal(-1);
      double covNpm =
              (Np * wplus + Nm * wminus - (model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge)
                                           * model->Volume() * pow(model->Parameters().T * xMath::GeVtoifm(), 3))) / 2.;

      double nudyn = wplus / Np + wminus / Nm - 2. * covNpm / Np / Nm - (Np + Nm) / Np / Nm;

      sprintf(cc, "%-25s = ", "Final <Nch>*nu_dyn");
      ret += QString(cc);
      ret += QString::number((Np + Nm) * nudyn) + "\r\n";

      double Dtot = 4. * model->Susc(ConservedCharge::ElectricCharge, ConservedCharge::ElectricCharge)
                     * model->Volume() * pow(model->Parameters().T * xMath::GeVtoifm(), 3) /
                     model->ChargedMultiplicityFinal(0);

      sprintf(cc, "%-25s = ", "Final D");
      ret += QString(cc);
      ret += QString::number(Dtot) + "\r\n";

      ret += "\r\n";
    }

    if (flucts != NULL && flucts->flag) {
      ret += "\r\n";

      sprintf(cc, "%-25s = ", "chi3B");
      ret += QString(cc);
      ret += QString::number(flucts->chi3B) + "\r\n";

      sprintf(cc, "%-25s = ", "chi4B");
      ret += QString(cc);
      ret += QString::number(flucts->chi4B) + "\r\n";

      sprintf(cc, "%-25s = ", "chi3B/chi2B");
      ret += QString(cc);
      ret += QString::number(flucts->chi3B / flucts->chi2B) + "\r\n";

      sprintf(cc, "%-25s = ", "chi4B/chi2B");
      ret += QString(cc);
      ret += QString::number(flucts->chi4B / flucts->chi2B) + "\r\n";


      ret += "\r\n";

      sprintf(cc, "%-25s = ", "chi3Q");
      ret += QString(cc);
      ret += QString::number(flucts->chi3Q) + "\r\n";

      sprintf(cc, "%-25s = ", "chi4Q");
      ret += QString(cc);
      ret += QString::number(flucts->chi4Q) + "\r\n";

      sprintf(cc, "%-25s = ", "chi3Q/chi2Q");
      ret += QString(cc);
      ret += QString::number(flucts->chi3Q / flucts->chi2Q) + "\r\n";

      sprintf(cc, "%-25s = ", "chi4Q/chi2Q");
      ret += QString(cc);
      ret += QString::number(flucts->chi4Q / flucts->chi2Q) + "\r\n";


      ret += "\r\n";

      if (model->TPS()->hasStrange()) {
        sprintf(cc, "%-25s = ", "chi3S");
        ret += QString(cc);
        ret += QString::number(flucts->chi3S) + "\r\n";

        sprintf(cc, "%-25s = ", "chi4S");
        ret += QString(cc);
        ret += QString::number(flucts->chi4S) + "\r\n";

        sprintf(cc, "%-25s = ", "chi3S/chi2S");
        ret += QString(cc);
        ret += QString::number(flucts->chi3S / flucts->chi2S) + "\r\n";

        sprintf(cc, "%-25s = ", "chi4S/chi2S");
        ret += QString(cc);
        ret += QString::number(flucts->chi4S / flucts->chi2S) + "\r\n";
      }


      if (model->TPS()->hasCharmed()) {
        ret += "\r\n";


        sprintf(cc, "%-25s = ", "chi3C");
        ret += QString(cc);
        ret += QString::number(flucts->chi3C) + "\r\n";

        sprintf(cc, "%-25s = ", "chi4C");
        ret += QString(cc);
        ret += QString::number(flucts->chi4C) + "\r\n";

        sprintf(cc, "%-25s = ", "chi3C/chi2C");
        ret += QString(cc);
        ret += QString::number(flucts->chi3C / flucts->chi2C) + "\r\n";

        sprintf(cc, "%-25s = ", "chi4C/chi2C");
        ret += QString(cc);
        ret += QString::number(flucts->chi4C / flucts->chi2C) + "\r\n";
      }
    }

  }

  return ret;
}

void ResultDialog::checkFixTableSize() {
}
