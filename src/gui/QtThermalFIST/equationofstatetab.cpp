/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "equationofstatetab.h"

#include <algorithm>
#include <fstream>
#include <sstream>

#include <QLayout>
#include <QGroupBox>
#include <QLabel>
#include <QTextStream>

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGEV/ThermalModelEVCrossterms.h"
#include "HRGVDW/ThermalModelVDW.h"
#include "HRGBase/xMath.h"
#include "HRGRealGas/ThermalModelRealGas.h"

using namespace std;
using namespace thermalfist;

void EoSWorker::run() {
  int i = 0;
  for (double param = Tmin; param <= Tmax + 1.e-10 && !(*stop); param += dT, ++i) {
    if (i < paramsTD->size()) {
      // T dependence at const. muB
      if (mode == 0) {
        model->SetTemperature(param * 1.e-3);
      }
      // T dependence at const. muB/T
      else if (mode == 1) {
        model->SetTemperature(param * 1.e-3);
        if (!model->ConstrainMuB())
          model->SetBaryonChemicalPotential(param * 1.e-3 * cParams[0]);
        if (!model->ConstrainMuQ())
          model->SetElectricChemicalPotential(param * 1.e-3 * cParams[1]);
        if (!model->ConstrainMuS())
          model->SetStrangenessChemicalPotential(param * 1.e-3 * cParams[2]);
      }
      // Isotherm
      else if (mode == 2) {
        model->SetBaryonChemicalPotential(param * 1.e-3);
        if (!model->ConstrainMuQ())
          model->SetElectricChemicalPotential(cParams[1] * 1.e-3);
        if (!model->ConstrainMuS())
          model->SetStrangenessChemicalPotential(cParams[2] * 1.e-3);
      }
      // Magnetic field
      else if (mode == 3) {
        model->SetTemperature(cParams[3] * 1.e-3);
        if (!model->ConstrainMuB())
          model->SetBaryonChemicalPotential(cParams[0] * 1.e-3);
        if (!model->ConstrainMuQ())
          model->SetElectricChemicalPotential(cParams[1] * 1.e-3);
        if (!model->ConstrainMuS())
          model->SetStrangenessChemicalPotential(cParams[2] * 1.e-3);
        model->SetMagneticField(param);
      }

      if ((mode == 0 || mode == 1) && config.ConstrainMuB && config.ConstrainMuBType == 1) {
        double totB = model->Volume() * config.RhoB;
        double totQ = config.QoverB * totB;
        double totS = 0.;
        double totC = 0.;
        double muBinit = model->Parameters().muB;
        double muQinit = model->Parameters().muQ;
        double muSinit = model->Parameters().muS;
        double muCinit = model->Parameters().muC;
        if (param == Tmin)
          muBinit = muQinit = muSinit = muCinit = 0.;
        model->SolveChemicalPotentials(
                totB, totQ, totS, totC,
                muBinit, muQinit, muSinit, muCinit,
                config.ConstrainMuB, config.ConstrainMuQ, config.ConstrainMuS, config.ConstrainMuC
        );
      } else {
        if (param == Tmin)
          model->ConstrainChemicalPotentials(true);
        else
          model->ConstrainChemicalPotentials(false);
      }

      model->CalculateDensities();

      varvalues->operator [](i) = param;
      double T = model->Parameters().T * 1.e3;
      Thermodynamics tTD;
      tTD.p = model->Pressure();
      tTD.e = model->EnergyDensity();
      tTD.s = model->EntropyDensity();
      tTD.I = tTD.e - 3. * tTD.p;
      tTD.pT4 = model->Pressure() / pow(T * 1.e-3, 4) / thermalfist::xMath::GeVtoifm3();
      tTD.eT4 = model->EnergyDensity() / pow(T * 1.e-3, 4) / thermalfist::xMath::GeVtoifm3();
      tTD.sT3 = model->EntropyDensity() / pow(T * 1.e-3, 3) / thermalfist::xMath::GeVtoifm3();
      tTD.IT4 = tTD.eT4 - 3. * tTD.pT4;
      tTD.nh  = model->HadronDensity();
      tTD.T = model->Parameters().T;
      tTD.muB = model->Parameters().muB;
      tTD.muQ = model->Parameters().muQ;
      tTD.muS = model->Parameters().muS;
      tTD.muC = model->Parameters().muC;
      tTD.densities = model->AllDensities();
      tTD.rhoB = model->BaryonDensity();
      tTD.rhoQ = model->ElectricChargeDensity();
      tTD.rhoS = model->StrangenessDensity();
      tTD.rhoC = model->CharmDensity();

      model->CalculateTwoParticleCorrelations();
      model->CalculateSusceptibilityMatrix();
      model->CalculateTemperatureDerivatives();

      tTD.cs2  = model->cs2();
      tTD.cVT3 = model->HeatCapacity() / pow(model->Parameters().T, 3) / thermalfist::xMath::GeVtoifm3();
      tTD.flag = true;
      paramsTD->operator [](i) = tTD;

      ChargesFluctuations flucts;

      // Higher-order fluctuations
      std::vector<double> chargesB(model->Densities().size()), chargesQ(model->Densities().size()), chargesS(model->Densities().size()), chargesC(model->Densities().size());
      std::vector<double> chchis;
      // Baryon number
      if (model->TPS()->hasBaryons()) {
        for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
          chargesB[i] = model->TPS()->Particles()[i].BaryonCharge();
        }

        chchis = model->CalculateChargeFluctuations(chargesB, 4);
        flucts.chi1B = chchis[0];
        flucts.chi2B = chchis[1];
        flucts.chi3B = chchis[2];
        flucts.chi4B = chchis[3];
      }

      // Electric charge
      if (model->TPS()->hasCharged()) {
        for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
          chargesQ[i] = model->TPS()->Particles()[i].ElectricCharge();
        }

        chchis = model->CalculateChargeFluctuations(chargesQ, 4);
        flucts.chi1Q = chchis[0];
        flucts.chi2Q = chchis[1];
        flucts.chi3Q = chchis[2];
        flucts.chi4Q = chchis[3];
      }

      // Strangeness
      if (model->TPS()->hasStrange()) {
        for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
          chargesS[i] = model->TPS()->Particles()[i].Strangeness();
        }

        chchis = model->CalculateChargeFluctuations(chargesS, 4);
        flucts.chi1S = chchis[0];
        flucts.chi2S = chchis[1];
        flucts.chi3S = chchis[2];
        flucts.chi4S = chchis[3];
      }

      // Charm
      if (model->TPS()->hasCharmed()) {
        for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
          chargesC[i] = model->TPS()->Particles()[i].Charm();
        }

        chchis = model->CalculateChargeFluctuations(chargesC, 4);
        flucts.chi1C = chchis[0];
        flucts.chi2C = chchis[1];
        flucts.chi3C = chchis[2];
        flucts.chi4C = chchis[3];
      }

      // Correlations
      flucts.chi11BQ = model->Susc(thermalfist::ConservedCharge::BaryonCharge,
        thermalfist::ConservedCharge::ElectricCharge);
      flucts.chi11BS = model->Susc(thermalfist::ConservedCharge::BaryonCharge,
        thermalfist::ConservedCharge::StrangenessCharge);
      flucts.chi11QS = model->Susc(thermalfist::ConservedCharge::ElectricCharge,
        thermalfist::ConservedCharge::StrangenessCharge);

      flucts.flag = true;

      paramsFl->operator [](i) = flucts;

      (*currentSize)++;
    }
  }
  emit calculated();
}

EquationOfStateTab::EquationOfStateTab(QWidget *parent, ThermalModelBase *modelop) :
    QWidget(parent)
{
    TPS = modelop->TPS();
    model = new ThermalModelIdeal(modelop->TPS());
    *model = *modelop;

    int index = 0;
    QString tname;

    tname = "p/T⁴";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "e/T⁴";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "s/T³";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "(e-3P)/T⁴";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "Δ";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "cs²";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "CV/T³";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "T[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "μB[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "μQ[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "μS[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "μC[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρB/T³";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρQ/T³";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρS/T³";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρC/T³";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ntot/T³";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ni/T³";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "p[MeV/fm^3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "e[MeV/fm^3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "(e-3P)[MeV/fm^3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "s[fm^-3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ntot[fm^-3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ni[fm^-3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρB[fm^-3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρQ[fm^-3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρS[fm^-3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρC[fm^-3]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₁B";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₁Q";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₁S";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₂B";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₂Q";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₂S";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₁₁BQ";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "-χ₁₁BS";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₁₁QS";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "CBS";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₃B";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₃Q";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₃S";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₄B";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₄Q";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "χ₄S";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;


    QHBoxLayout *mainLayout = new QHBoxLayout();

    QVBoxLayout *layLeft = new QVBoxLayout();

    QHBoxLayout *layLeftTop = new QHBoxLayout();
    layLeftTop->setAlignment(Qt::AlignLeft);
    CBratio = new QCheckBox(tr("Ratio"));
    CBratio->setChecked(false);
    connect(CBratio, SIGNAL(toggled(bool)), this, SLOT(replot()));
    connect(CBratio, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));

    CBxAxis = new QCheckBox(tr("Show quantity 2 on x-axis"));
    CBxAxis->setChecked(false);
    connect(CBxAxis, SIGNAL(toggled(bool)), this, SLOT(replot()));
    connect(CBxAxis, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));


    CBflipAxes = new QCheckBox(tr("Flip axes"));
    CBflipAxes->setChecked(false);
    connect(CBflipAxes, SIGNAL(toggled(bool)), this, SLOT(replot()));

    layLeftTop->addWidget(CBratio);
    layLeftTop->addWidget(CBxAxis);
    layLeftTop->addWidget(CBflipAxes);

    QGridLayout *selLay = new QGridLayout();
    selLay->setAlignment(Qt::AlignLeft);
    QLabel *labelSel = new QLabel(tr("Quantity:"));
    comboQuantity = new QComboBox();
    for(int i=0;i<index;++i) comboQuantity->addItem(paramnames[i]);
    comboQuantity->setCurrentIndex(0);
    connect(comboQuantity, SIGNAL(currentIndexChanged(int)), this, SLOT(replot()));
    connect(comboQuantity, SIGNAL(currentIndexChanged(int)), this, SLOT(modelChanged()));
    selLay->addWidget(labelSel, 0, 0);
    selLay->addWidget(comboQuantity, 0, 1);

    labelParticle = new QLabel(tr("Particle:"));
    comboParticle = new QComboBox();
    connect(comboParticle, SIGNAL(currentIndexChanged(int)), this, SLOT(replot()));
    selLay->addWidget(labelParticle, 0, 2);
    selLay->addWidget(comboParticle, 0, 3);

    labelFeeddown = new QLabel(tr("Feeddown:"));
    comboFeeddown= new QComboBox();
    comboFeeddown->addItem(tr("Primordial"));
    comboFeeddown->addItem(tr("Stability flags"));
    comboFeeddown->addItem(tr("Strong+EM+weak"));
    comboFeeddown->addItem(tr("Strong+EM"));
    comboFeeddown->addItem(tr("Strong"));
    comboFeeddown->setCurrentIndex(0);
    connect(comboFeeddown, SIGNAL(currentIndexChanged(int)), this, SLOT(replot()));
    selLay->addWidget(labelFeeddown, 0, 4);
    selLay->addWidget(comboFeeddown, 0, 5);

    QLabel *labelSel2 = new QLabel(tr("Quantity 2:"));
    comboQuantity2 = new QComboBox();
    for (int i = 0; i<index; ++i) comboQuantity2->addItem(paramnames[i]);
    comboQuantity2->setCurrentIndex(0);
    connect(comboQuantity2, SIGNAL(currentIndexChanged(int)), this, SLOT(replot()));
    connect(comboQuantity2, SIGNAL(currentIndexChanged(int)), this, SLOT(modelChanged()));
    selLay->addWidget(labelSel2, 1, 0);
    selLay->addWidget(comboQuantity2, 1, 1);

    labelParticle2 = new QLabel(tr("Particle:"));
    comboParticle2 = new QComboBox();
    connect(comboParticle2, SIGNAL(currentIndexChanged(int)), this, SLOT(replot()));
    selLay->addWidget(labelParticle2, 1, 2);
    selLay->addWidget(comboParticle2, 1, 3);

    labelFeeddown2 = new QLabel(tr("Feeddown:"));
    comboFeeddown2 = new QComboBox();
    comboFeeddown2->addItem(tr("Primordial"));
    comboFeeddown2->addItem(tr("Stability flags"));
    comboFeeddown2->addItem(tr("Strong+EM+weak"));
    comboFeeddown2->addItem(tr("Strong+EM"));
    comboFeeddown2->addItem(tr("Strong"));
    comboFeeddown2->setCurrentIndex(0);
    connect(comboFeeddown2, SIGNAL(currentIndexChanged(int)), this, SLOT(replot()));
    selLay->addWidget(labelFeeddown2, 1, 4);
    selLay->addWidget(comboFeeddown2, 1, 5);

    plotDependence = new QCustomPlot();
    plotDependence->xAxis->setLabel("T (MeV)");
    plotDependence->yAxis->setLabel(comboQuantity->currentText());

    plotDependence->axisRect()->setupFullAxesBox();

    plotDependence->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(plotDependence, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequest(QPoint)));

    buttonEoSTable = new QPushButton(tr("Tabulated EoS..."));
    connect(buttonEoSTable, SIGNAL(clicked()), this, SLOT(showEoSTable()));

    //layLeft->addWidget(CBratio, 0, Qt::AlignLeft);
    layLeft->addLayout(layLeftTop);
    layLeft->addLayout(selLay);
    layLeft->addWidget(plotDependence);
    layLeft->addWidget(buttonEoSTable, 0, Qt::AlignLeft);

    QVBoxLayout *layRight = new QVBoxLayout();
    layRight->setContentsMargins(15, 0, 0, 0);
    layRight->setAlignment(Qt::AlignTop);

    QGroupBox *grModelConfig = new QGroupBox(tr("HRG model configuration:"));

    configWidget = new ModelConfigWidget(NULL, model);
    connect(configWidget, SIGNAL(changed()), this, SLOT(modelChanged()));

    QHBoxLayout *layModelConfig = new QHBoxLayout();
    layModelConfig->setAlignment(Qt::AlignLeft);
    layModelConfig->addWidget(configWidget);

    grModelConfig->setLayout(layModelConfig);

    QHBoxLayout* layMode = new QHBoxLayout();
    layMode->setAlignment(Qt::AlignLeft);
    comboMode = new QComboBox();
    comboMode->addItem(tr("Const. μ"));
    comboMode->addItem(tr("Const. μ/T"));
    comboMode->addItem(tr("Const. T"));
    comboMode->addItem(tr("Magnetic field dependence"));
    comboMode->setCurrentIndex(0);
    connect(comboMode, SIGNAL(currentIndexChanged(int)), this, SLOT(modelChanged()));
    layMode->addWidget(new QLabel(tr("Mode:")));
    layMode->addWidget(comboMode);

    QGroupBox *grRange = new QGroupBox(tr("Parameter range:"));

    QVBoxLayout *layParamBounds = new QVBoxLayout();
    QGridLayout *layBounds = new QGridLayout();
    layBounds->setAlignment(Qt::AlignLeft);

    labelmuB = new QLabel("μ<sub>B</sub> (MeV):");
    spinmuB = new QDoubleSpinBox();
    spinmuB->setDecimals(3);
    spinmuB->setMinimum(-10000.);
    spinmuB->setMaximum(10000.);
    spinmuB->setValue(0.);

    labelmuQ = new QLabel("μ<sub>Q</sub> (MeV):");
    spinmuQ = new QDoubleSpinBox();
    spinmuQ->setDecimals(3);
    spinmuQ->setMinimum(-10000.);
    spinmuQ->setMaximum(10000.);
    spinmuQ->setValue(0.);

    labelmuS = new QLabel("μ<sub>S</sub> (MeV):");
    spinmuS = new QDoubleSpinBox();
    spinmuS->setDecimals(3);
    spinmuS->setMinimum(-10000.);
    spinmuS->setMaximum(10000.);
    spinmuS->setValue(0.);

    labelTMin = new QLabel(tr("T<sub>min</sub> (MeV):"));
    spinTMin = new QDoubleSpinBox();
    spinTMin->setDecimals(3);
    spinTMin->setMinimum(0.);
    spinTMin->setMaximum(10000.);
    spinTMin->setValue(100.);

    labelTMax = new QLabel(tr("T<sub>max</sub> (MeV):"));
    spinTMax = new QDoubleSpinBox();
    spinTMax->setDecimals(3);
    spinTMax->setMinimum(0.01);
    spinTMax->setMaximum(10000.);
    spinTMax->setValue(200.);

    labeldT = new QLabel(tr("∆T (MeV):"));
    spindT = new QDoubleSpinBox();
    spindT->setDecimals(3);
    spindT->setMinimum(0.);
    spindT->setMaximum(10000.);
    spindT->setValue(5.);

    labelEMin = new QLabel(tr("eB<sub>min</sub> (GeV<sup>2</sup>):"));
    spinEMin = new QDoubleSpinBox();
    spinEMin->setDecimals(3);
    spinEMin->setMinimum(0.);
    spinEMin->setMaximum(100.);
    spinEMin->setSingleStep(0.01);
    spinEMin->setValue(0.);

    labelEMax = new QLabel(tr("eB<sub>max</sub> (GeV<sup>2</sup>):"));
    spinEMax = new QDoubleSpinBox();
    spinEMax->setDecimals(3);
    spinEMax->setMinimum(0.);
    spinEMax->setMaximum(100.);
    spinEMax->setSingleStep(0.01);
    spinEMax->setValue(0.16);

    labeldE = new QLabel(tr("∆eB (GeV<sup>2</sup>):"));
    spindE = new QDoubleSpinBox();
    spindE->setDecimals(3);
    spindE->setMinimum(0.);
    spindE->setMaximum(100.);
    spindE->setSingleStep(0.01);
    spindE->setValue(0.01);

    labelTaux = new QLabel(tr("T (MeV):"));
    spinTaux = new QDoubleSpinBox();
    spinTaux->setMinimum(0.01);
    spinTaux->setMaximum(10000.);
    spinTaux->setValue(155.);



    layBounds->addWidget(labelTMin, 0, 0);
    layBounds->addWidget(spinTMin, 0, 1);
    layBounds->addWidget(labelTMax, 0, 2);
    layBounds->addWidget(spinTMax, 0, 3);
    layBounds->addWidget(labeldT, 0, 4);
    layBounds->addWidget(spindT, 0, 5);
    layBounds->addWidget(labelEMin, 1, 0);
    layBounds->addWidget(spinEMin, 1, 1);
    layBounds->addWidget(labelEMax, 1, 2);
    layBounds->addWidget(spinEMax, 1, 3);
    layBounds->addWidget(labeldE, 1, 4);
    layBounds->addWidget(spindE, 1, 5);
    layBounds->addWidget(labelmuB, 2, 0);
    layBounds->addWidget(spinmuB, 2, 1);
    layBounds->addWidget(labelmuQ, 2, 2);
    layBounds->addWidget(spinmuQ, 2, 3);
    layBounds->addWidget(labelmuS, 2, 4);
    layBounds->addWidget(spinmuS, 2, 5);
    layBounds->addWidget(labelTaux, 3, 0);
    layBounds->addWidget(spinTaux, 3, 1);

    labelConstr = new QLabel(tr(""));


    layParamBounds->addLayout(layBounds);
    layParamBounds->addWidget(labelConstr, 0, Qt::AlignLeft);

    grRange->setLayout(layParamBounds);

    buttonCalculate = new QPushButton(tr("Calculate"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));

    layRight->addWidget(grModelConfig);
    layRight->addLayout(layMode);
    layRight->addWidget(grRange);
    layRight->addWidget(buttonCalculate, 0, Qt::AlignLeft);

    mainLayout->addLayout(layLeft, 1);
    mainLayout->addLayout(layRight);

    setLayout(mainLayout);

    fCurrentSize = 0;
    fStop = 1;
    fRunning = false;

    calcTimer = new QTimer(this);
    connect(calcTimer, SIGNAL(timeout()), this, SLOT(replot()));

    readLatticeData();

    fillParticleLists();

    modelChanged();
    
    replot();

    cpath = QApplication::applicationDirPath();
}

EquationOfStateTab::~EquationOfStateTab() {
    if (model!=NULL) delete model;
}

ThermalModelConfig EquationOfStateTab::getConfig()
{
  ThermalModelConfig ret = configWidget->updatedConfig();

  ret.T   = spinTMin->value() * 1.e-3;
  ret.muB = spinmuB->value()  * 1.e-3;
  ret.muQ = spinmuQ->value()  * 1.e-3;
  ret.muS = spinmuS->value()  * 1.e-3;
  ret.muC = 0.;
  ret.gq = 1.;
  ret.gS = 1.;
  ret.gC = 1.;
  ret.VolumeR = ret.VolumeRSC = 8.;
  ret.B = 0;
  ret.Q = 0;
  ret.S = 0;
  ret.C = 0;

  ret.ComputeFluctations = false; //checkFluctuations->isChecked();
  ret.ResetMus = true;

  return ret;
}

void EquationOfStateTab::calculate() {
    if (!fRunning) {
      ThermalModelBase *modelnew;

      ThermalModelConfig config = getConfig();

      if (config.ModelType == ThermalModelConfig::DiagonalEV) {
        modelnew = new ThermalModelEVDiagonal(model->TPS());
      }
      else if (config.ModelType == ThermalModelConfig::CrosstermsEV) {
        modelnew = new ThermalModelEVCrossterms(model->TPS());
      }
      else if (config.ModelType == ThermalModelConfig::QvdW) {
        modelnew = new ThermalModelVDW(model->TPS());
      }
      else if (config.ModelType == ThermalModelConfig::RealGas) {
        modelnew = new ThermalModelRealGas(model->TPS());
      }
      else {
        modelnew = new ThermalModelIdeal(model->TPS());
      }

      if (model != NULL)
        modelnew->SetNormBratio(config.RenormalizeBR);

      configWidget->setModel(modelnew);

      if (model != NULL)
        delete model;

      model = modelnew;

      model->SetTemperature(config.T);
      model->SetBaryonChemicalPotential(config.muB);
      model->SetElectricChemicalPotential(config.muQ);
      model->SetStrangenessChemicalPotential(config.muS);
      model->SetCharmChemicalPotential(config.muC);
      model->SetGammaq(config.gq);
      model->SetGammaS(config.gS);
      model->SetGammaC(config.gC);
      model->SetVolumeRadius(config.VolumeR);

      SetThermalModelConfiguration(model, config);

      SetThermalModelInteraction(model, config);

      if (comboMode->currentIndex() == 2) {
        model->ConstrainMuB(false);
      }

      paramsTD.resize(0);
      paramsFl.resize(0);
      varvalues.resize(0);

      double pmin = spinTMin->value(), pmax = spinTMax->value(), dp = spindT->value();
      if (comboMode->currentIndex() == 3) {
        pmin = spinEMin->value();
        pmax = spinEMax->value();
        dp = spindE->value();
      }

      int iters = static_cast<int>((pmax - pmin) / dp) + 1;

      fCurrentSize = 0;

      paramsTD.resize(iters);
      paramsFl.resize(iters);
      varvalues.resize(iters);

      fStop = 0;

      std::vector<double> mus = {spinmuB->value(), spinmuQ->value(), spinmuS->value(), spinTaux->value()};

      EoSWorker *wrk = new EoSWorker(model, config, pmin, pmax, dp,
                                     mus, comboMode->currentIndex(),
                                    &paramsTD, &paramsFl, &varvalues, &fCurrentSize, &fStop, this);
      connect(wrk, SIGNAL(calculated()), this, SLOT(finalize()));
      connect(wrk, SIGNAL(finished()), wrk, SLOT(deleteLater()));

      wrk->start();

      buttonCalculate->setText(tr("Stop"));
      fRunning = true;

      calcTimer->start(10);
    }
    else {
      fStop = 1;
    }
}

void EquationOfStateTab::finalize() {
    calcTimer->stop();
    buttonCalculate->setText(tr("Calculate"));
    fRunning = false;
    replot();
}

void EquationOfStateTab::modelChanged()
{
  if (CBratio->isChecked() || CBxAxis->isChecked()) {
    comboQuantity2->setEnabled(true);
  }
  else {
    comboQuantity2->setEnabled(false);
  }

  if (comboQuantity->currentText() == "ni/T³" || comboQuantity->currentText() == "ni[fm^-3]") {
    labelParticle->setVisible(true);
    comboParticle->setVisible(true);
    labelFeeddown->setVisible(true);
    comboFeeddown->setVisible(true);
  }
  else {
    labelParticle->setVisible(false);
    comboParticle->setVisible(false);
    labelFeeddown->setVisible(false);
    comboFeeddown->setVisible(false);
  }

  if ((comboQuantity2->currentText() == "ni/T³" || comboQuantity2->currentText() == "ni[fm^-3]")
    && (CBratio->isChecked() || CBxAxis->isChecked())) {
    labelParticle2->setVisible(true);
    comboParticle2->setVisible(true);
    labelFeeddown2->setVisible(true);
    comboFeeddown2->setVisible(true);
  }
  else {
    labelParticle2->setVisible(false);
    comboParticle2->setVisible(false);
    labelFeeddown2->setVisible(false);
    comboFeeddown2->setVisible(false);
  }

  if (configWidget->currentConfig.ConstrainMuB && comboMode->currentIndex() != 2) {
    spinmuB->setEnabled(false);
  }
  else {
    spinmuB->setEnabled(true);
  }

  spinmuQ->setEnabled(!configWidget->currentConfig.ConstrainMuQ);
  spinmuS->setEnabled(!configWidget->currentConfig.ConstrainMuS);

  if (comboMode->currentIndex() == 2) {
    labelmuB->setText("T (MeV)");
    labelTMin->setText("μB<sub>min</sub> (MeV)");
    labelTMax->setText("μB<sub>max</sub> (MeV)");
    labeldT->setText("∆μB (MeV)");
    labelmuQ->setText("μ<sub>Q</sub> (MeV)");
    labelmuS->setText("μ<sub>S</sub> (MeV)");
  }
  else {
    labelmuB->setText("μ<sub>B</sub> (MeV)");
    labelmuQ->setText("μ<sub>Q</sub> (MeV)");
    labelmuS->setText("μ<sub>S</sub> (MeV)");
    if (comboMode->currentIndex() == 1) {
      labelmuB->setText("μ<sub>B</sub>/T");
      labelmuQ->setText("μ<sub>Q</sub>/T");
      labelmuS->setText("μ<sub>S</sub>/T");
    }
    labelTMin->setText("T<sub>min</sub> (MeV)");
    labelTMax->setText("T<sub>max</sub> (MeV)");
    labeldT->setText("∆T (MeV)");
  }


  if (comboMode->currentIndex() != 3) {
    labelEMin->setVisible(false);
    labelEMax->setVisible(false);
    labeldE->setVisible(false);
    spinEMin->setVisible(false);
    spinEMax->setVisible(false);
    spindE->setVisible(false);
    labelTMin->setVisible(true);
    labelTMax->setVisible(true);
    labeldT->setVisible(true);
    spinTMin->setVisible(true);
    spinTMax->setVisible(true);
    spindT->setVisible(true);
    labelTaux->setVisible(false);
    spinTaux->setVisible(false);
  } else {
    labelEMin->setVisible(true);
    labelEMax->setVisible(true);
    labeldE->setVisible(true);
    spinEMin->setVisible(true);
    spinEMax->setVisible(true);
    spindE->setVisible(true);
    labelTMin->setVisible(false);
    labelTMax->setVisible(false);
    labeldT->setVisible(false);
    spinTMin->setVisible(false);
    spinTMax->setVisible(false);
    spindT->setVisible(false);
    labelTaux->setVisible(true);
    spinTaux->setVisible(true);
  }

  QString constr = "";

  if (configWidget->currentConfig.ConstrainMuB && model->TPS()->hasBaryons() && comboMode->currentIndex() != 2) {
    constr += "S/B = " + QString::number(configWidget->currentConfig.SoverB, 'g', 2);
  }

  if (configWidget->currentConfig.ConstrainMuQ && model->TPS()->hasCharged()) {
    if (constr != "")
      constr += ", ";
    constr += "Q/B = " + QString::number(configWidget->currentConfig.QoverB, 'g', 2);
  }

  if (configWidget->currentConfig.ConstrainMuS && model->TPS()->hasStrange()) {
    if (constr != "")
      constr += ", ";
    constr += "S = 0";
  }

  if (configWidget->currentConfig.ConstrainMuC && model->TPS()->hasCharmed()) {
    if (constr != "")
      constr += ", ";
    constr += "C = 0";
  }

  if (constr != "") {
    labelConstr->setText(constr);
    labelConstr->setVisible(true);
  }
  else {
    labelConstr->setVisible(false);
  }

  labelmuB->setVisible(model->TPS()->hasBaryons());
  spinmuB->setVisible(model->TPS()->hasBaryons());
  labelmuQ->setVisible(model->TPS()->hasCharged());
  spinmuQ->setVisible(model->TPS()->hasCharged());
  labelmuS->setVisible(model->TPS()->hasStrange());
  spinmuS->setVisible(model->TPS()->hasStrange());

  configWidget->comboEnsemble->setEnabled(false);
}

void EquationOfStateTab::resetTPS()
{
  model->ChangeTPS(model->TPS());
  fillParticleLists();

  configWidget->setModel(model);
  configWidget->currentConfig.vdWparams = QvdWParameters::GetParameters(model->TPS(), &configWidget->currentConfig);
}

void EquationOfStateTab::plotLatticeData()
{
  int index = comboQuantity->currentIndex();
  if (index >= 0 && index < paramnames.size()) {
    int graphStart = plotDependence->graphCount();
    QString paramname;

    if (CBratio->isChecked()) {
      int index2 = comboQuantity2->currentIndex();
      if (!(index2 >= 0 && index2 < paramnames.size()))
        return;
      paramname = paramnames[index] + "/" + paramnames[index2];
    }
    else {
      paramname = paramnames[index];
    }

    // WB data
    if (mapWB.count(paramname) > 0
    || (CBratio->isChecked()
    && mapWB.count(paramnames[index]) && mapWB.count(paramnames[comboQuantity2->currentIndex()])
    && dataWBx[mapWB[paramnames[index]]].size() == dataWBx[mapWB[paramnames[comboQuantity2->currentIndex()]]].size())
    ) {
      
      int indWB = 0;
      if (mapWB.count(paramname))
        indWB = mapWB[paramname];
      int indWB2 = 0;
      if (CBratio->isChecked()) {
        int index2 = comboQuantity2->currentIndex();
        indWB  = mapWB[paramnames[index]];
        indWB2 = mapWB[paramnames[index2]];
        if (dataWBx[indWB].size() != dataWBx[indWB2].size())
          return;
      }
      plotDependence->addGraph();
      plotDependence->graph(graphStart)->setName("LQCD (Wuppertal-Budapest)");
      plotDependence->graph(graphStart)->setPen(QPen(QColor(255, 255, 255, 0)));
      plotDependence->graph(graphStart)->setBrush(QBrush(QColor(0, 0, 255, 120)));
      plotDependence->addGraph();
      plotDependence->legend->removeItem(plotDependence->legend->itemCount() - 1); // don't show two confidence band graphs in legend
      plotDependence->graph(graphStart+1)->setName("LQCD (Wuppertal-Budapest)");
      plotDependence->graph(graphStart+1)->setPen(QPen(QColor(255, 255, 255, 0)));
      plotDependence->graph(graphStart)->setChannelFillGraph(plotDependence->graph(graphStart+1));

      double tmin = 1.e5, tmax = 0.;
      QVector<double> x0(dataWBx[indWB].size());
      QVector<double> yConfUpper(dataWBx[indWB].size()), yConfLower(dataWBx[indWB].size());
      for (int i = 0; i < x0.size(); ++i) {
        x0[i] = dataWBx[indWB][i];
        if (mapWB.count(paramname)) {
          yConfUpper[i] = dataWBy[indWB][i] + dataWByerrp[indWB][i];
          yConfLower[i] = dataWBy[indWB][i] - dataWByerrm[indWB][i];
        }
        else if (CBratio->isChecked()) {
          double mean = dataWBy[indWB][i] / dataWBy[indWB2][i];
          double error = mean * sqrt((dataWByerrp[indWB][i]/dataWBy[indWB][i])*(dataWByerrp[indWB][i]/dataWBy[indWB][i])
                  + (dataWByerrp[indWB2][i]/dataWBy[indWB2][i])*(dataWByerrp[indWB2][i]/dataWBy[indWB2][i]));

          yConfUpper[i] = mean + error;
          yConfLower[i] = mean - error;
        }
        if (x0[i] >= spinTMin->value() && x0[i] <= spinTMax->value()) {
          tmin = std::min(tmin, yConfLower[i]);
          tmax = std::max(tmax, yConfUpper[i]);
        }
      }

      plotDependence->graph(graphStart)->setData(x0, yConfUpper);
      plotDependence->graph(graphStart+1)->setData(x0, yConfLower);

      plotDependence->yAxis->setRange(0.*tmin, 1.1*tmax);

      graphStart += 2;
    }

    // HotQCD data
    if (mapHotQCD.count(paramname) > 0 || (CBratio->isChecked()
      && mapHotQCD.count(paramnames[index]) && mapHotQCD.count(paramnames[comboQuantity2->currentIndex()])
      && dataHotQCDx[mapHotQCD[paramnames[index]]].size() == dataHotQCDx[mapWB[paramnames[comboQuantity2->currentIndex()]]].size())
      )
    {
      int indHotQCD = 0;
      if (mapHotQCD.count(paramname))
        indHotQCD = mapHotQCD[paramname];
      int indHotQCD2 = 0;
      if (CBratio->isChecked()) {
        int index2 = comboQuantity2->currentIndex();
        indHotQCD  = mapHotQCD[paramnames[index]];
        indHotQCD2 = mapHotQCD[paramnames[index2]];
        if (dataHotQCDx[indHotQCD].size() != dataHotQCDx[indHotQCD2].size())
          return;
      }

      plotDependence->addGraph();
      plotDependence->graph(graphStart)->setName("LQCD (HotQCD)");
      plotDependence->graph(graphStart)->setPen(QPen(QColor(255, 255, 255, 0)));
      plotDependence->graph(graphStart)->setBrush(QBrush(QColor(0, 255, 0, 120)));
      plotDependence->addGraph();
      plotDependence->legend->removeItem(plotDependence->legend->itemCount() - 1); // don't show two confidence band graphs in legend
      plotDependence->graph(graphStart + 1)->setName("LQCD (HotQCD)");
      plotDependence->graph(graphStart + 1)->setPen(QPen(QColor(255, 255, 255, 0)));
      plotDependence->graph(graphStart)->setChannelFillGraph(plotDependence->graph(graphStart + 1));

      double tmin = 1.e5, tmax = 0.;
      QVector<double> x0(dataHotQCDx[indHotQCD].size());
      QVector<double> yConfUpper(dataHotQCDx[indHotQCD].size()), yConfLower(dataHotQCDx[indHotQCD].size());
      for (int i = 0; i < x0.size(); ++i) {
        x0[i] = dataHotQCDx[indHotQCD][i];

        if (mapHotQCD.count(paramname)) {
          yConfUpper[i] = dataHotQCDy[indHotQCD][i] + dataHotQCDyerrp[indHotQCD][i];
          yConfLower[i] = dataHotQCDy[indHotQCD][i] - dataHotQCDyerrm[indHotQCD][i];
        }
        else if (CBratio->isChecked()) {
          double mean = dataHotQCDy[indHotQCD][i] / dataHotQCDy[indHotQCD2][i];
          double error = mean * sqrt((dataHotQCDyerrp[indHotQCD][i]/dataHotQCDy[indHotQCD][i])*(dataHotQCDyerrp[indHotQCD][i]/dataHotQCDy[indHotQCD][i])
                                     + (dataHotQCDyerrp[indHotQCD2][i]/dataHotQCDy[indHotQCD2][i])*(dataHotQCDyerrp[indHotQCD2][i]/dataHotQCDy[indHotQCD2][i]));

          yConfUpper[i] = mean + error;
          yConfLower[i] = mean - error;
        }

        if (x0[i] >= spinTMin->value() && x0[i] <= spinTMax->value()) {
          tmin = std::min(tmin, yConfLower[i]);
          tmax = std::max(tmax, yConfUpper[i]);
        }
      }

      plotDependence->graph(graphStart)->setData(x0, yConfUpper);
      plotDependence->graph(graphStart + 1)->setData(x0, yConfLower);

      plotDependence->yAxis->setRange(0.*tmin, 1.1*tmax);
    }
  }
}

void EquationOfStateTab::showEoSTable() {
  recomputeCalcTable();

  CalculationTableDialog dialog(this, calcTable);
  dialog.setWindowFlags(Qt::Window);
  dialog.showMaximized();
  dialog.exec();
}

void EquationOfStateTab::fillParticleLists()
{
  for (int ii = 0; ii < 2; ++ii) {
    QComboBox *combo;
    if (ii == 0)
      combo = comboParticle;
    else
      combo = comboParticle2;

    int tind = combo->currentIndex();
    combo->clear();
    for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
      const ThermalParticle &part = model->TPS()->Particles()[i];
      combo->addItem(QString::number(part.PdgId()) + " - " + QString::fromStdString(part.Name()));
    }
    if (tind >= 0 && tind < combo->count())
      combo->setCurrentIndex(tind);
    else
      combo->setCurrentIndex(0);
  }
}

void EquationOfStateTab::contextMenuRequest(QPoint pos)
{
  QMenu *menu = new QMenu(this);
  menu->setAttribute(Qt::WA_DeleteOnClose);

  menu->addAction("Save as pdf...", this, SLOT(saveAsPdf()));
  menu->addAction("Save as png...", this, SLOT(saveAsPng()));
  menu->addAction("Write computed values to file...", this, SLOT(saveAsAscii()));

  menu->popup(plotDependence->mapToGlobal(pos));
}

void EquationOfStateTab::saveAsPdf()
{
  saveAs(0);
}

void EquationOfStateTab::saveAsPng()
{
  saveAs(1);
}

void EquationOfStateTab::saveAsAscii()
{
  saveAs(2);
}

void EquationOfStateTab::saveAs(int type)
{
  QString tname = plotDependence->yAxis->label();
  for (int i = 0; i < tname.size(); ++i) {
    if (tname[i] == '/')
      tname[i] = '.';
  }

  QVector<QString> exts;
  exts.push_back("pdf");
  exts.push_back("png");
  exts.push_back("dat");

  QString listpathprefix = cpath + "/" + tname + "." + exts[type];
  QString path = QFileDialog::getSaveFileName(this, "Save plot as " + exts[type], listpathprefix, "*." + exts[type]);
  if (path.length()>0)
  {
    if (type == 0)
      plotDependence->savePdf(path, plotDependence->width(), plotDependence->height());
    else if (type == 1)
      plotDependence->savePng(path, plotDependence->width(), plotDependence->height());
    else {
      std::vector<double> yvalues;

      if (CBratio->isChecked() && !CBxAxis->isChecked()) {
        yvalues = getValuesRatio(comboQuantity->currentIndex(), comboQuantity2->currentIndex());
      }
      else {
        yvalues = getValues(comboQuantity->currentIndex());
      }

      std::vector<double> xvalues = varvalues;
      if (CBxAxis->isChecked()) {
        xvalues = getValues(comboQuantity2->currentIndex(), 1);
      }

      QFile fout(path);

      if (yvalues.size() > 0 && fout.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&fout);
        out.setFieldWidth(15);
        out.setFieldAlignment(QTextStream::AlignLeft);
        if (!CBflipAxes->isChecked()) {
          out << plotDependence->xAxis->label();
          out << plotDependence->yAxis->label();
        } else {
          out << plotDependence->yAxis->label();
          out << plotDependence->xAxis->label();
        }
        out << qSetFieldWidth(0) << Qt::endl << qSetFieldWidth(15);
        for (int i = 0; i < yvalues.size(); ++i) {
          out.setFieldWidth(15);
          out.setFieldAlignment(QTextStream::AlignLeft);
          out << xvalues[i] << yvalues[i];
          out << qSetFieldWidth(0) << Qt::endl << qSetFieldWidth(15);
        }
      }
    }
    
    QFileInfo saved(path);
    cpath = saved.absolutePath();
  }
}

std::vector<double> EquationOfStateTab::getValues(int index, int num)
{
  int tsize = fCurrentSize;
  std::vector<double> ret(tsize, 0.);

  if (index >= 0 && index < paramnames.size() && paramnames[index] == "ni/T³") {
    int pid = comboParticle->currentIndex();
    int feed = comboFeeddown->currentIndex();
    if (num == 1) {
      pid = comboParticle2->currentIndex();
      feed = comboFeeddown2->currentIndex();
    }

    for (int j = 0; j < tsize; ++j) {
      ret[j] = paramsTD[j].densities[feed][pid];
      ret[j] *= 1. / pow(paramsTD[j].T, 3) / xMath::GeVtoifm3();
    }

    return ret;
  }

  if (index >= 0 && index < paramnames.size() && paramnames[index] == "ni[fm^-3]") {
    int pid = comboParticle->currentIndex();
    int feed = comboFeeddown->currentIndex();
    if (num == 1) {
      pid = comboParticle2->currentIndex();
      feed = comboFeeddown2->currentIndex();
    }

    for (int j = 0; j < tsize; ++j) {
      ret[j] = paramsTD[j].densities[feed][pid];
    }

    return ret;
  }
  
  //for (int i = 0; i < paramnames.size(); ++i) {
    int tind = parammap["p/T⁴"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].pT4;
      }
    }

    tind = parammap["e/T⁴"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].eT4;
      }
    }

    tind = parammap["(e-3P)/T⁴"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].IT4;
      }
    }

    tind = parammap["s/T³"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].sT3;
      }
    }

    tind = parammap["T[MeV]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].T * 1.e3;
      }
    }

    tind = parammap["μB[MeV]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].muB * 1.e3;
      }
    }

    tind = parammap["μQ[MeV]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].muQ * 1.e3;
      }
    }

    tind = parammap["μS[MeV]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].muS * 1.e3;
      }
    }

    tind = parammap["μC[MeV]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].muC * 1.e3;
      }
    }

    tind = parammap["ntot/T³"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].nh;
        ret[j] *= 1. / pow(paramsTD[j].T, 3) / xMath::GeVtoifm3();
      }
    }

    tind = parammap["ρB/T³"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoB;
        ret[j] *= 1. / pow(paramsTD[j].T, 3) / xMath::GeVtoifm3();
      }
    }

    tind = parammap["ρQ/T³"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoQ;
        ret[j] *= 1. / pow(paramsTD[j].T, 3) / xMath::GeVtoifm3();
      }
    }

    tind = parammap["ρS/T³"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoS;
        ret[j] *= 1. / pow(paramsTD[j].T, 3) / xMath::GeVtoifm3();
      }
    }

    tind = parammap["ρC/T³"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoC;
        ret[j] *= 1. / pow(paramsTD[j].T, 3) / xMath::GeVtoifm3();
      }
    }

    tind = parammap["p[MeV/fm^3]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].p * 1.e3;
      }
    }

    tind = parammap["e[MeV/fm^3]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].e * 1.e3;
      }
    }

    tind = parammap["(e-3P)[MeV/fm^3]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].I * 1.e3;
      }
    }

    tind = parammap["Δ"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].I / paramsTD[j].e / 3.;
      }
    }

    tind = parammap["cs²"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].cs2;
      }
    }

    tind = parammap["CV/T³"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].cVT3;
      }
    }

    tind = parammap["s[fm^-3]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].s;
      }
    }

    tind = parammap["ntot[fm^-3]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].nh;
      }
    }

    tind = parammap["ρB[fm^-3]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoB;
      }
    }

    tind = parammap["ρQ[fm^-3]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoQ;
      }
    }

    tind = parammap["ρS[fm^-3]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoS;
      }
    }

    tind = parammap["ρC[fm^-3]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoC;
      }
    }

    tind = parammap["χ₁B"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi1B;
      }
    }

    tind = parammap["χ₁Q"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi1Q;
      }
    }

    tind = parammap["χ₁S"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi1S;
      }
    }

    tind = parammap["χ₂B"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi2B;
      }
    }

    tind = parammap["χ₂Q"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi2Q;
      }
    }

    tind = parammap["χ₂S"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi2S;
      }
    }

    tind = parammap["χ₁₁BQ"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi11BQ;
      }
    }

    tind = parammap["-χ₁₁BS"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = -paramsFl[j].chi11BS;
      }
    }

    tind = parammap["χ₁₁QS"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi11QS;
      }
    }

    tind = parammap["CBS"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = -3. * paramsFl[j].chi11BS / paramsFl[j].chi2S;
      }
    }

    tind = parammap["χ₃B"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi3B;
      }
    }

    tind = parammap["χ₃Q"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi3Q;
      }
    }

    tind = parammap["χ₃S"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi3S;
      }
    }

    tind = parammap["χ₄B"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi4B;
      }
    }

    tind = parammap["χ₄Q"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi4Q;
      }
    }

    tind = parammap["χ₄S"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsFl[j].chi4S;
      }
    }
  //}

  return ret;
}

std::vector<double> EquationOfStateTab::getValuesRatio(int index, int index2)
{
  std::vector<double> ret1 = getValues(index, 0), ret2 = getValues(index2, 1);
  std::vector<double> ret;
  for (int i = 0; i < ret1.size(); ++i)
    ret.push_back(ret1[i] / ret2[i]);
  return ret;
}

void EquationOfStateTab::readLatticeData()
{
  int sizeWB = 0, sizeHotQCD = 0;
  
  ifstream fin;
  
  // WB thermodynamics
  fin.open((string(ThermalFIST_INPUT_FOLDER) + "/lqcd/WB-EoS.dat").c_str());

  if (fin.is_open()) {
    QString tname;

    tname = "p/T⁴";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "e/T⁴";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "(e-3P)/T⁴";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "s/T³";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "cs²";
    mapWB[tname] = sizeWB;
    sizeWB++;

    dataWBx.resize(sizeWB);
    dataWBy.resize(sizeWB);
    dataWByerrp.resize(sizeWB);
    dataWByerrm.resize(sizeWB);

    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1 || elems[0].size() == 0)
        continue;

      istringstream iss(elems[0]);

      double T, IT4, IT4err, pT4, pT4err, eT4, eT4err, sT3, sT3err, cs2, cs2err;

      if (iss >> T >> IT4 >> IT4err >> pT4 >> pT4err >> eT4 >> eT4err >> sT3 >> sT3err >> cs2 >> cs2err) {
        
        dataWBx[mapWB["p/T⁴"]].push_back(T);
        dataWBy[mapWB["p/T⁴"]].push_back(pT4);
        dataWByerrp[mapWB["p/T⁴"]].push_back(pT4err);
        dataWByerrm[mapWB["p/T⁴"]].push_back(pT4err);

        dataWBx[mapWB["e/T⁴"]].push_back(T);
        dataWBy[mapWB["e/T⁴"]].push_back(eT4);
        dataWByerrp[mapWB["e/T⁴"]].push_back(eT4err);
        dataWByerrm[mapWB["e/T⁴"]].push_back(eT4err);

        dataWBx[mapWB["(e-3P)/T⁴"]].push_back(T);
        dataWBy[mapWB["(e-3P)/T⁴"]].push_back(IT4);
        dataWByerrp[mapWB["(e-3P)/T⁴"]].push_back(IT4err);
        dataWByerrm[mapWB["(e-3P)/T⁴"]].push_back(IT4err);

        dataWBx[mapWB["s/T³"]].push_back(T);
        dataWBy[mapWB["s/T³"]].push_back(sT3);
        dataWByerrp[mapWB["s/T³"]].push_back(sT3err);
        dataWByerrm[mapWB["s/T³"]].push_back(sT3err);

        dataWBx[mapWB["cs²"]].push_back(T);
        dataWBy[mapWB["cs²"]].push_back(cs2);
        dataWByerrp[mapWB["cs²"]].push_back(cs2err);
        dataWByerrm[mapWB["cs²"]].push_back(cs2err);
      }
    }

    fin.close();
  }

  // WB chi2
  fin.open((string(ThermalFIST_INPUT_FOLDER) + "/lqcd/WB-chi2-1112.4416.dat").c_str());

  if (fin.is_open()) {
    QString tname;

    tname = "χ₂B";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "χ₂Q";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "χ₂S";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "CBS";
    mapWB[tname] = sizeWB;
    sizeWB++;

    dataWBx.resize(sizeWB);
    dataWBy.resize(sizeWB);
    dataWByerrp.resize(sizeWB);
    dataWByerrm.resize(sizeWB);

    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1 || elems[0].size() == 0)
        continue;

      istringstream iss(elems[0]);

      double T;
      std::string chi2I, chi2U, chi2S, chi2Q, chi2B, chi11us, CBS;

      if (iss >> T >> chi2I >> chi2U >> chi2Q >> chi2S >> chi2B >> chi11us >> CBS) {

        vector<string> abc;

        abc = CuteHRGHelper::split(chi2B, '(');
        if (abc.size() >= 2) {
          dataWBx[mapWB["χ₂B"]].push_back(T);
          dataWBy[mapWB["χ₂B"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataWByerrp[mapWB["χ₂B"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 1000.);
          dataWByerrm[mapWB["χ₂B"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 1000.);
        }


        abc = CuteHRGHelper::split(chi2Q, '(');
        if (abc.size() >= 2) {
          dataWBx[mapWB["χ₂Q"]].push_back(T);
          dataWBy[mapWB["χ₂Q"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataWByerrp[mapWB["χ₂Q"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 1000.);
          dataWByerrm[mapWB["χ₂Q"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 1000.);
        }


        abc = CuteHRGHelper::split(chi2S, '(');
        if (abc.size() >= 2) {
          dataWBx[mapWB["χ₂S"]].push_back(T);
          dataWBy[mapWB["χ₂S"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataWByerrp[mapWB["χ₂S"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 1000.);
          dataWByerrm[mapWB["χ₂S"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 1000.);
        }

        abc = CuteHRGHelper::split(CBS, '(');
        if (abc.size() >= 2) {
          dataWBx[mapWB["CBS"]].push_back(T);
          dataWBy[mapWB["CBS"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataWByerrp[mapWB["CBS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
          dataWByerrm[mapWB["CBS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        }
      }
    }

    fin.close();
  }

  // WB chi11
  fin.open((string(ThermalFIST_INPUT_FOLDER) + "/lqcd/WB-chi11-1910.14592.dat").c_str());

  if (fin.is_open()) {
    QString tname;

    tname = "χ₁₁BQ";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "χ₁₁QS";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "-χ₁₁BS";
    mapWB[tname] = sizeWB;
    sizeWB++;

    //tname = "CBS";
    //mapWB[tname] = sizeWB;
    //sizeWB++;

    dataWBx.resize(sizeWB);
    dataWBy.resize(sizeWB);
    dataWByerrp.resize(sizeWB);
    dataWByerrm.resize(sizeWB);

    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1 || elems[0].size() == 0)
        continue;

      istringstream iss(elems[0]);

      double T;
      std::string chi11BQ, chi11QS, chi11BS, CBS;

      if (iss >> T >> chi11BQ >> chi11QS >> chi11BS >> CBS) {

        vector<string> abc;

        abc = CuteHRGHelper::split(chi11BQ, '(');
        if (abc.size() >= 2) {
          dataWBx[mapWB["χ₁₁BQ"]].push_back(T);
          dataWBy[mapWB["χ₁₁BQ"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataWByerrp[mapWB["χ₁₁BQ"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
          dataWByerrm[mapWB["χ₁₁BQ"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        }


        abc = CuteHRGHelper::split(chi11QS, '(');
        if (abc.size() >= 2) {
          dataWBx[mapWB["χ₁₁QS"]].push_back(T);
          dataWBy[mapWB["χ₁₁QS"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataWByerrp[mapWB["χ₁₁QS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
          dataWByerrm[mapWB["χ₁₁QS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        }


        abc = CuteHRGHelper::split(chi11BS, '(');
        if (abc.size() >= 2) {
          dataWBx[mapWB["-χ₁₁BS"]].push_back(T);
          dataWBy[mapWB["-χ₁₁BS"]].push_back(-QString::fromStdString(abc[0]).toDouble());
          dataWByerrp[mapWB["-χ₁₁BS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
          dataWByerrm[mapWB["-χ₁₁BS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        }

        //abc = CuteHRGHelper::split(CBS, '(');
        //if (abc.size() >= 2) {
        //  dataWBx[mapWB["CBS"]].push_back(T);
        //  dataWBy[mapWB["CBS"]].push_back(QString::fromStdString(abc[0]).toDouble());
        //  dataWByerrp[mapWB["CBS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        //  dataWByerrm[mapWB["CBS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        //}
      }
    }

    fin.close();
  }

  // HotQCD thermodynamics
  fin.open((string(ThermalFIST_INPUT_FOLDER) + "/lqcd/HotQCD-EoS.dat").c_str());

  if (fin.is_open()) {
    QString tname;

    tname = "p/T⁴";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "e/T⁴";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "(e-3P)/T⁴";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "s/T³";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "CV/T³";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "cs²";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    dataHotQCDx.resize(sizeHotQCD);
    dataHotQCDy.resize(sizeHotQCD);
    dataHotQCDyerrp.resize(sizeHotQCD);
    dataHotQCDyerrm.resize(sizeHotQCD);

    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1 || elems[0].size() == 0)
        continue;

      istringstream iss(elems[0]);

      double T;// , IT4, IT4err, pT4, pT4err, eT4, eT4err, sT3, sT3err, cs2, cs2err;
      std::string IT4, pT4, eT4, sT3, CVT3, cs2;

      if (iss >> T >> IT4 >> pT4 >> eT4 >> sT3 >> CVT3 >> cs2) {

        vector<string> abc;

        abc = CuteHRGHelper::split(pT4, '(');
        if (abc.size() >= 3) {
          dataHotQCDx[mapHotQCD["p/T⁴"]].push_back(T);
          dataHotQCDy[mapHotQCD["p/T⁴"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["p/T⁴"]].push_back(QString::fromStdString(abc[2].substr(1, abc[2].size() - 2)).toDouble() / 1000.);
          dataHotQCDyerrm[mapHotQCD["p/T⁴"]].push_back(QString::fromStdString(abc[1].substr(1, abc[1].size() - 2)).toDouble() / 1000.);
        }


        abc = CuteHRGHelper::split(eT4, '(');
        if (abc.size() >= 3) {
          dataHotQCDx[mapHotQCD["e/T⁴"]].push_back(T);
          dataHotQCDy[mapHotQCD["e/T⁴"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["e/T⁴"]].push_back(QString::fromStdString(abc[2].substr(1, abc[2].size() - 2)).toDouble() / 100.);
          dataHotQCDyerrm[mapHotQCD["e/T⁴"]].push_back(QString::fromStdString(abc[1].substr(1, abc[1].size() - 2)).toDouble() / 100.);
        }


        abc = CuteHRGHelper::split(IT4, '(');
        if (abc.size() >= 3) {
          dataHotQCDx[mapHotQCD["(e-3P)/T⁴"]].push_back(T);
          dataHotQCDy[mapHotQCD["(e-3P)/T⁴"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["(e-3P)/T⁴"]].push_back(QString::fromStdString(abc[2].substr(1, abc[2].size() - 2)).toDouble() / 100.);
          dataHotQCDyerrm[mapHotQCD["(e-3P)/T⁴"]].push_back(QString::fromStdString(abc[1].substr(1, abc[1].size() - 2)).toDouble() / 100.);
        }

        abc = CuteHRGHelper::split(sT3, '(');
        if (abc.size() >= 3) {
          dataHotQCDx[mapHotQCD["s/T³"]].push_back(T);
          dataHotQCDy[mapHotQCD["s/T³"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["s/T³"]].push_back(QString::fromStdString(abc[2].substr(1, abc[2].size() - 2)).toDouble() / 100.);
          dataHotQCDyerrm[mapHotQCD["s/T³"]].push_back(QString::fromStdString(abc[1].substr(1, abc[1].size() - 2)).toDouble() / 100.);
        }

        abc = CuteHRGHelper::split(CVT3, '(');
        if (abc.size() >= 3) {
          dataHotQCDx[mapHotQCD["CV/T³"]].push_back(T);
          dataHotQCDy[mapHotQCD["CV/T³"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["CV/T³"]].push_back(QString::fromStdString(abc[2].substr(1, abc[2].size() - 2)).toDouble());
          dataHotQCDyerrm[mapHotQCD["CV/T³"]].push_back(QString::fromStdString(abc[1].substr(1, abc[1].size() - 2)).toDouble());
        }

        abc = CuteHRGHelper::split(cs2, '(');
        if (abc.size() >= 3) {
          dataHotQCDx[mapHotQCD["cs²"]].push_back(T);
          dataHotQCDy[mapHotQCD["cs²"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["cs²"]].push_back(QString::fromStdString(abc[2].substr(1, abc[2].size() - 2)).toDouble() / 1000.);
          dataHotQCDyerrm[mapHotQCD["cs²"]].push_back(QString::fromStdString(abc[1].substr(1, abc[1].size() - 2)).toDouble() / 1000.);
        }
      }
    }

    fin.close();
  }

  // HotQCD chi2
  fin.open((string(ThermalFIST_INPUT_FOLDER) + "/lqcd/HotQCD-chi2-1203.0784.dat").c_str());

  if (fin.is_open()) {
    QString tname;

    //tname = "χ₂B";
    //mapHotQCD[tname] = sizeHotQCD;
    //sizeHotQCD++;

    tname = "χ₂Q";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "χ₂S";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "χ₁₁BQ";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "-χ₁₁BS";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "χ₁₁QS";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    dataHotQCDx.resize(sizeHotQCD);
    dataHotQCDy.resize(sizeHotQCD);
    dataHotQCDyerrp.resize(sizeHotQCD);
    dataHotQCDyerrm.resize(sizeHotQCD);

    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1 || elems[0].size() == 0)
        continue;

      istringstream iss(elems[0]);

      double T;
      std::string chi2B, chi2Q, chi2S, chi11BQ, chi11BSm, chi11QS;

      if (iss >> T >> chi2B >> chi2Q >> chi2S >> chi11BSm >> chi11BQ >> chi11QS) {

        vector<string> abc;

        abc = CuteHRGHelper::split(chi2B, '(');
        //if (abc.size() >= 2) {
        //  dataHotQCDx[mapHotQCD["χ₂B"]].push_back(T);
        //  dataHotQCDy[mapHotQCD["χ₂B"]].push_back(QString::fromStdString(abc[0]).toDouble());
        //  dataHotQCDyerrp[mapHotQCD["χ₂B"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        //  dataHotQCDyerrm[mapHotQCD["χ₂B"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        //}


        abc = CuteHRGHelper::split(chi2Q, '(');
        if (abc.size() >= 2) {
          dataHotQCDx[mapHotQCD["χ₂Q"]].push_back(T);
          dataHotQCDy[mapHotQCD["χ₂Q"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["χ₂Q"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
          dataHotQCDyerrm[mapHotQCD["χ₂Q"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        }


        abc = CuteHRGHelper::split(chi2S, '(');
        if (abc.size() >= 2) {
          dataHotQCDx[mapHotQCD["χ₂S"]].push_back(T);
          dataHotQCDy[mapHotQCD["χ₂S"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["χ₂S"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
          dataHotQCDyerrm[mapHotQCD["χ₂S"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        }

        abc = CuteHRGHelper::split(chi11BQ, '(');
        if (abc.size() >= 2) {
          dataHotQCDx[mapHotQCD["χ₁₁BQ"]].push_back(T);
          dataHotQCDy[mapHotQCD["χ₁₁BQ"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["χ₁₁BQ"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
          dataHotQCDyerrm[mapHotQCD["χ₁₁BQ"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        }

        abc = CuteHRGHelper::split(chi11BSm, '(');
        if (abc.size() >= 2) {
          dataHotQCDx[mapHotQCD["-χ₁₁BS"]].push_back(T);
          dataHotQCDy[mapHotQCD["-χ₁₁BS"]].push_back(-QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["-χ₁₁BS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
          dataHotQCDyerrm[mapHotQCD["-χ₁₁BS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        }

        abc = CuteHRGHelper::split(chi11QS, '(');
        if (abc.size() >= 2) {
          dataHotQCDx[mapHotQCD["χ₁₁QS"]].push_back(T);
          dataHotQCDy[mapHotQCD["χ₁₁QS"]].push_back(QString::fromStdString(abc[0]).toDouble());
          dataHotQCDyerrp[mapHotQCD["χ₁₁QS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
          dataHotQCDyerrm[mapHotQCD["χ₁₁QS"]].push_back(QString::fromStdString(abc[1].substr(0, abc[1].size() - 1)).toDouble() / 10000.);
        }
      }
    }

    fin.close();
  }

  // HotQCD chi2-chi8
  fin.open((string(ThermalFIST_INPUT_FOLDER) + "/lqcd/HotQCD-chi2-chi8-2001.08530.dat").c_str());

  if (fin.is_open()) {
    QString tname;

    tname = "χ₂B";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "χ₄B";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "χ₆B";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    tname = "χ₈B";
    mapHotQCD[tname] = sizeHotQCD;
    sizeHotQCD++;

    dataHotQCDx.resize(sizeHotQCD);
    dataHotQCDy.resize(sizeHotQCD);
    dataHotQCDyerrp.resize(sizeHotQCD);
    dataHotQCDyerrm.resize(sizeHotQCD);

    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1 || elems[0].size() == 0)
        continue;

      istringstream iss(elems[0]);

      double T;
      //std::string chi2B, chi2Q, chi2S, chi11BQ, chi11BSm, chi11QS;
      double chi2B, chi2Berr, chi4B, chi4Berr, chi6B, chi6Berr, chi8B, chi8Berr;

      if (iss >> T >> chi2B >> chi2Berr >> chi4B >> chi4Berr >> chi6B >> chi6Berr >> chi8B >> chi8Berr) {

        vector<string> abc;

        dataHotQCDx[mapHotQCD["χ₂B"]].push_back(T);
        dataHotQCDy[mapHotQCD["χ₂B"]].push_back(chi2B);
        dataHotQCDyerrp[mapHotQCD["χ₂B"]].push_back(chi2Berr);
        dataHotQCDyerrm[mapHotQCD["χ₂B"]].push_back(chi2Berr);


        dataHotQCDx[mapHotQCD["χ₄B"]].push_back(T);
        dataHotQCDy[mapHotQCD["χ₄B"]].push_back(chi4B);
        dataHotQCDyerrp[mapHotQCD["χ₄B"]].push_back(chi4Berr);
        dataHotQCDyerrm[mapHotQCD["χ₄B"]].push_back(chi4Berr);


        dataHotQCDx[mapHotQCD["χ₆B"]].push_back(T);
        dataHotQCDy[mapHotQCD["χ₆B"]].push_back(chi6B);
        dataHotQCDyerrp[mapHotQCD["χ₆B"]].push_back(chi6Berr);
        dataHotQCDyerrm[mapHotQCD["χ₆B"]].push_back(chi6Berr);

        dataHotQCDx[mapHotQCD["χ₈B"]].push_back(T);
        dataHotQCDy[mapHotQCD["χ₈B"]].push_back(chi8B);
        dataHotQCDyerrp[mapHotQCD["χ₈B"]].push_back(chi8Berr);
        dataHotQCDyerrm[mapHotQCD["χ₈B"]].push_back(chi8Berr);
      }
    }

    fin.close();
  }

  // WB chi2-chi8
  fin.open((string(ThermalFIST_INPUT_FOLDER) + "/lqcd/WB-chiB-1805.04445.dat").c_str());

  if (fin.is_open()) {
    QString tname;

    tname = "χ₂B";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "χ₄B";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "χ₆B";
    mapWB[tname] = sizeWB;
    sizeWB++;

    tname = "χ₈B";
    mapWB[tname] = sizeWB;
    sizeWB++;

    dataWBx.resize(sizeWB);
    dataWBy.resize(sizeWB);
    dataWByerrp.resize(sizeWB);
    dataWByerrm.resize(sizeWB);

    char cc[2000];
    while (!fin.eof()) {
      fin.getline(cc, 2000);
      string tmp = string(cc);
      vector<string> elems = CuteHRGHelper::split(tmp, '#');
      if (elems.size() < 1 || elems[0].size() == 0)
        continue;

      istringstream iss(elems[0]);

      double T;
      //std::string chi2B, chi2Q, chi2S, chi11BQ, chi11BSm, chi11QS;
      double chi2B, chi2Berr, chi4B, chi4Berr, chi6B, chi6Berr, chi8B, chi8Berr;

      if (iss >> T >> chi2B >> chi2Berr >> chi4B >> chi4Berr >> chi6B >> chi6Berr >> chi8B >> chi8Berr) {

        vector<string> abc;

        dataWBx[mapWB["χ₂B"]].push_back(T);
        dataWBy[mapWB["χ₂B"]].push_back(chi2B);
        dataWByerrp[mapWB["χ₂B"]].push_back(chi2Berr);
        dataWByerrm[mapWB["χ₂B"]].push_back(chi2Berr);


        dataWBx[mapWB["χ₄B"]].push_back(T);
        dataWBy[mapWB["χ₄B"]].push_back(chi4B);
        dataWByerrp[mapWB["χ₄B"]].push_back(chi4Berr);
        dataWByerrm[mapWB["χ₄B"]].push_back(chi4Berr);


        dataWBx[mapWB["χ₆B"]].push_back(T);
        dataWBy[mapWB["χ₆B"]].push_back(chi6B);
        dataWByerrp[mapWB["χ₆B"]].push_back(chi6Berr);
        dataWByerrm[mapWB["χ₆B"]].push_back(chi6Berr);

        dataWBx[mapWB["χ₈B"]].push_back(T);
        dataWBy[mapWB["χ₈B"]].push_back(chi8B);
        dataWByerrp[mapWB["χ₈B"]].push_back(chi8Berr);
        dataWByerrm[mapWB["χ₈B"]].push_back(chi8Berr);
      }
    }

    fin.close();
  }
}

QString EquationOfStateTab::getParameterName() const {
  if (comboMode->currentIndex() == 2)
    return "μB[MeV]";
  else if (comboMode->currentIndex() == 3)
    return "eB[GeV^2]";
  else
    return "T[MeV]";
}

void EquationOfStateTab::replot() {
    bool flipAxes = false;
    flipAxes = CBflipAxes->isChecked();

    QCPAxis* xAxis = plotDependence->xAxis;
    QCPAxis* yAxis = plotDependence->yAxis;

    if (flipAxes) {
      xAxis = plotDependence->yAxis;
      yAxis = plotDependence->xAxis;
    }

    xAxis->setLabel(getParameterName());

//    if (comboMode->currentIndex() == 2) {
//      xAxis->setLabel("μB (MeV)");
//    }
//    else {
//      xAxis->setLabel("T (MeV)");
//    }

    int index = comboQuantity->currentIndex();

//    plotDependence->clearGraphs();
    plotDependence->clearPlottables();

    xAxis->setLabelFont(QFont("Arial", font().pointSize() + 4));
    yAxis->setLabelFont(QFont("Arial", font().pointSize() + 4));

    if (index >= 0 && index < paramnames.size()) {
      std::vector<double> yvalues;
      std::vector<double> xvalues = varvalues;

      if (CBratio->isChecked() && !CBxAxis->isChecked()) {
        int index2 = comboQuantity2->currentIndex();
        if (!(index2 >= 0 && index2 < paramnames.size()))
          return;
        yAxis->setLabel(paramnames[index] + "/" + paramnames[index2]);
        if ((paramnames[index] == "ni/T³" && paramnames[index2] == "ni/T³") ||
          (paramnames[index] == "ni[fm^-3]" && paramnames[index2] == "ni[fm^-3]")) {
          QString tname = paramnames[index] + "/" + paramnames[index2];
          int pid = comboParticle->currentIndex(), pid2 = comboParticle2->currentIndex();
          if (pid >= 0 && pid < model->TPS()->Particles().size()
            && pid2 >= 0 && pid2 < model->TPS()->Particles().size())
            tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) + ")/n(" 
                  + QString::fromStdString(model->TPS()->Particles()[pid2].Name()) + ")";
          yAxis->setLabel(tname);
        }
        yvalues = getValuesRatio(index, index2);
      }
      else {
        yAxis->setLabel(paramnames[index]);
        if (paramnames[index] == "ni/T³") {
          QString tname = paramnames[index];
          int pid = comboParticle->currentIndex();
          if (pid >= 0 && pid < model->TPS()->Particles().size())
            tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) +")/T³";
          yAxis->setLabel(tname);
        }
        if (paramnames[index] == "ni[fm^-3]") {
          QString tname = paramnames[index];
          int pid = comboParticle->currentIndex();
          if (pid >= 0 && pid < model->TPS()->Particles().size())
            tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) + ")[fm^-3]³";
          yAxis->setLabel(tname);
        }
        yvalues = getValues(index);

        if (CBxAxis->isChecked()) {
          xvalues = getValues(comboQuantity2->currentIndex(), 1);

          int index2 = comboQuantity2->currentIndex();
          xAxis->setLabel(paramnames[index2]);
          if (paramnames[index2] == "ni/T³") {
            QString tname = paramnames[index2];
            int pid = comboParticle2->currentIndex();
            if (pid >= 0 && pid < model->TPS()->Particles().size())
              tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) +")/T³";
            xAxis->setLabel(tname);
          }
          if (paramnames[index2] == "ni[fm^-3]") {
            QString tname = paramnames[index2];
            int pid = comboParticle2->currentIndex();
            if (pid >= 0 && pid < model->TPS()->Particles().size())
              tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) + ")[fm^-3]³";
            xAxis->setLabel(tname);
          }
        }
      }
      double tmin = 0., tmax = 0.;
      for (int i = 0; i<yvalues.size(); ++i) {
        tmin = std::min(tmin, yvalues[i]);
        tmax = std::max(tmax, yvalues[i]);
      }

      if (tmin > 0.)
        tmin *= 0.;

      if (yvalues.size() > 0)
        yAxis->setRange(tmin, tmax);

      if (CBxAxis->isChecked()) {
        double xmin = 1.e12, xmax = -1.e12;
        for (int i = 0; i<xvalues.size(); ++i) {
          xmin = std::min(xmin, xvalues[i]);
          xmax = std::max(xmax, xvalues[i]);
        }

//        if (xmin > 0.)
//          xmin *= 0.;

        if (xvalues.size() > 0)
          xAxis->setRange(xmin, xmax);
      }

      // Lattice data for muB = 0 only
      if (comboMode->currentIndex() < 2 && spinmuB->value() == 0.0 && !model->ConstrainMuB() && !flipAxes && !CBxAxis->isChecked())
        plotLatticeData();

      int graphNumber = plotDependence->graphCount();

//      plotDependence->addGraph();
//      plotDependence->graph(graphNumber)->setName("Model");
//      plotDependence->graph(graphNumber)->setPen(QPen(Qt::black, 2, Qt::SolidLine));
//
//      plotDependence->graph(graphNumber)->data()->clear();
      QCPCurve *curve = new QCPCurve(plotDependence->xAxis, plotDependence->yAxis);
      curve->setName("Model");
      curve->setPen(QPen(Qt::black, 2, Qt::SolidLine));
      curve->data()->clear();

      for(int i=0;i<yvalues.size();++i) {
        if (flipAxes)
          curve->addData(yvalues[i], xvalues[i]);
        else
          curve->addData(xvalues[i], yvalues[i]);
      }

      tmin = std::min(tmin, plotDependence->yAxis->range().lower);
      tmax = std::max(tmax, plotDependence->yAxis->range().upper);

      double pmin = spinTMin->value(), pmax = spinTMax->value();
      if (comboMode->currentIndex() == 3) {
        pmin = spinEMin->value();
        pmax = spinEMax->value();
      }

      if (!CBxAxis->isChecked())
        xAxis->setRange(pmin, pmax);
      if (yvalues.size() > 0)
        yAxis->setRange(1.1*tmin, 1.1*tmax);

      plotDependence->legend->setFont(QFont("Arial", font().pointSize() + 2));
      plotDependence->legend->setVisible(true);

      plotDependence->replot();
    }
}

void EquationOfStateTab::recomputeCalcTable() {
  calcTable.clear();

  calcTable.parameter_name = getParameterName();
  for(int i = 0; i < varvalues.size(); ++i) {
    calcTable.parameter_values.push_back(varvalues[i]);
    calcTable.temperature_values.push_back(paramsTD[i].T);
  }

  for(int ir = 0; ir < paramnames.size(); ++ir) {
    if (paramnames[ir] == "ni/T³" || paramnames[ir] == "ni[fm^-3]")
      continue;

    calcTable.quantities_names.push_back(paramnames[ir]);
    calcTable.quantities_values.push_back(getValues(ir));
  }

  for(int ic = 0; ic < model->TPS()->Particles().size(); ++ic) {
    calcTable.densities_names.push_back(QString::fromStdString(model->TPS()->Particles()[ic].Name()));
  }

  for(int ir = 0; ir < varvalues.size(); ++ir) {
    calcTable.densities_values.push_back(paramsTD[ir].densities);
  }
}
