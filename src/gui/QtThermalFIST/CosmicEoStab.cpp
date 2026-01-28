/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "CosmicEoStab.h"

#include <algorithm>
#include <fstream>
#include <sstream>

#include <QLayout>
#include <QGroupBox>
#include <QLabel>
#include <QTextStream>
#include <QBuffer>

#include "ThermalFISTConfig.h"
#include "WasmFileIO.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGEV/ThermalModelEVCrossterms.h"
#include "HRGVDW/ThermalModelVDW.h"
#include "HRGBase/xMath.h"
#include "HRGRealGas/ThermalModelRealGas.h"

using namespace std;
using namespace thermalfist;

void CosmicEoSWorker::run() {
  double Tprev = 0., Tprev2 = 0.;
  vector<double> chemsprev, chemsprev2;

  int i = 0;
  vector<double> params;
  if (reverseDir) {
    for (double param = Tmax; param >= Tmin - 1.e-10; param -= dT) {
      params.push_back(param);
    }
  }
  else {
    for (double param = Tmin; param <= Tmax + 1.e-10; param += dT) {
      params.push_back(param);
    }
  }
  for (double param : params) {
    if (i < paramsTD->size()) {
      if (i > 1) {
        double dTcur = param - Tprev2;
        double dTder = Tprev - Tprev2;
        for(int i = 0; i < 5; i++) {
          chems[i] = chemsprev2[i] + (chemsprev[i] - chemsprev2[i]) * dTcur / dTder;
        }
      }

      chems = cosmos->SolveChemicalPotentials(param * 1.e-3, chems);
      // printf("Param: %f %E %E %E %E %E\n", param, chems[0], chems[1], chems[2], chems[3], chems[4]);
      Tprev2 = Tprev;
      Tprev = param;
      chemsprev2 = chemsprev;
      chemsprev = chems;

      varvalues->operator [](i) = param;
      double T = param;
      ThermodynamicsCosmic tTD;

      tTD.T      = cosmos->Temperature();
      tTD.muB    = cosmos->ChemicalPotentials()[0];
      tTD.muQ    = cosmos->ChemicalPotentials()[1];
      tTD.mu_e   = cosmos->ChemicalPotentials()[2];
      tTD.mu_mu  = cosmos->ChemicalPotentials()[3];
      tTD.mu_tau = cosmos->ChemicalPotentials()[4];

      tTD.p = cosmos->Pressure();
      tTD.e = cosmos->EnergyDensity();
      tTD.s = cosmos->EntropyDensity();
      tTD.I = tTD.e - 3. * tTD.p;
      tTD.pT4 = cosmos->Pressure() / pow(T * 1.e-3, 4) / thermalfist::xMath::GeVtoifm3();
      tTD.eT4 = cosmos->EnergyDensity() / pow(T * 1.e-3, 4) / thermalfist::xMath::GeVtoifm3();
      tTD.sT3 = cosmos->EntropyDensity() / pow(T * 1.e-3, 3) / thermalfist::xMath::GeVtoifm3();
      tTD.IT4 = tTD.eT4 - 3. * tTD.pT4;
      tTD.rhoB = cosmos->BaryonDensity();
      tTD.rhoQ = cosmos->ElectricChargeDensity();
      tTD.rhoE = cosmos->LeptonFlavorDensity(LeptonFlavor::Electron);
      tTD.rhoMu = cosmos->LeptonFlavorDensity(LeptonFlavor::Muon);
      tTD.rhoTau = cosmos->LeptonFlavorDensity(LeptonFlavor::Tau);
      std::vector<double> NonQCDDensities(13, 0.);
      for(int i = 0; i < 13; ++i) {
        NonQCDDensities[i] = cosmos->GetDensity(i);
      }
      tTD.densities = NonQCDDensities;
      tTD.flag = true;
      paramsTD->operator [](i) = tTD;

      ThermalModelBase *model = cosmos->HRGModel();

      model->CalculateFeeddown();

      Thermodynamics tTDHRG;
      tTDHRG.p = model->Pressure();
      tTDHRG.e = model->EnergyDensity();
      tTDHRG.s = model->EntropyDensity();
      tTDHRG.I = tTD.e - 3. * tTD.p;
      tTDHRG.pT4 = model->Pressure() / pow(T * 1.e-3, 4) / thermalfist::xMath::GeVtoifm3();
      tTDHRG.eT4 = model->EnergyDensity() / pow(T * 1.e-3, 4) / thermalfist::xMath::GeVtoifm3();
      tTDHRG.sT3 = model->EntropyDensity() / pow(T * 1.e-3, 3) / thermalfist::xMath::GeVtoifm3();
      tTDHRG.IT4 = tTDHRG.eT4 - 3. * tTDHRG.pT4;
      tTDHRG.nh  = model->HadronDensity();
      tTDHRG.T = model->Parameters().T;
      tTDHRG.muB = model->Parameters().muB;
      tTDHRG.muQ = model->Parameters().muQ;
      tTDHRG.muS = model->Parameters().muS;
      tTDHRG.muC = model->Parameters().muC;
      tTDHRG.densities = model->AllDensities();
      tTDHRG.rhoB = model->BaryonDensity();
      tTDHRG.rhoQ = model->ElectricChargeDensity();
      tTDHRG.rhoS = model->StrangenessDensity();
      tTDHRG.rhoC = model->CharmDensity();
      tTDHRG.flag = true;
      paramsTDHRG->operator [](i) = tTDHRG;

      (*currentSize)++;
    }
    i++;
  }
  if (emitSignal) {
    emit calculated();
  }
}

CosmicEoSTab::CosmicEoSTab(QWidget *parent, ThermalModelBase *modelop) :
    QWidget(parent)
{
    TPS = modelop->TPS();
    model = new ThermalModelIdeal(modelop->TPS());
    *model = *modelop;
    cosmos = new CosmicEoS(model, false);

    int index = 0;
    QString tname;

    tname = "μB[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "μQ[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "μ_e[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "μ_μ[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "μ_τ[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "T[MeV]";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

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

    tname = "p/T⁴_hrg";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "e/T⁴_hrg";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "s/T³_hrg";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "(e-3P)/T⁴_hrg";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "Δ_hrg";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρB/s";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρQ/s";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρ_e/s";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρ_μ/s";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "ρ_τ/s";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "n_i^lep/T³";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "n_i^had/T³";
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
    CBflipAxes->setChecked(true);
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
    connect(comboQuantity, SIGNAL(currentIndexChanged(int)), this, SLOT(fillParticleLists()));
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
    connect(comboQuantity2, SIGNAL(currentIndexChanged(int)), this, SLOT(fillParticleLists()));
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
    plotDependence->xAxis->setLabel("T[MeV]");
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
    layLeft->addWidget(CBflipAxes, 0, Qt::AlignLeft);
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

//    QHBoxLayout* layMode = new QHBoxLayout();
//    layMode->setAlignment(Qt::AlignLeft);
//    comboMode = new QComboBox();
//    comboMode->addItem(tr("Const. μB"));
//    comboMode->addItem(tr("Const. μB/T"));
//    comboMode->addItem(tr("Const. T"));
//    comboMode->setCurrentIndex(0);
//    connect(comboMode, SIGNAL(currentIndexChanged(int)), this, SLOT(modelChanged()));
//    layMode->addWidget(new QLabel(tr("Mode:")));
//    layMode->addWidget(comboMode);

    QGroupBox *grCons = new QGroupBox(tr("Conservation laws:"));
    QGridLayout *layCons = new QGridLayout();
    layCons->setAlignment(Qt::AlignLeft);

    QLabel *labelnB = new QLabel(tr("n<sub>B</sub>/s (x10<sup>-11</sup>):"));
    spinnB = new QDoubleSpinBox();
    spinnB->setDecimals(5);
    spinnB->setMinimum(-1.e12);
    spinnB->setMaximum(1.e12);
    spinnB->setValue(8.6);
    spinnB->setSingleStep(0.1);

    QLabel *labelnQ = new QLabel(tr("n<sub>Q</sub>/s (x10<sup>-11</sup>):"));
    spinnQ = new QDoubleSpinBox();
    spinnQ->setDecimals(5);
    spinnQ->setMinimum(-10000.);
    spinnQ->setMaximum(10000.);
    spinnQ->setValue(0.);
    spinnQ->setSingleStep(0.1);

    QLabel *labelnE = new QLabel(tr("ℓ<sub>e</sub>/s:"));
    spinnE = new QDoubleSpinBox();
    spinnE->setDecimals(12);
    spinnE->setMinimum(-10000.);
    spinnE->setMaximum(10000.);
    spinnE->setValue(0.);
    spinnE->setSingleStep(0.1);

    QLabel *labelnMu = new QLabel(tr("ℓ<sub>μ</sub>/s:"));
    spinnMu = new QDoubleSpinBox();
    spinnMu->setDecimals(12);
    spinnMu->setMinimum(-10000.);
    spinnMu->setMaximum(10000.);
    spinnMu->setValue(0.);
    spinnMu->setSingleStep(0.1);

    QLabel *labelnTau = new QLabel(tr("ℓ<sub>τ</sub>/s:"));
    spinnTau = new QDoubleSpinBox();
    spinnTau->setDecimals(12);
    spinnTau->setMinimum(-10000.);
    spinnTau->setMaximum(10000.);
    spinnTau->setValue(0.);
    spinnTau->setSingleStep(0.1);

    layCons->addWidget(labelnB, 0, 0);
    layCons->addWidget(spinnB, 0, 1);
    layCons->addWidget(labelnQ, 0, 2);
    layCons->addWidget(spinnQ, 0, 3);
    layCons->addWidget(labelnE, 1, 0);
    layCons->addWidget(spinnE, 1, 1);
    layCons->addWidget(labelnMu, 1, 2);
    layCons->addWidget(spinnMu, 1, 3);
    layCons->addWidget(labelnTau, 1, 4);
    layCons->addWidget(spinnTau, 1, 5);

    grCons->setLayout(layCons);

    QGroupBox *grRange = new QGroupBox(tr("Temperature range:"));

    QVBoxLayout *layParamBounds = new QVBoxLayout();
    QGridLayout *layBounds = new QGridLayout();
    layBounds->setAlignment(Qt::AlignLeft);

    labelTMin = new QLabel(tr("T<sub>min</sub> (MeV):"));
    spinTMin = new QDoubleSpinBox();
    spinTMin->setDecimals(3);
    spinTMin->setMinimum(0.01);
    spinTMin->setMaximum(10000.);
    spinTMin->setValue(10.);

    labelTMax = new QLabel(tr("T<sub>max</sub> (MeV):"));
    spinTMax = new QDoubleSpinBox();
    spinTMax->setDecimals(3);
    spinTMax->setMinimum(0.01);
    spinTMax->setMaximum(10000.);
    spinTMax->setValue(180.);

    labeldT = new QLabel(tr("∆T (MeV):"));
    spindT = new QDoubleSpinBox();
    spindT->setDecimals(3);
    spindT->setMinimum(0.);
    spindT->setMaximum(10000.);
    spindT->setValue(1.);

    CBreverseDir = new QCheckBox(tr("Reverse direction"));
    CBreverseDir->setChecked(true);


    layBounds->addWidget(labelTMin, 0, 0);
    layBounds->addWidget(spinTMin, 0, 1);
    layBounds->addWidget(labelTMax, 0, 2);
    layBounds->addWidget(spinTMax, 0, 3);
    layBounds->addWidget(labeldT, 0, 4);
    layBounds->addWidget(spindT, 0, 5);
    // layBounds->addWidget(CBreverseDir, 0, 6);

//    labelConstr = new QLabel(tr(""));


    layParamBounds->addLayout(layBounds);
    layParamBounds->addWidget(CBreverseDir, 0, Qt::AlignLeft);
    // layParamBounds->addWidget(labelConstr, 0, Qt::AlignLeft);

    grRange->setLayout(layParamBounds);

    buttonCalculate = new QPushButton(tr("Calculate"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));

    layRight->addWidget(grModelConfig);
    layRight->addWidget(grCons);
    layRight->addWidget(grRange);
    layRight->addWidget(buttonCalculate, 0, Qt::AlignLeft);

    mainLayout->addLayout(layLeft, 1);
    mainLayout->addLayout(layRight);

    setLayout(mainLayout);

    fCurrentSize = 0;
    fStop = 1;
    fRunning = false;
    fExpectedIterations = 0;

    calcTimer = new QTimer(this);
    connect(calcTimer, SIGNAL(timeout()), this, SLOT(checkProgress()));

    fillParticleLists();

    modelChanged();
    
    replot();

    cpath = QApplication::applicationDirPath();
}

CosmicEoSTab::~CosmicEoSTab() {
    if (model!=NULL) delete model;
}

ThermalModelConfig CosmicEoSTab::getConfig()
{
  ThermalModelConfig ret = configWidget->updatedConfig();

  ret.T = spinTMin->value() * 1.e-3;
  ret.muB = 0.;
  ret.muQ = 0.;
  ret.muS = 0.;
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

void CosmicEoSTab::calculate() {
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

      if (cosmos != NULL)
        delete cosmos;
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



      cosmos = new CosmicEoS(model, false);
      if (config.UseEMMPions)
        cosmos->SetPionsInteracting(true, config.EMMPionFPi);
      if (config.UseEMMKaons)
        cosmos->SetKaonsInteracting(true, config.EMMKaonFKa);

      paramsTD.resize(0);
      paramsTDHRG.resize(0);
      varvalues.resize(0);

      int iters = static_cast<int>((spinTMax->value() - spinTMin->value()) / spindT->value()) + 1;

      fCurrentSize = 0;
      fExpectedIterations = iters;  // Store for completion detection

      paramsTD.resize(iters);
      paramsTDHRG.resize(iters);
      varvalues.resize(iters);

      fStop = 0;

      vector<double> asymmetries = {
              spinnB->value() * 1.e-11,
              spinnQ->value() * 1.e-11,
              spinnE->value(),
              spinnMu->value(),
              spinnTau->value()
      };

      bool useThreading = WasmFileIO::isThreadingAvailable();
      // Don't emit signal from worker thread in WASM - causes memory errors
      bool emitSignalFromWorker = !useThreading;

      CosmicEoSWorker *wrk = new CosmicEoSWorker(cosmos, spinTMin->value(), spinTMax->value(), spindT->value(),
                                                 asymmetries, &paramsTD, &paramsTDHRG,
                                                 &varvalues, &fCurrentSize, &fStop, CBreverseDir->isChecked(),
                                                 emitSignalFromWorker, this);

      connect(wrk, SIGNAL(finished()), wrk, SLOT(deleteLater()));

      buttonCalculate->setText(tr("Stop"));
      fRunning = true;

      if (useThreading) {
        // Threading available - run in background thread
        // Don't connect calculated() signal - use timer-based completion detection instead
        // to avoid WASM threading issues with cross-thread signal emission
        wrk->start();
        calcTimer->start(10);
      } else {
        // No threading (single-threaded WASM) - run synchronously
        // Signal emission is safe here since everything is on main thread
        connect(wrk, SIGNAL(calculated()), this, SLOT(finalize()));
        wrk->run();
        finalize();
        wrk->deleteLater();
      }
    }
    else {
      fStop = 1;
    }
}

void CosmicEoSTab::finalize() {
    calcTimer->stop();
    buttonCalculate->setText(tr("Calculate"));
    fRunning = false;
    replot();
}

void CosmicEoSTab::checkProgress() {
    // Update the plot with current progress
    replot();

    // Check if calculation is complete or stopped (timer-based completion detection for WASM threading)
    // This avoids cross-thread signal emission which can cause memory errors in WASM
    if (fRunning && (fCurrentSize >= fExpectedIterations || fStop)) {
        finalize();
    }
}

void CosmicEoSTab::modelChanged()
{
  if (CBratio->isChecked() || CBxAxis->isChecked()) {
    comboQuantity2->setEnabled(true);
  }
  else {
    comboQuantity2->setEnabled(false);
  }

  if (comboQuantity->currentText() == "n_i^lep/T³" || comboQuantity->currentText() == "n_i^had/T³") {
    labelParticle->setVisible(true);
    comboParticle->setVisible(true);
    labelFeeddown->setVisible(comboQuantity->currentText() == "n_i^had/T³");
    comboFeeddown->setVisible(comboQuantity->currentText() == "n_i^had/T³");
  }
  else {
    labelParticle->setVisible(false);
    comboParticle->setVisible(false);
    labelFeeddown->setVisible(false);
    comboFeeddown->setVisible(false);
  }

  if ((comboQuantity2->currentText() == "n_i^lep/T³" || comboQuantity2->currentText() == "n_i^had/T³")
    && (CBratio->isChecked() || CBxAxis->isChecked())) {
    labelParticle2->setVisible(true);
    comboParticle2->setVisible(true);
    labelFeeddown2->setVisible(comboQuantity2->currentText() == "n_i^had/T³");
    comboFeeddown2->setVisible(comboQuantity2->currentText() == "n_i^had/T³");
  }
  else {
    labelParticle2->setVisible(false);
    comboParticle2->setVisible(false);
    labelFeeddown2->setVisible(false);
    comboFeeddown2->setVisible(false);
  }

  configWidget->buttonConservationLaws->setEnabled(false);
  configWidget->comboEnsemble->setEnabled(false);
}

void CosmicEoSTab::resetTPS()
{
  model->ChangeTPS(model->TPS());
  fillParticleLists();

  configWidget->setModel(model);
  configWidget->currentConfig.vdWparams = QvdWParameters::GetParameters(model->TPS(), &configWidget->currentConfig);
}

void CosmicEoSTab::fillParticleLists()
{
  for (int ii = 0; ii < 2; ++ii) {
    QComboBox *combo;
    if (ii == 0)
      combo = comboParticle;
    else
      combo = comboParticle2;

    int tind = combo->currentIndex();
    combo->clear();
    if (comboQuantity->currentText() == "n_i^had/T³") {
      for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
        const ThermalParticle &part = model->TPS()->Particles()[i];
        combo->addItem(QString::number(part.PdgId()) + " - " + QString::fromStdString(part.Name()));
      }
    }
    else {
      for (int i = 0; i < cosmos->NumberOfElectroWeakSpecies(); ++i) {
        combo->addItem(QString::fromStdString(cosmos->GetSpeciesName(i)));
      }
    }
    if (tind >= 0 && tind < combo->count())
      combo->setCurrentIndex(tind);
    else
      combo->setCurrentIndex(0);
  }
}

void CosmicEoSTab::contextMenuRequest(QPoint pos)
{
  QMenu *menu = new QMenu(this);
  menu->setAttribute(Qt::WA_DeleteOnClose);

  menu->addAction("Save as pdf...", this, SLOT(saveAsPdf()));
  menu->addAction("Save as png...", this, SLOT(saveAsPng()));
  menu->addAction("Write computed values to file...", this, SLOT(saveAsAscii()));

  menu->popup(plotDependence->mapToGlobal(pos));
}

void CosmicEoSTab::saveAsPdf()
{
  saveAs(0);
}

void CosmicEoSTab::saveAsPng()
{
  saveAs(1);
}

void CosmicEoSTab::saveAsAscii()
{
  saveAs(2);
}

void CosmicEoSTab::saveAs(int type)
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

#ifdef Q_OS_WASM
  // WASM: Save to buffer and trigger browser download
  QString fileName = tname + "." + exts[type];

  if (type == 0) {
    QString tempPath = WasmFileIO::getSandboxTempDir() + "/" + fileName;
    plotDependence->savePdf(tempPath, plotDependence->width(), plotDependence->height());
    QFile file(tempPath);
    if (file.open(QIODevice::ReadOnly)) {
      WasmFileIO::saveFile(file.readAll(), fileName);
      file.close();
    }
  }
  else if (type == 1) {
    QString tempPath = WasmFileIO::getSandboxTempDir() + "/" + fileName;
    plotDependence->savePng(tempPath, plotDependence->width(), plotDependence->height());
    QFile file(tempPath);
    if (file.open(QIODevice::ReadOnly)) {
      WasmFileIO::saveFile(file.readAll(), fileName);
      file.close();
    }
  }
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

    if (yvalues.size() > 0) {
      QByteArray data;
      QBuffer buffer(&data);
      buffer.open(QIODevice::WriteOnly | QIODevice::Text);
      QTextStream out(&buffer);
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
      for (size_t i = 0; i < yvalues.size(); ++i) {
        out.setFieldWidth(15);
        out.setFieldAlignment(QTextStream::AlignLeft);
        out << xvalues[i] << yvalues[i];
        out << qSetFieldWidth(0) << Qt::endl << qSetFieldWidth(15);
      }
      buffer.close();
      WasmFileIO::saveFile(data, fileName);
    }
  }
#else
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
        for (size_t i = 0; i < yvalues.size(); ++i) {
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
#endif
}

std::vector<double> CosmicEoSTab::getValues(int index, int num)
{
  int tsize = fCurrentSize;
  std::vector<double> ret(tsize, 0.);

  if (index >= 0 && index < paramnames.size() && paramnames[index] == "n_i^lep/T³") {
    int pid = comboParticle->currentIndex();
    if (num == 1) {
      pid = comboParticle2->currentIndex();
    }

    for (int j = 0; j < tsize; ++j) {
      ret[j] = paramsTD[j].densities[pid];
      ret[j] *= 1. / pow(paramsTD[j].T, 3) / xMath::GeVtoifm3();
    }

    return ret;
  }

  if (index >= 0 && index < paramnames.size() && paramnames[index] == "n_i^had/T³") {
    int pid = comboParticle->currentIndex();
    int feed = comboFeeddown->currentIndex();
    if (num == 1) {
      pid = comboParticle2->currentIndex();
      feed = comboFeeddown2->currentIndex();
    }

    for (int j = 0; j < tsize; ++j) {
      ret[j] = paramsTDHRG[j].densities[feed][pid];
      ret[j] *= 1. / pow(paramsTD[j].T, 3) / xMath::GeVtoifm3();
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

    tind = parammap["Δ"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].I / paramsTD[j].e / 3.;
      }
    }

    tind = parammap["p/T⁴_hrg"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTDHRG[j].pT4;
      }
    }

    tind = parammap["e/T⁴_hrg"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTDHRG[j].eT4;
      }
    }

    tind = parammap["(e-3P)/T⁴_hrg"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTDHRG[j].IT4;
      }
    }

    tind = parammap["s/T³_hrg"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTDHRG[j].sT3;
      }
    }

    tind = parammap["Δ_hrg"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTDHRG[j].I / paramsTDHRG[j].e / 3.;
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

    tind = parammap["μ_e[MeV]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].mu_e * 1.e3;
      }
    }

    tind = parammap["μ_μ[MeV]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].mu_mu * 1.e3;
      }
    }

    tind = parammap["μ_τ[MeV]"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].mu_tau * 1.e3;
      }
    }

    tind = parammap["ρB/s"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoB / paramsTD[j].s;
      }
    }

    tind = parammap["ρQ/s"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoQ / paramsTD[j].s;
      }
    }

    tind = parammap["ρ_e/s"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoE / paramsTD[j].s;
      }
    }

    tind = parammap["ρ_μ/s"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoMu / paramsTD[j].s;
      }
    }

    tind = parammap["ρ_τ/s"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].rhoTau / paramsTD[j].s;
      }
    }

  return ret;
}

std::vector<double> CosmicEoSTab::getValuesRatio(int index, int index2)
{
  std::vector<double> ret1 = getValues(index, 0), ret2 = getValues(index2, 1);
  std::vector<double> ret;
  for (int i = 0; i < ret1.size(); ++i)
    ret.push_back(ret1[i] / ret2[i]);
  return ret;
}

void CosmicEoSTab::replot() {
    bool flipAxes = false;
    flipAxes = CBflipAxes->isChecked();

    QCPAxis* xAxis = plotDependence->xAxis;
    QCPAxis* yAxis = plotDependence->yAxis;

    if (flipAxes) {
      xAxis = plotDependence->yAxis;
      yAxis = plotDependence->xAxis;
    }

    xAxis->setLabel(getParameterName());

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
        if ((paramnames[index] == "n_i^had/T³" && paramnames[index2] == "n_i^had/T³") ||
          (paramnames[index] == "n_i^had[fm^-3]" && paramnames[index2] == "n_i^had[fm^-3]")) {
          QString tname = paramnames[index] + "/" + paramnames[index2];
          int pid = comboParticle->currentIndex(), pid2 = comboParticle2->currentIndex();
          if (pid >= 0 && pid < model->TPS()->Particles().size()
            && pid2 >= 0 && pid2 < model->TPS()->Particles().size())
            tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) + ")/n(" 
                  + QString::fromStdString(model->TPS()->Particles()[pid2].Name()) + ")";
          yAxis->setLabel(tname);
        }
        if ((paramnames[index] == "n_i^lep/T³" && paramnames[index2] == "n_i^had/T³") ||
            (paramnames[index] == "n_i^lep[fm^-3]" && paramnames[index2] == "n_i^had[fm^-3]")) {
          QString tname = paramnames[index] + "/" + paramnames[index2];
          int pid = comboParticle->currentIndex(), pid2 = comboParticle2->currentIndex();
          if (pid >= 0 && pid < cosmos->NumberOfElectroWeakSpecies()
              && pid2 >= 0 && pid2 < model->TPS()->Particles().size())
            tname = "n(" + QString::fromStdString(cosmos->GetSpeciesName(pid)) + ")/n("
                    + QString::fromStdString(model->TPS()->Particles()[pid2].Name()) + ")";
          yAxis->setLabel(tname);
        }
        if ((paramnames[index] == "n_i^had/T³" && paramnames[index2] == "n_i^lep/T³") ||
            (paramnames[index] == "n_i^had[fm^-3]" && paramnames[index2] == "n_i^lep[fm^-3]")) {
          QString tname = paramnames[index] + "/" + paramnames[index2];
          int pid = comboParticle->currentIndex(), pid2 = comboParticle2->currentIndex();
          if (pid >= 0 && pid < model->TPS()->Particles().size()
              && pid2 >= 0 && pid2 < cosmos->NumberOfElectroWeakSpecies())
            tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) + ")/n("
                    + QString::fromStdString(cosmos->GetSpeciesName(pid2)) + ")";
          yAxis->setLabel(tname);
        }
        if ((paramnames[index] == "n_i^lep/T³" && paramnames[index2] == "n_i^lep/T³") ||
            (paramnames[index] == "n_i^lep[fm^-3]" && paramnames[index2] == "n_i^lep[fm^-3]")) {
          QString tname = paramnames[index] + "/" + paramnames[index2];
          int pid = comboParticle->currentIndex(), pid2 = comboParticle2->currentIndex();
          if (pid >= 0 && pid < cosmos->NumberOfElectroWeakSpecies()
              && pid2 >= 0 && pid2 < cosmos->NumberOfElectroWeakSpecies())
            tname = "n(" + QString::fromStdString(cosmos->GetSpeciesName(pid)) + ")/n("
                    + QString::fromStdString(cosmos->GetSpeciesName(pid2)) + ")";
          yAxis->setLabel(tname);
        }
        yvalues = getValuesRatio(index, index2);
      }
      else {
        yAxis->setLabel(paramnames[index]);
        if (paramnames[index] == "n_i^had/T³" || paramnames[index] == "n_i^had[fm^-3]") {
          QString tname = paramnames[index];
          int pid = comboParticle->currentIndex();
          if (pid >= 0 && pid < model->TPS()->Particles().size())
            tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) +")";
          if (paramnames[index] == "n_i^had/T³")
            tname += "/T³";
          else
            tname += "[fm^-3]";
          yAxis->setLabel(tname);
        }
        if (paramnames[index] == "n_i^lep/T³" || paramnames[index] == "n_i^lep[fm^-3]") {
          QString tname = paramnames[index];
          int pid = comboParticle->currentIndex();
          if (pid >= 0 && pid < cosmos->NumberOfElectroWeakSpecies())
            tname = "n(" + QString::fromStdString(cosmos->GetSpeciesName(pid)) +")";
          if (paramnames[index] == "n_i^lep/T³")
            tname += "/T³";
          else
            tname += "[fm^-3]";
          yAxis->setLabel(tname);
        }
        yvalues = getValues(index);

        if (CBxAxis->isChecked()) {
          xvalues = getValues(comboQuantity2->currentIndex(), 1);

          int index2 = comboQuantity2->currentIndex();
          xAxis->setLabel(paramnames[index2]);
          if (paramnames[index2] == "n_i^had/T³" || paramnames[index2] == "n_i^had[fm^-3]") {
            QString tname = paramnames[index2];
            int pid = comboParticle2->currentIndex();
            if (pid >= 0 && pid < model->TPS()->Particles().size())
              tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) +")";
            if (paramnames[index2] == "n_i^had/T³")
              tname += "/T³";
            else
              tname += "[fm^-3]";
            xAxis->setLabel(tname);
          }
          if (paramnames[index2] == "n_i^lep/T³" || paramnames[index2] == "n_i^lep[fm^-3]") {
            QString tname = paramnames[index2];
            int pid = comboParticle2->currentIndex();
            if (pid >= 0 && pid < cosmos->NumberOfElectroWeakSpecies())
              tname = "n(" + QString::fromStdString(cosmos->GetSpeciesName(pid)) +")";
            if (paramnames[index2] == "n_i^lep/T³")
              tname += "/T³";
            else
              tname += "[fm^-3]";
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

          if (xvalues.size() > 0)
            xAxis->setRange(xmin, xmax);
        }

        QCPCurve *curve = new QCPCurve(plotDependence->xAxis, plotDependence->yAxis);
        curve->setName("Model");
        curve->setPen(QPen(Qt::black, 2, Qt::SolidLine));
        curve->data()->clear();

      for(int i=0;i<yvalues.size();++i) {
          if (flipAxes)
              curve->addData(yvalues[i], xvalues[i]);
              //plotDependence->graph(graphNumber)->addData(yvalues[i], varvalues[i]);
          else
              curve->addData(xvalues[i], yvalues[i]);
            //plotDependence->graph(graphNumber)->addData(varvalues[i], yvalues[i]);
      }

      tmin = std::min(tmin, yAxis->range().lower);
      tmax = std::max(tmax, yAxis->range().upper);

      if (!CBxAxis->isChecked())
        xAxis->setRange(spinTMin->value(), spinTMax->value());
      if (yvalues.size() > 0)
          yAxis->setRange(1.1 * tmin, 1.1 * tmax);

      plotDependence->legend->setFont(QFont("Arial", font().pointSize() + 2));
      plotDependence->legend->setVisible(true);

      plotDependence->replot();
    }
}

QString CosmicEoSTab::getParameterName() const {
  return "T[MeV]";
}

void CosmicEoSTab::showEoSTable() {
  recomputeCalcTable();

#ifdef Q_OS_WASM
  CalculationTableDialog *dialog = new CalculationTableDialog(this, calcTable);
  dialog->setAttribute(Qt::WA_DeleteOnClose);
  dialog->setModal(true);
  dialog->showMaximized();
#else
  CalculationTableDialog dialog(this, calcTable);
  dialog.setWindowFlags(Qt::Window);
  dialog.showMaximized();
  dialog.exec();
#endif
}


void CosmicEoSTab::recomputeCalcTable() {
  calcTable.clear();

  calcTable.parameter_name = getParameterName();
  for(int i = 0; i < varvalues.size(); ++i) {
    calcTable.parameter_values.push_back(varvalues[i]);
    calcTable.temperature_values.push_back(paramsTD[i].T);
  }

  for(int ir = 0; ir < paramnames.size(); ++ir) {
    if (paramnames[ir] == "n_i^had/T³" || paramnames[ir] == "n_i^had[fm^-3]")
      continue;

    if (paramnames[ir] == "n_i^lep/T³" || paramnames[ir] == "n_i^lep[fm^-3]")
      continue;

    calcTable.quantities_names.push_back(paramnames[ir]);
    calcTable.quantities_values.push_back(getValues(ir));
  }

  for(int ic = 0; ic < model->TPS()->Particles().size(); ++ic) {
    calcTable.densities_names.push_back(QString::fromStdString(model->TPS()->Particles()[ic].Name()));
  }

  for(int ir = 0; ir < varvalues.size(); ++ir) {
    calcTable.densities_values.push_back(paramsTDHRG[ir].densities);
  }

  for(int ilep = 0; ilep < cosmos->NumberOfElectroWeakSpecies(); ++ilep) {
    calcTable.quantities_names.push_back(QString::fromStdString(
            "n(" + cosmos->GetSpeciesName(ilep) + ")/T³"));
    std::vector<double> lepton_densities;
    for(int ival = 0; ival < varvalues.size(); ++ival) {
      lepton_densities.push_back(paramsTD[ival].densities[ilep] / pow(paramsTD[ival].T * xMath::GeVtoifm(), 3.));
    }
    calcTable.quantities_values.push_back(lepton_densities);
  }
}