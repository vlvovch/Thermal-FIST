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

using namespace std;
using namespace thermalfist;

void EoSWorker::run() {
  int i = 0;
  for (double T = Tmin; T <= Tmax + 1.e-10 && !(*stop); T += dT, ++i) {
    if (i < paramsTD->size()) {
      model->SetTemperature(T * 1.e-3);

      if (T == Tmin)
        model->ConstrainChemicalPotentials(true);
      else
        model->ConstrainChemicalPotentials(false);

      model->CalculateDensities();

      varvalues->operator [](i) = T;
      Thermodynamics tTD;
      tTD.pT4 = model->Pressure() / pow(T * 1.e-3, 4) / thermalfist::xMath::GeVtoifm3();
      tTD.eT4 = model->EnergyDensity() / pow(T * 1.e-3, 4) / thermalfist::xMath::GeVtoifm3();
      tTD.sT3 = model->EntropyDensity() / pow(T * 1.e-3, 3) / thermalfist::xMath::GeVtoifm3();
      tTD.IT4 = tTD.eT4 - 3. * tTD.pT4;
      tTD.nhT3 = model->HadronDensity() / pow(T * 1.e-3, 3) / thermalfist::xMath::GeVtoifm3();
      tTD.T = model->Parameters().T;
      tTD.muB = model->Parameters().muB;
      tTD.muQ = model->Parameters().muQ;
      tTD.muS = model->Parameters().muS;
      tTD.muC = model->Parameters().muC;
      tTD.densities = model->AllDensities();
      tTD.flag = true;
      paramsTD->operator [](i) = tTD;

      model->CalculateTwoParticleCorrelations();
      model->CalculateSusceptibilityMatrix();

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

    tname = "(e-3P)/T⁴";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "s/T³";
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

    CBratio = new QCheckBox(tr("Ratio"));
    CBratio->setChecked(false);
    connect(CBratio, SIGNAL(toggled(bool)), this, SLOT(replot()));
    connect(CBratio, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));

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

    layLeft->addWidget(CBratio, 0, Qt::AlignLeft);
    layLeft->addLayout(selLay);
    layLeft->addWidget(plotDependence);

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


    QGroupBox *grRange = new QGroupBox(tr("Parameter range:"));

    QVBoxLayout *layParamBounds = new QVBoxLayout();
    QGridLayout *layBounds = new QGridLayout();
    layBounds->setAlignment(Qt::AlignLeft);

    QLabel *labelmuB = new QLabel("μ<sub>B</sub> (MeV):");
    spinmuB = new QDoubleSpinBox();
    spinmuB->setDecimals(3);
    spinmuB->setMinimum(-10000.);
    spinmuB->setMaximum(10000.);
    spinmuB->setValue(0.);

    QLabel *labelTMin = new QLabel(tr("T<sub>min</sub> (MeV):"));
    spinTMin = new QDoubleSpinBox();
    spinTMin->setDecimals(3);
    spinTMin->setMinimum(0.01);
    spinTMin->setMaximum(10000.);
    spinTMin->setValue(100.);

    QLabel *labelTMax = new QLabel(tr("T<sub>max</sub> (MeV):"));
    spinTMax = new QDoubleSpinBox();
    spinTMax->setDecimals(3);
    spinTMax->setMinimum(0.01);
    spinTMax->setMaximum(10000.);
    spinTMax->setValue(200.);

    QLabel *labeldT = new QLabel(tr("∆T (MeV):"));
    spindT = new QDoubleSpinBox();
    spindT->setDecimals(3);
    spindT->setMinimum(0.);
    spindT->setMaximum(10000.);
    spindT->setValue(5.);


    layBounds->addWidget(labelTMin, 0, 0);
    layBounds->addWidget(spinTMin, 0, 1);
    layBounds->addWidget(labelTMax, 0, 2);
    layBounds->addWidget(spinTMax, 0, 3);
    layBounds->addWidget(labeldT, 0, 4);
    layBounds->addWidget(spindT, 0, 5);
    layBounds->addWidget(labelmuB, 1, 0);
    layBounds->addWidget(spinmuB, 1, 1);

    labelConstr = new QLabel(tr(""));


    layParamBounds->addLayout(layBounds);
    layParamBounds->addWidget(labelConstr, 0, Qt::AlignLeft);

    grRange->setLayout(layParamBounds);

    buttonCalculate = new QPushButton(tr("Calculate"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));

    layRight->addWidget(grModelConfig);
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

  ret.T = spinTMin->value() * 1.e-3;
  ret.muB = spinmuB->value() * 1.e-3;
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

      paramsTD.resize(0);
      paramsFl.resize(0);
      varvalues.resize(0);

      int iters = static_cast<int>((spinTMax->value() - spinTMin->value()) / spindT->value()) + 1;

      fCurrentSize = 0;

      paramsTD.resize(iters);
      paramsFl.resize(iters);
      varvalues.resize(iters);

      fStop = 0;

      EoSWorker *wrk = new EoSWorker(model, spinTMin->value(), spinTMax->value(), spindT->value(),
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
  if (CBratio->isChecked()) {
    comboQuantity2->setEnabled(true);
  }
  else {
    comboQuantity2->setEnabled(false);
  }

  if (comboQuantity->currentText() == "ni/T³") {
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

  if (comboQuantity2->currentText() == "ni/T³"
    && CBratio->isChecked()) {
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

  if (configWidget->currentConfig.ConstrainMuB) {
    spinmuB->setEnabled(false);
  }
  else {
    spinmuB->setEnabled(true);
  }

  QString constr = "";

  if (configWidget->currentConfig.ConstrainMuB && model->TPS()->hasBaryons()) {
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
}

void EquationOfStateTab::resetTPS()
{
  model->ChangeTPS(model->TPS());
  fillParticleLists();

  configWidget->setModel(model);
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

  menu->addAction("Save as pdf", this, SLOT(saveAsPdf()));
  menu->addAction("Save as png", this, SLOT(saveAsPng()));
  menu->addAction("Write computed values to file", this, SLOT(saveAsAscii()));

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

      if (CBratio->isChecked()) {
        yvalues = getValuesRatio(comboQuantity->currentIndex(), comboQuantity2->currentIndex());
      }
      else {
        yvalues = getValues(comboQuantity->currentIndex());
      }

      QFile fout(path);

      if (yvalues.size() > 0 && fout.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&fout);
        out.setFieldWidth(15);
        out.setFieldAlignment(QTextStream::AlignLeft);
        out << plotDependence->xAxis->label();
        out << plotDependence->yAxis->label();
        out << qSetFieldWidth(0) << endl << qSetFieldWidth(15);
        for (int i = 0; i < yvalues.size(); ++i) {
          out.setFieldWidth(15);
          out.setFieldAlignment(QTextStream::AlignLeft);
          out << varvalues[i] << yvalues[i];
          out << qSetFieldWidth(0) << endl << qSetFieldWidth(15);
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
  
  for (int i = 0; i < paramnames.size(); ++i) {
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

    tind = parammap["ntot/T³"];
    if (tind == index) {
      for (int j = 0; j < tsize; ++j) {
        ret[j] = paramsTD[j].nhT3;
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
  }

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

void EquationOfStateTab::replot() {
    int index = comboQuantity->currentIndex();

    plotDependence->clearGraphs();

    plotDependence->xAxis->setLabelFont(QFont("Arial", font().pointSize() + 4));
    plotDependence->yAxis->setLabelFont(QFont("Arial", font().pointSize() + 4));

    if (index >= 0 && index < paramnames.size()) {
      std::vector<double> yvalues;

      if (CBratio->isChecked()) {
        int index2 = comboQuantity2->currentIndex();
        if (!(index2 >= 0 && index2 < paramnames.size()))
          return;
        plotDependence->yAxis->setLabel(paramnames[index] + "/" + paramnames[index2]);
        if (paramnames[index] == "ni/T³" && paramnames[index2] == "ni/T³") {
          QString tname = paramnames[index] + "/" + paramnames[index2];
          int pid = comboParticle->currentIndex(), pid2 = comboParticle2->currentIndex();
          if (pid >= 0 && pid < model->TPS()->Particles().size()
            && pid2 >= 0 && pid2 < model->TPS()->Particles().size())
            tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) + ")/n(" 
                  + QString::fromStdString(model->TPS()->Particles()[pid2].Name()) + ")";
          plotDependence->yAxis->setLabel(tname);
        }
        yvalues = getValuesRatio(index, index2);
      }
      else {
        plotDependence->yAxis->setLabel(paramnames[index]);
        if (paramnames[index] == "ni/T³") {
          QString tname = paramnames[index];
          int pid = comboParticle->currentIndex();
          if (pid >= 0 && pid < model->TPS()->Particles().size())
            tname = "n(" + QString::fromStdString(model->TPS()->Particles()[pid].Name()) +")/T³";
          plotDependence->yAxis->setLabel(tname);
        }
        yvalues = getValues(index);
      }
      double tmin = 0., tmax = 0.;
      for (int i = 0; i<yvalues.size(); ++i) {
        tmin = std::min(tmin, yvalues[i]);
        tmax = std::max(tmax, yvalues[i]);
      }

      if (tmin > 0.)
        tmin *= 0.;

      if (yvalues.size() > 0)
        plotDependence->yAxis->setRange(tmin, tmax);

      // Lattice data for muB = 0 only
      if (spinmuB->value() == 0.0 && !model->ConstrainMuB())
        plotLatticeData();

      int graphNumber = plotDependence->graphCount();

      plotDependence->addGraph();
      plotDependence->graph(graphNumber)->setName("Model");
      plotDependence->graph(graphNumber)->setPen(QPen(Qt::black, 2, Qt::SolidLine));

      plotDependence->graph(graphNumber)->data()->clear();

      for(int i=0;i<yvalues.size();++i) {
          plotDependence->graph(graphNumber)->addData(varvalues[i], yvalues[i]);
      }

      

      tmin = std::min(tmin, plotDependence->yAxis->range().lower);
      tmax = std::max(tmax, plotDependence->yAxis->range().upper);

      plotDependence->xAxis->setRange(spinTMin->value(), spinTMax->value());
      if (yvalues.size() > 0)
        plotDependence->yAxis->setRange(1.1*tmin, 1.1*tmax);

      plotDependence->legend->setFont(QFont("Arial", font().pointSize() + 2));
      plotDependence->legend->setVisible(true);

      plotDependence->replot();
    }
}
