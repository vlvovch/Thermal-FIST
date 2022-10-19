/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "modeltab.h"

#include <QLayout>
#include <QLabel>
#include <QHeaderView>
#include <QGroupBox>
#include <QApplication>
#include <QMessageBox>
#include <QElapsedTimer>
#include <QFileDialog>
#include <QDebug>
#include <cstdio>
#include <vector>

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGEV/ThermalModelEVCrossterms.h"
#include "HRGVDW/ThermalModelVDW.h"
#include "HRGBase/ThermalModelCanonical.h"
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGEV/ThermalModelEVCanonicalStrangeness.h"
#include "HRGVDW/ThermalModelVDWCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonicalCharm.h"
#include "HRGPCE/ThermalModelPCE.h"

#include "DebugText.h"
#include "tablemodel.h"
#include "particledialog.h"
#include "resultdialog.h"
#include "fittoexperimenttab.h"
#include "correlationsdialog.h"

using namespace thermalfist;

ModelTab::ModelTab(QWidget *parent, ThermalModelBase *modelop)
    : QWidget(parent)
{
    dbgstrm.setString(&dbgstr);

    model = new ThermalModelIdeal(modelop->TPS());
    *model = *modelop;

  

    QHBoxLayout *DataEditLay = new QHBoxLayout();

    QVBoxLayout *dataLayv = new QVBoxLayout();

    myModel = new TableModel(0);
    myModel->setModel(model);
    tableParticles = new QTableView();
    tableParticles->setModel(myModel);
    tableParticles->setSelectionBehavior(QAbstractItemView::SelectRows);
    tableParticles->setSelectionMode(QAbstractItemView::SingleSelection);
    tableParticles->resizeColumnsToContents();
    connect(tableParticles->selectionModel(), SIGNAL(selectionChanged(QItemSelection,QItemSelection)), this, SLOT(changedRow()));
    connect(tableParticles, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(particleInfoDoubleClick(QModelIndex)));

    QHBoxLayout *layMisc = new QHBoxLayout();
    layMisc->setAlignment(Qt::AlignLeft);

    checkOnlyStable = new QCheckBox(tr("Show only stable particles"));
    connect(checkOnlyStable, SIGNAL(toggled(bool)), this, SLOT(switchStability(bool)));

    buttonResults = new QPushButton(tr("Equation of state..."));
    connect(buttonResults, SIGNAL(clicked()), this, SLOT(showResults()));

    buttonCorrelations = new QPushButton(tr("Correlations..."));
    connect(buttonCorrelations, SIGNAL(clicked()), this, SLOT(showCorrelations()));
    buttonCorrelations->setEnabled(false);

    labelHint = new QLabel(tr("Hint: double-click on particle for more info"));
    QFont tmpf = QApplication::font();
    tmpf.setPointSize(tmpf.pointSize() - 1);
    labelHint->setFont(tmpf);

    labelValid = new QPushButton(tr("Calculation valid!"));
    labelValid->setFlat(true);
    labelValid->setVisible(false);
    connect(labelValid, SIGNAL(clicked()), this, SLOT(showValidityCheckLog()));


    layMisc->addWidget(checkOnlyStable);
    layMisc->addWidget(buttonResults);
    layMisc->addWidget(buttonCorrelations);
    layMisc->addWidget(labelHint);
    layMisc->addStretch(1);
    layMisc->addWidget(labelValid, 0, Qt::AlignRight);

    dataLayv->addWidget(tableParticles);
    dataLayv->addLayout(layMisc);


    QVBoxLayout *editorLay = new QVBoxLayout();
    editorLay->setContentsMargins(15, 0, 0, 0);
    //editorLay->setMargin(10);
    editorLay->setAlignment(Qt::AlignTop);

    QGroupBox *grModelConfig = new QGroupBox(tr("HRG model configuration:"));

    configWidget = new ModelConfigWidget(NULL, model);
    connect(configWidget, SIGNAL(changed()), this, SLOT(modelChanged()));

    QHBoxLayout *layModelConfig = new QHBoxLayout();
    layModelConfig->setAlignment(Qt::AlignLeft);
    layModelConfig->addWidget(configWidget);


    grModelConfig->setLayout(layModelConfig);


    QGroupBox *grParameters = new QGroupBox(tr("Parameters:"));

    QGridLayout *layParameters = new QGridLayout();
    layParameters->setAlignment(Qt::AlignLeft);
    QLabel *labelTemperature = new QLabel(tr("T (MeV):"));
    spinTemperature = new QDoubleSpinBox();
    spinTemperature->setMinimum(0.);
    spinTemperature->setMaximum(10000.);
    spinTemperature->setValue(model->Parameters().T * 1e3);
    spinTemperature->setToolTip(tr("Temperature"));
    QLabel *labelmuB = new QLabel(tr("μ<sub>B</sub> (MeV):"));
    spinmuB = new QDoubleSpinBox();
    spinmuB->setMinimum(-10000.);
    spinmuB->setMaximum(10000.);
    spinmuB->setValue(model->Parameters().muB * 1e3);
    spinmuB->setToolTip(tr("Baryochemical potential"));
    QLabel *labelgammaq = new QLabel(tr("γ<sub>q</sub>:"));
    spingammaq = new QDoubleSpinBox();
    spingammaq->setMinimum(0.);
    spingammaq->setMaximum(10.);
    spingammaq->setDecimals(4);
    spingammaq->setValue(model->Parameters().gammaq);
    spingammaq->setToolTip(tr("Chemical non-equilibrium factor for light quarks"));
    labelgammaS = new QLabel(tr("γ<sub>S</sub>:"));
    spingammaS = new QDoubleSpinBox();
    spingammaS->setMinimum(0.);
    spingammaS->setMaximum(10.);
    spingammaS->setDecimals(4);
    spingammaS->setValue(model->Parameters().gammaS);
    spingammaS->setToolTip(tr("Chemical non-equilibrium factor for strange quarks"));
    labelgammaC = new QLabel(tr("γ<sub>C</sub>:"));
    spingammaC = new QDoubleSpinBox();
    spingammaC->setMinimum(0.);
    spingammaC->setMaximum(50.);
    spingammaC->setDecimals(4);
    spingammaC->setValue(model->Parameters().gammaC);
    spingammaC->setToolTip(tr("Chemical non-equilibrium factor for charm quarks"));

    labelmuS = new QLabel(tr("μ<sub>S</sub> (MeV):"));
    spinmuS = new QDoubleSpinBox();
    spinmuS->setMinimum(-1000.);
    spinmuS->setMaximum(1000.);
    spinmuS->setValue(model->Parameters().muS * 1e3);
    spinmuS->setToolTip(tr("Strangeness chemical potential"));
    QLabel *labelmuQ = new QLabel(tr("μ<sub>Q</sub> (MeV):"));
    spinmuQ = new QDoubleSpinBox();
    spinmuQ->setMinimum(-1000.);
    spinmuQ->setMaximum(1000.);
    spinmuQ->setValue(model->Parameters().muQ * 1e3);
    spinmuQ->setToolTip(tr("Electric charge chemical potential"));
    labelmuC = new QLabel(tr("μ<sub>C</sub> (MeV):"));
    spinmuC = new QDoubleSpinBox();
    spinmuC->setMinimum(-1000.);
    spinmuC->setMaximum(1000.);
    spinmuC->setValue(model->Parameters().muC * 1e3);
    spinmuC->setToolTip(tr("Charm chemical potential"));
    QLabel *labelVolumeR = new QLabel(tr("R (fm):"));
    spinVolumeR = new QDoubleSpinBox();
    spinVolumeR->setMinimum(0.);
    spinVolumeR->setMaximum(25.);
    spinVolumeR->setDecimals(4);
    spinVolumeR->setValue(8.);
    spinVolumeR->setToolTip(tr("System radius: the system volume is a sphere of this radius"));
    connect(spinVolumeR, SIGNAL(valueChanged(double)), this, SLOT(changeVolumeRSC(double)));
    QLabel *labelVolumeRSC = new QLabel(tr("R<sub>C</sub>:"));
    spinVolumeRSC = new QDoubleSpinBox();
    spinVolumeRSC->setDecimals(4);
    spinVolumeRSC->setMinimum(0.);
    spinVolumeRSC->setMaximum(25.);
    spinVolumeRSC->setValue(spinVolumeR->value());
    spinVolumeRSC->setEnabled(false);
    spinVolumeRSC->setToolTip(tr("Correlation radius: the (canonical) correlation volume is a sphere of this radius"));

    QLabel* labelVolume = new QLabel(tr("V (fm<sup>3</sup>):"));
    labelVolumeVal = new QLabel("4000");

    changeVolumeRSC(spinVolumeR->value());


    labelB = new QLabel(tr("B:"));
    spinB = new QSpinBox();
    spinB->setMinimum(-1000);
    spinB->setMaximum(1000);
    spinB->setValue(0);
    labelS = new QLabel(tr("S:"));
    spinS = new QSpinBox();
    spinS->setMinimum(-1000);
    spinS->setMaximum(1000);
    spinS->setValue(0);
    labelQ = new QLabel(tr("Q:"));
    spinQ = new QSpinBox();
    spinQ->setMinimum(-1000);
    spinQ->setMaximum(1000);
    spinQ->setValue(0);
    labelC = new QLabel(tr("C:"));
    spinC = new QSpinBox();
    spinC->setMinimum(-1000);
    spinC->setMaximum(1000);
    spinC->setValue(0);

    spinB->setToolTip(tr("Total baryon number in CE calculation"));
    spinQ->setToolTip(tr("Total electric charge in CE calculation"));
    spinS->setToolTip(tr("Total strangeness in CE calculation"));
    spinC->setToolTip(tr("Total charm in CE calculation"));

    layParameters->addWidget(labelTemperature, 0, 0, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinTemperature, 0, 1);
    layParameters->addWidget(labelmuB, 1, 0, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinmuB, 1, 1);
    layParameters->addWidget(labelgammaq, 0, 2, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spingammaq, 0, 3); 
    layParameters->addWidget(labelgammaS, 0, 4, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spingammaS, 0, 5);
    layParameters->addWidget(labelgammaC, 0, 6, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spingammaC, 0, 7);

    layParameters->addWidget(labelmuS, 1, 4, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinmuS, 1, 5);
    layParameters->addWidget(labelmuQ, 1, 2, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinmuQ, 1, 3);
    layParameters->addWidget(labelmuC, 1, 6, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinmuC, 1, 7); 
    layParameters->addWidget(labelVolumeR, 2, 0, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinVolumeR, 2, 1);
    layParameters->addWidget(labelVolumeRSC, 2, 2, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinVolumeRSC, 2, 3);
    layParameters->addWidget(labelVolume, 2, 4, 1, 1, Qt::AlignRight);
    layParameters->addWidget(labelVolumeVal, 2, 5);

    layParameters->addWidget(labelB, 3, 0, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinB, 3, 1);
    layParameters->addWidget(labelS, 3, 4, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinS, 3, 5);
    layParameters->addWidget(labelQ, 3, 2, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinQ, 3, 3);
    layParameters->addWidget(labelC, 3, 6, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinC, 3, 7);

    grParameters->setLayout(layParameters);


    QHBoxLayout *layFlags = new QHBoxLayout();
    layFlags->setAlignment(Qt::AlignLeft);
    checkFluctuations = new QCheckBox(tr("Compute fluctuations and correlations"));
    checkFluctuations->setChecked(false);
    checkFluctuations->setToolTip(tr("Compute fluctuation observables (viewable in \"Equation of state\" and \"Correlations\" dialogs)"));

    checkMuInitials = new QCheckBox(tr("Reset mu's"));
    checkMuInitials->setChecked(true);
    checkMuInitials->setToolTip(tr("Whether the current values of chemical potentials are used as initial conditions for conservation laws or they are reset"));

    layFlags->addWidget(checkFluctuations);
    layFlags->addWidget(checkMuInitials);

    QHBoxLayout *layButtons = new QHBoxLayout();
    layButtons->setAlignment(Qt::AlignLeft);

    buttonCalculate = new QPushButton(tr("Calculate"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));

    buttonCalculateFitted = new QPushButton(tr("Calculate from fit tab"));
    buttonCalculateFitted->setToolTip(tr("Use parameters found from a thermal fit in the other tab (faster and avoids the rounding errors due to spin boxes)"));
    connect(buttonCalculateFitted, SIGNAL(clicked()), this, SLOT(calculateFitted()));

    buttonWriteToFile = new QPushButton(tr("Write to file..."));
    connect(buttonWriteToFile, SIGNAL(clicked()), this, SLOT(writetofile()));
    buttonWriteToFile->setToolTip(tr("Writes total and primordial yields of all particles to file"));

    layButtons->addWidget(buttonCalculate);
    layButtons->addWidget(buttonCalculateFitted);
    layButtons->addWidget(buttonWriteToFile);

    buttonBenchmark = new QPushButton(tr("Benchmark"));
    connect(buttonBenchmark, SIGNAL(clicked()), this, SLOT(benchmark()));

    teDebug = new QTextEdit;
    teDebug->setReadOnly(true);

    editorLay->addWidget(grModelConfig);
    editorLay->addWidget(grParameters);
    editorLay->addLayout(layFlags);
    editorLay->addLayout(layButtons);
    editorLay->addWidget(teDebug);

    DataEditLay->addLayout(dataLayv, 1);
    DataEditLay->addLayout(editorLay, 0);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(DataEditLay);
    setLayout(mainLayout);

    modelChanged();

    tabFit = NULL;
}

ModelTab::~ModelTab()
{
}

int ModelTab::getCurrentRow()
{
    QModelIndexList selectedList = tableParticles->selectionModel()->selectedRows();
    if (selectedList.count()==0) return -1;
    return selectedList.at(0).row();
}

void ModelTab::changedRow()
{
    QModelIndexList selectedList = tableParticles->selectionModel()->selectedRows();
}

void ModelTab::particleInfoDoubleClick(const QModelIndex & index) {
    int row = index.row();
    if (row>=0) {
      labelHint->setVisible(false);
      ParticleDialog dialog(this, model, myModel->GetRowToParticle()[row]);
      dialog.setWindowFlags(Qt::Window);
      dialog.exec();
    }
}

void ModelTab::changeVolumeRSC(double VRSC)
{
  spinVolumeRSC->setValue(VRSC);

  double R = spinVolumeR->value();
  double V = 4. / 3. * xMath::Pi() * R * R * R;
  labelVolumeVal->setText(QString("%1").arg(V));
}


ThermalModelConfig ModelTab::getConfig()
{
  ThermalModelConfig ret = configWidget->updatedConfig();

  ret.T = spinTemperature->value() * 1.e-3;
  ret.muB = spinmuB->value() * 1.e-3;
  ret.muQ = spinmuQ->value() * 1.e-3;
  ret.muS = spinmuS->value() * 1.e-3;
  ret.muC = spinmuC->value() * 1.e-3;
  ret.gq = spingammaq->value();
  ret.gS = spingammaS->value();
  ret.gC = spingammaC->value();
  ret.VolumeR = spinVolumeR->value();
  ret.VolumeRSC = spinVolumeRSC->value();
  ret.B = spinB->value();
  ret.Q = spinQ->value();
  ret.S = spinS->value();
  ret.C = spinC->value();

  ret.ComputeFluctations = checkFluctuations->isChecked();
  ret.ResetMus = checkMuInitials->isChecked();

  return ret;
}

ThermalModelConfig ModelTab::getConfigFromFit(thermalfist::ThermalModelFit * fit, const ThermalModelConfig & configfit)
{
  ThermalModelConfig ret = configfit;

  ret.T = fit->Parameters().T.value;
  ret.muB = fit->Parameters().muB.value;
  ret.muQ = fit->Parameters().muQ.value;
  ret.muS = fit->Parameters().muS.value;
  ret.muC = fit->Parameters().muC.value;
  ret.gq = fit->Parameters().gammaq.value;
  ret.gS = fit->Parameters().gammaS.value;
  ret.gC = fit->Parameters().gammaC.value;
  ret.VolumeR = fit->Parameters().R.value;
  ret.VolumeRSC = fit->Parameters().Rc.value;
  ret.Tkin = fit->Parameters().Tkin.value;

  ret.ComputeFluctations = checkFluctuations->isChecked();
  ret.ResetMus = false;

  return ret;
}

void ModelTab::updateControlsWithConfig(const ThermalModelConfig & config)
{
  configWidget->setNewConfig(config);
  
  spinTemperature->setValue(config.T * 1.e3);
  spinmuB->setValue(config.muB * 1.e3);
  spinmuQ->setValue(config.muQ * 1.e3);
  spinmuS->setValue(config.muS * 1.e3);
  spinmuC->setValue(config.muC * 1.e3);
  spingammaq->setValue(config.gq);
  spingammaS->setValue(config.gS);
  spingammaC->setValue(config.gC);
  spinVolumeR->setValue(config.VolumeR);
  spinVolumeRSC->setValue(config.VolumeRSC);
  spinB->setValue(config.B);
  spinQ->setValue(config.Q);
  spinS->setValue(config.S);
  spinC->setValue(config.C);
  checkFluctuations->setChecked(config.ComputeFluctations);

  modelChanged();
}

void ModelTab::performCalculation(const ThermalModelConfig & config)
{
  QElapsedTimer timerc;
  timerc.start();
  
  ThermalModelBase *modelnew;

  if (config.ModelType == ThermalModelConfig::DiagonalEV) {
    modelnew = new ThermalModelEVDiagonal(model->TPS());
  }
  else if (config.ModelType == ThermalModelConfig::CrosstermsEV)
    modelnew = new ThermalModelEVCrossterms(model->TPS());  
  else if (config.ModelType == ThermalModelConfig::CE) {
    modelnew = new ThermalModelCanonical(model->TPS());
  }
  else if (config.ModelType == ThermalModelConfig::SCE)
    modelnew = new ThermalModelCanonicalStrangeness(model->TPS());
  else if (config.ModelType == ThermalModelConfig::EVSCE)
    modelnew = new ThermalModelEVCanonicalStrangeness(model->TPS());
  else if (config.ModelType == ThermalModelConfig::VDWSCE)
    modelnew = new ThermalModelVDWCanonicalStrangeness(model->TPS());
  else if (config.ModelType == ThermalModelConfig::QvdW) {
    modelnew = new ThermalModelVDW(model->TPS());
  }
  else if (config.ModelType == ThermalModelConfig::CCE)
    modelnew = new ThermalModelCanonicalCharm(model->TPS());
  else
    modelnew = new ThermalModelIdeal(model->TPS());

  myModel->setModel(modelnew);
  configWidget->setModel(modelnew);
  if (model != NULL) 
    delete model;

  model = modelnew;

  QElapsedTimer timer;
  timer.start();


  model->SetTemperature(config.T);
  model->SetBaryonChemicalPotential(config.muB);
  model->SetElectricChemicalPotential(config.muQ);
  model->SetStrangenessChemicalPotential(config.muS);
  model->SetCharmChemicalPotential(config.muC);
  model->SetGammaq(config.gq);
  model->SetGammaS(config.gS);
  model->SetGammaC(config.gC);
  model->SetVolumeRadius(config.VolumeR);
  model->SetCanonicalVolumeRadius(config.VolumeRSC);

  dbgstrm << "T\t= " << model->Parameters().T * 1.e3 << " MeV" << endl;

  SetThermalModelConfiguration(model, config);

  if (config.InteractionModel != ThermalModelConfig::InteractionIdeal)
    SetThermalModelInteraction(model, config);


  // If fluctuations are calculated within the CE one needs a twice larger range of quantum numbers
  if (config.ModelType == ThermalModelConfig::CE) {
    static_cast<ThermalModelCanonical*>(model)->CalculateQuantumNumbersRange(config.ComputeFluctations);
  }

  printf("Initialization time = %ld ms\n", static_cast<long int>(timerc.elapsed()));

  timerc.restart();
  model->ConstrainChemicalPotentials(checkMuInitials->isChecked());
  model->CalculatePrimordialDensities();
  printf("Densities time = %ld ms\n", static_cast<long int>(timerc.elapsed()));

  if (config.UsePCE) {
    dbgstrm << "Chemical freeze-out:" << endl;
  }
  if (model->TPS()->hasBaryons() && !model->IsConservedChargeCanonical(ConservedCharge::BaryonCharge))
    dbgstrm << "muB\t= " << model->Parameters().muB * 1.e3 << " MeV" << endl;

  if (model->TPS()->hasCharged() && !model->IsConservedChargeCanonical(ConservedCharge::ElectricCharge))
    dbgstrm << "muQ\t= " << model->Parameters().muQ * 1.e3 << " MeV" << endl;

  if (model->TPS()->hasStrange() && !model->IsConservedChargeCanonical(ConservedCharge::StrangenessCharge))
    dbgstrm << "muS\t= " << model->Parameters().muS * 1.e3 << " MeV" << endl;

  if (model->TPS()->hasCharmed() && !model->IsConservedChargeCanonical(ConservedCharge::CharmCharge))
    dbgstrm << "muC\t= " << model->Parameters().muC * 1.e3 << " MeV" << endl;

  if (config.ModelType == ThermalModelConfig::CE)
  {
    if (model->TPS()->hasBaryons() && model->IsConservedChargeCanonical(ConservedCharge::BaryonCharge))
      dbgstrm << "B\t= " << model->CalculateBaryonDensity() * model->Volume() << endl;
    if (model->TPS()->hasStrange() && model->IsConservedChargeCanonical(ConservedCharge::StrangenessCharge))
      dbgstrm << "S\t= " << model->CalculateStrangenessDensity() * model->Volume() << endl;
    if (model->TPS()->hasCharged() && model->IsConservedChargeCanonical(ConservedCharge::ElectricCharge))
      dbgstrm << "Q\t= " << model->CalculateChargeDensity() * model->Volume() << endl;
    if (model->TPS()->hasCharmed() && model->IsConservedChargeCanonical(ConservedCharge::CharmCharge))
      dbgstrm << "C\t= " << model->CalculateCharmDensity() * model->Volume() << endl;
  }
  dbgstrm << "gammaq\t= " << model->Parameters().gammaq << endl;
  dbgstrm << "gammaS\t= " << model->Parameters().gammaS << endl;
  if (model->TPS()->hasCharmed())
    dbgstrm << "gammaC\t= " << model->Parameters().gammaC << endl;
  dbgstrm << "V\t= " << model->Volume() << " fm^3" << endl;
  dbgstrm << endl;
  dbgstrm << "Particle density\t= " << model->CalculateHadronDensity() << " fm^-3" << endl;
  double nb = model->CalculateBaryonDensity();
  dbgstrm << "Net baryon density\t= " << nb << " fm^-3" << endl;
  dbgstrm << "Net baryon number\t= " << nb * model->Volume() << endl;
  if (model->TPS()->hasCharged())
    dbgstrm << "Net electric charge\t= " << model->CalculateChargeDensity() * model->Volume() << endl;
  if (model->TPS()->hasStrange())
    dbgstrm << "Net strangeness\t= " << model->CalculateStrangenessDensity() * model->Volume() << endl;
  if (model->TPS()->hasCharmed())
    dbgstrm << "Net charm\t= " << model->CalculateCharmDensity() * model->Volume() << endl;
  dbgstrm << "Absolute baryon number\t= " << model->AbsoluteBaryonDensity() * model->Volume() << endl;
  dbgstrm << "E/N\t\t= " << model->CalculateEnergyDensity() / model->CalculateHadronDensity() << endl;
  if (fabs(nb) > 1.e-10)
    dbgstrm << "E/Nb\t\t= " << model->CalculateEnergyDensity() / nb << endl;
  if (fabs(nb) > 1.e-10)
    dbgstrm << "S/B\t\t= " << model->CalculateEntropyDensity() / nb << endl;
  dbgstrm << "S/N\t\t= " << model->CalculateEntropyDensity() / model->CalculateHadronDensity() << endl;
  if (fabs(nb) > 1.e-10)
    dbgstrm << "Q/B\t\t= " << model->CalculateChargeDensity() / model->CalculateBaryonDensity() << endl;
  if (model->TPS()->hasStrange())
    dbgstrm << "S/|S|\t\t= " << model->CalculateStrangenessDensity() / model->CalculateAbsoluteStrangenessDensity() << endl;
  if (model->TPS()->hasCharmed())
    dbgstrm << "C/|C|\t\t= " << model->CalculateCharmDensity() / model->CalculateAbsoluteCharmDensity() << endl;
  if (model->InteractionModel() == ThermalModelBase::DiagonalEV)
    dbgstrm << "EV/V\t\t= " << model->CalculateEigenvolumeFraction() << endl;
  dbgstrm << endl;

  if (config.UsePCE) {
    timerc.restart();

    ThermalModelPCE modelpce(model);
    modelpce.SetStabilityFlags(ThermalModelPCE::ComputePCEStabilityFlags(model->TPS(), config.PCESahaForNuclei, config.PCEFreezeLongLived, config.PCEWidthCut));
    modelpce.SetChemicalFreezeout(model->Parameters(), model->ChemicalPotentials());
    modelpce.CalculatePCE(config.Tkin);

    printf("PCE time = %ld ms\n", static_cast<long int>(timerc.elapsed()));
  }

  timerc.restart();
  model->CalculateFeeddown();
  printf("Feeddown time = %ld ms\n", static_cast<long int>(timerc.elapsed()));

  timerc.restart();

  if (config.ComputeFluctations) {
    model->CalculateFluctuations();

    computeHigherOrderFluctuations();

    buttonCorrelations->setEnabled(true);
  }
  else {
    buttonCorrelations->setEnabled(false);
  }

  printf("Fluctuations time = %ld ms\n", static_cast<long int>(timerc.elapsed()));

  timerc.restart();

  

  if (config.UsePCE) {
    dbgstrm << "Kinetic freeze-out:" << endl;
    dbgstrm << "Tkin\t\t= " << config.Tkin * 1.e3 << " MeV" << endl;
    dbgstrm << "Vkin\t\t= " << model->Volume() << " fm^3" << endl;
    dbgstrm << "Vkin/Vch\t\t= " << model->Volume() / (4./3.*xMath::Pi()*pow(config.VolumeR,3)) << endl;
    if (model->InteractionModel() == ThermalModelBase::DiagonalEV)
      dbgstrm << "EV/V\t\t= " << model->CalculateEigenvolumeFraction() << endl;
    dbgstrm << endl;
  }

  dbgstrm << endl;
  dbgstrm << "Calculation time = " << timer.elapsed() << " ms" << endl;
  dbgstrm << "----------------------------------------------------------" << endl;
  teDebug->append(dbgstr);
  dbgstr.clear();

  printf("Finalizing time = %ld ms\n", static_cast<long int>(timerc.elapsed()));

  if (model->IsLastSolutionOK()) {
    labelValid->setText(tr("Calculation valid!"));
    labelValid->setStyleSheet("border : none; background-color : lightgreen;");
  }
  else {
    labelValid->setText(tr("Calculation NOT valid!"));
    labelValid->setStyleSheet("border : none; background-color : red;");
  }
  labelValid->setVisible(true);

  myModel->updateAll();
}

void ModelTab::calculate() {
  ThermalModelConfig config = getConfig();
  performCalculation(config);

  if (!(config.Ensemble == ThermalModelConfig::EnsembleCE)) {
    spinmuB->setValue(model->Parameters().muB * 1.e3);
    spinmuQ->setValue(model->Parameters().muQ * 1.e3);
    if (!(config.Ensemble == ThermalModelConfig::EnsembleSCE))
      spinmuS->setValue(model->Parameters().muS * 1.e3);
    if (!(config.Ensemble == ThermalModelConfig::EnsembleSCE) && !(config.Ensemble == ThermalModelConfig::EnsembleCCE))
      spinmuC->setValue(model->Parameters().muC * 1.e3);
  }
  return;
}

void ModelTab::calculateFitted()
{
  ThermalModelConfig config = getConfigFromFit(tabFit->Fit(), tabFit->LastUsedConfig());
  updateControlsWithConfig(config);
  performCalculation(config);
}

void ModelTab::writetofile() {
    if (model->IsCalculated()) {
        QString path = QFileDialog::getSaveFileName(this, tr("Save data to file"), QApplication::applicationDirPath() + "/output.dat", tr("*.dat"));
        if (path.length()>0)
        {
            FILE *f = fopen(path.toStdString().c_str(), "w");
            fprintf(f, "%20s%15s%8s%8s%8s%8s%8s%15s%15s\n", "Name", "PDGID", "Stable", "B", "Q", "S", "C", "Primary yield", "Total yield");
            for(int i=0;i<model->TPS()->Particles().size();++i) {
                fprintf(f, "%20s%15lld%8d%8d%8d%8d%8d%15E%15E\n",
                        model->TPS()->Particles()[i].Name().c_str(),
                        model->TPS()->Particles()[i].PdgId(),
                        model->TPS()->Particles()[i].IsStable(),
                        model->TPS()->Particles()[i].BaryonCharge(),
                        model->TPS()->Particles()[i].ElectricCharge(),
                        model->TPS()->Particles()[i].Strangeness(),
                        model->TPS()->Particles()[i].Charm(),
                        model->Densities()[i] * model->Volume(),
                        model->TotalDensities()[i] * model->Volume());
            }
            fclose(f);
        }
    }
}


void ModelTab::benchmark() {
}

void ModelTab::switchStability(bool showStable) {
  myModel->setOnlyStable(showStable);
}

void ModelTab::showResults() {
  ResultDialog dialog(this, model, &flucts);
  dialog.setWindowFlags(Qt::Window);
  dialog.exec();
}

void ModelTab::showCorrelations() {
  CorrelationsDialog dialog(this, model);
  dialog.setWindowFlags(Qt::Window);
  dialog.exec();
}

void ModelTab::showValidityCheckLog() {
  if (!model->IsLastSolutionOK()) {
    QMessageBox msgBox;
    msgBox.setText("There were some issues in calculation. See the log below:");
    msgBox.setDetailedText(model->ValidityCheckLog().c_str());
    msgBox.exec();
  }
}

void ModelTab::computeHigherOrderFluctuations()
{
  if (model->Ensemble() != ThermalModelBase::GCE) {
    flucts.flag = false;
    return;
  }
  
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

  flucts.flag = true;
}

void ModelTab::setModel(ThermalModelBase *modelop) {
    *model = *modelop;
    myModel->setModel(model);
    tableParticles->resizeColumnsToContents();
}

void ModelTab::updateModel() {
    myModel->reset();
}

void ModelTab::resetTPS() {
    model->ChangeTPS(model->TPS());
    myModel->reset();
    tableParticles->setColumnHidden(5, !model->TPS()->hasBaryons());
    tableParticles->setColumnHidden(6, !model->TPS()->hasCharged());
    tableParticles->setColumnHidden(7, !model->TPS()->hasStrange());
    tableParticles->setColumnHidden(8, !model->TPS()->hasCharmed());
    tableParticles->resizeColumnsToContents();

    labelValid->setVisible(false);

    configWidget->setModel(model);

    modelChanged();
}

void ModelTab::modelChanged()
{ 
  if (configWidget->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE) {
    spinmuB->setEnabled(true);
    spinmuS->setEnabled(true);
    spinmuQ->setEnabled(true);
    spinmuC->setEnabled(true);
    spinB->setEnabled(false);
    spinS->setEnabled(false);
    spinQ->setEnabled(false);
    spinC->setEnabled(false);
  }
  else {
    spinmuB->setEnabled(!configWidget->currentConfig.CanonicalB);
    spinmuS->setEnabled(!configWidget->currentConfig.CanonicalS);
    spinmuQ->setEnabled(!configWidget->currentConfig.CanonicalQ);
    spinmuC->setEnabled(!configWidget->currentConfig.CanonicalC);
    spinB->setEnabled(configWidget->currentConfig.CanonicalB);
    spinS->setEnabled(configWidget->currentConfig.CanonicalS);
    spinQ->setEnabled(configWidget->currentConfig.CanonicalQ);
    spinC->setEnabled(configWidget->currentConfig.CanonicalC);

    spinVolumeRSC->setEnabled(configWidget->currentConfig.CanonicalB ||
    configWidget->currentConfig.CanonicalS ||
    configWidget->currentConfig.CanonicalQ ||
    configWidget->currentConfig.CanonicalC);
  }
  
  if (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleSCE) {
    spinmuS->setEnabled(false);
    spinmuC->setEnabled(false);
    spinVolumeRSC->setEnabled(true);
  }
  else if (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleCCE) {
    spinmuC->setEnabled(false);
    spinVolumeRSC->setEnabled(true);
  }
  else if (configWidget->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE ||
    (!configWidget->currentConfig.CanonicalB &&
    !configWidget->currentConfig.CanonicalS &&
    !configWidget->currentConfig.CanonicalQ &&
    !configWidget->currentConfig.CanonicalC))
    spinVolumeRSC->setEnabled(false);
  
  spinmuS->setVisible(model->TPS()->hasStrange());
  spinmuC->setVisible(model->TPS()->hasCharmed());
  labelmuS->setVisible(model->TPS()->hasStrange());
  labelmuC->setVisible(model->TPS()->hasCharmed());
  spingammaS->setVisible(model->TPS()->hasStrange());
  spingammaC->setVisible(model->TPS()->hasCharmed());
  labelgammaS->setVisible(model->TPS()->hasStrange());
  labelgammaC->setVisible(model->TPS()->hasCharmed());

  labelQ->setVisible(model->TPS()->hasCharged());
  spinQ->setVisible(model->TPS()->hasCharged());
  labelS->setVisible(model->TPS()->hasStrange());
  spinS->setVisible(model->TPS()->hasStrange());
  labelC->setVisible(model->TPS()->hasCharmed());
  spinC->setVisible(model->TPS()->hasCharmed());
}

void ModelTab::updateFontSizes() {
  QFont tmpf = QApplication::font();
  tmpf.setPointSize(tmpf.pointSize() - 1);
  labelHint->setFont(tmpf);
}