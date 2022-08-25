/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "fittoexperimenttab.h"

#include <QLayout>
//#include <QFileDialog>
#include <QLabel>
#include <QHeaderView>
#include <QGroupBox>
#include <QApplication>
#include <QMessageBox>
#include <QTimer>
#include <QDebug>
#include <QFileDialog>
#include <algorithm>

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

#include "ItemDelegateCustom.h"

#include "quantitydialog.h"
#include "DebugText.h"
#include "quantitiesmodel.h"
#include "particledialog.h"
#include "resultdialog.h"
#include "chi2ProfileDialog.h"
#include "thermalfitplots.h"
#include "fitparametersmodel.h"

using namespace thermalfist;

void FitWorker::run()
{
     fTHMFit->PerformFit();
     emit calculated();
}

FitToExperimentTab::FitToExperimentTab(QWidget *parent, ThermalModelBase *modelop) :
    QWidget(parent)
{
    fitcopy = NULL;
    cpath = QString(ThermalFIST_INPUT_FOLDER) + "/data";
    dbgstrm.setString(&dbgstr);

    TPS = modelop->TPS();
    model = new ThermalModelIdeal(modelop->TPS());
    *model = *modelop;

    fitcopy = new ThermalModelFit(model);

    QHBoxLayout *DataEditLay = new QHBoxLayout();

    QVBoxLayout *dataLayv = new QVBoxLayout();

    QHBoxLayout *layoutTop = new QHBoxLayout;
    QLabel *labelQuantities = new QLabel(tr("Data to fit:"));
    labelHint = new QLabel(tr("Hint: double-click on yield to edit"));
    QFont tmpf = QApplication::font();
    tmpf.setPointSize(tmpf.pointSize() - 1);
    labelHint->setFont(tmpf);
    myModel = new QuantitiesModel(this, &quantities, fitcopy);
    tableQuantities = new QTableView();
    tableQuantities->setModel(myModel);
    tableQuantities->setSelectionBehavior(QAbstractItemView::SelectRows);
    tableQuantities->setSelectionMode(QAbstractItemView::SingleSelection);
    tableQuantities->setItemDelegate(new ItemDelegateCustom(this));
    tableQuantities->resizeColumnsToContents();
    connect(tableQuantities, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(quantityDoubleClick(QModelIndex)));
    
    layoutTop->addWidget(labelQuantities);
    layoutTop->addStretch(1);
    layoutTop->addWidget(labelHint, 0, Qt::AlignRight);

    QHBoxLayout *layEditQuantities = new QHBoxLayout();
    layEditQuantities->setAlignment(Qt::AlignLeft);
    buttonAddQuantity = new QPushButton(tr("Add quantity to fit..."));
    connect(buttonAddQuantity, SIGNAL(clicked()), this, SLOT(addQuantity()));
    buttonRemoveQuantity = new QPushButton(tr("Remove selected quantity from fit"));
    connect(buttonRemoveQuantity, SIGNAL(clicked()), this, SLOT(removeQuantityFromFit()));
    buttonLoadFromFile = new QPushButton(tr("Load data from file..."));
    connect(buttonLoadFromFile, SIGNAL(clicked()), this, SLOT(loadFromFile()));
    QPushButton *buttonSaveToFile = new QPushButton(tr("Save data to file..."));
    connect(buttonSaveToFile, SIGNAL(clicked()), this, SLOT(saveToFile()));

    layEditQuantities->addWidget(buttonAddQuantity);
    layEditQuantities->addWidget(buttonRemoveQuantity);
    layEditQuantities->addWidget(buttonLoadFromFile);
    layEditQuantities->addWidget(buttonSaveToFile);

    QLabel *labelParams = new QLabel(tr("Extracted parameters:"));
    tableParameters = new QTableWidget();

    QHBoxLayout *layPlots = new QHBoxLayout();
    layPlots->setAlignment(Qt::AlignLeft);

    QLabel *labelPlots = new QLabel(tr("Plots: "));
    buttonPlotYields = new QPushButton(tr("Yields"));
    connect(buttonPlotYields, SIGNAL(clicked()), this, SLOT(plotYields()));
    buttonPlotDeviations = new QPushButton(tr("Deviations"));
    connect(buttonPlotDeviations, SIGNAL(clicked()), this, SLOT(plotDeviations()));
    buttonPlotDataModel = new QPushButton(tr("Data/Model"));
    connect(buttonPlotDataModel, SIGNAL(clicked()), this, SLOT(plotDataModel()));
    buttonPlotDataVsModel = new QPushButton(tr("Data vs Model"));
    connect(buttonPlotDataVsModel, SIGNAL(clicked()), this, SLOT(plotDataVsModel()));

    layPlots->addWidget(labelPlots);
    layPlots->addWidget(buttonPlotYields);
    layPlots->addWidget(buttonPlotDeviations);
    layPlots->addWidget(buttonPlotDataModel);
    layPlots->addWidget(buttonPlotDataVsModel);

    QHBoxLayout *layMisc = new QHBoxLayout();
    layMisc->setAlignment(Qt::AlignLeft);

    buttonResults = new QPushButton(tr("Equation of state..."));
    connect(buttonResults, SIGNAL(clicked()), this, SLOT(showResults()));
    buttonChi2Profile = new QPushButton(tr("Chi2 profile..."));
    connect(buttonChi2Profile, SIGNAL(clicked()), this, SLOT(showChi2Profile()));

    labelValid = new QPushButton(tr("Calculation valid!"));
    labelValid->setFlat(true);
    labelValid->setVisible(false); 
    connect(labelValid, SIGNAL(clicked()), this, SLOT(showValidityCheckLog()));

    //layMisc->addWidget(checkOnlyStable);
    layMisc->addWidget(buttonResults);
    //layMisc->addWidget(buttonChi2Map);
    layMisc->addWidget(buttonChi2Profile);
    layMisc->addStretch(1);
    layMisc->addWidget(labelValid, 0, Qt::AlignRight);

    dataLayv->addLayout(layoutTop);
    dataLayv->addWidget(tableQuantities);
    dataLayv->addLayout(layEditQuantities);
    dataLayv->addWidget(labelParams);
    dataLayv->addWidget(tableParameters);
    dataLayv->addLayout(layPlots);
    dataLayv->addLayout(layMisc);


    QVBoxLayout *editorLay = new QVBoxLayout();
    editorLay->setContentsMargins(15, 0, 0, 0);
    editorLay->setAlignment(Qt::AlignTop);


    QGroupBox *grModelConfig = new QGroupBox(tr("HRG model configuration:"));

    configWidget = new ModelConfigWidget(NULL, model, false, true);
    connect(configWidget, SIGNAL(changed()), this, SLOT(modelChanged()));

    QHBoxLayout *layModelConfig = new QHBoxLayout();
    layModelConfig->setAlignment(Qt::AlignLeft);
    layModelConfig->addWidget(configWidget);


    grModelConfig->setLayout(layModelConfig);


    QHBoxLayout *layModelEnsemble = new QHBoxLayout();
    layModelEnsemble->setAlignment(Qt::AlignLeft);

    QGroupBox *grParameters = new QGroupBox(tr("Fit parameters:"));

    tableFitParameters = new QTableView();
    FitParametersModel *fitParametersModel = new FitParametersModel(0, &m_FitParameters);
    tableFitParameters->verticalHeader()->hide();
    tableFitParameters->setModel(fitParametersModel);
    tableFitParameters->setItemDelegate(new ItemDelegateCustom(this));
    tableFitParameters->resizeColumnsToContents();


    QHBoxLayout *layCE = new QHBoxLayout();
    layCE->setAlignment(Qt::AlignLeft);
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

    spinB->setToolTip(tr("Total baryon number in the CE calculation"));
    spinQ->setToolTip(tr("Total electric charge in the CE calculation"));
    spinS->setToolTip(tr("Total strangeness in the CE calculation"));
    spinC->setToolTip(tr("Total charm in the CE calculation"));


    QHBoxLayout *layFlags = new QHBoxLayout();
    layFlags->setAlignment(Qt::AlignLeft);

    checkFixRc = new QCheckBox(tr("Fix Vc/V:"));
    checkFixRc->setChecked(false);
    checkFixRc->setToolTip(tr("Fix the correlation volume to be proportional to the total volume, otherwise treat as fit parameter"));
    connect(checkFixRc, &QCheckBox::toggled, this, &FitToExperimentTab::modelChanged);

    spinVcV = new QDoubleSpinBox();
    spinVcV->setMinimum(0.);
    spinVcV->setMaximum(1000.);
    spinVcV->setDecimals(3);
    spinVcV->setValue(1.);

    layCE->addWidget(checkFixRc);
    layCE->addWidget(spinVcV);
    layCE->addSpacing(30);
    layCE->addWidget(labelB);
    layCE->addWidget(spinB);
    layCE->addWidget(labelQ);
    layCE->addWidget(spinQ);
    layCE->addWidget(labelS);
    layCE->addWidget(spinS);
    layCE->addWidget(labelC);
    layCE->addWidget(spinC);

    QVBoxLayout *layParams = new QVBoxLayout();
    layParams->addWidget(tableFitParameters);
    layParams->addLayout(layCE);
    grParameters->setLayout(layParams);


    QHBoxLayout *layButtons = new QHBoxLayout();
    layButtons->setAlignment(Qt::AlignLeft);

    buttonCalculate = new QPushButton(tr("Perform fit"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));

    QPushButton *buttonWriteToFile = new QPushButton(tr("Write to file..."));
    buttonWriteToFile->setToolTip(tr("Write three files with (i) all yields in tex format; (ii) all yields in ascii table format; (iii) fit log"));
    connect(buttonWriteToFile, SIGNAL(clicked()), this, SLOT(writetofile()));

    layButtons->addWidget(buttonCalculate);
    layButtons->addWidget(buttonWriteToFile);


    teDebug = new QTextEdit;
    teDebug->setReadOnly(true);

    editorLay->addWidget(grModelConfig);
    editorLay->addWidget(grParameters);
    editorLay->addLayout(layButtons);
    editorLay->addWidget(teDebug);

    DataEditLay->addLayout(dataLayv, 1);
    DataEditLay->addLayout(editorLay, 0);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(DataEditLay);
    setLayout(mainLayout);

    calcTimer = new QTimer(this);
    connect(calcTimer, SIGNAL(timeout()), this, SLOT(updateProgress()));

    modelChanged();

    QString datapathprefix = QString(ThermalFIST_INPUT_FOLDER) + "/data";
    quantities = ThermalModelFit::loadExpDataFromFile((QString(ThermalFIST_INPUT_FOLDER) + "/data/ALICE-PbPb2.76TeV-0-10-all.dat").toStdString());
    myModel->setQuantities(&quantities);
    tableQuantities->resizeColumnsToContents();

    lastconfig = getConfig();

    m_FitParameters = ThermalModelFitParameters(model->Parameters());
    m_FitParameters.R.value = 8.;
    m_FitParameters.Rc.value = 2.;

    modelChanged();
}



FitToExperimentTab::~FitToExperimentTab()
{
    if (fitcopy!=NULL) { delete fitcopy; }
    delete model;
}


ThermalModelConfig FitToExperimentTab::getConfig()
{
  ThermalModelConfig ret = configWidget->updatedConfig();

  ret.T = m_FitParameters.T.value;
  ret.muB = m_FitParameters.muB.value;
  ret.muQ = m_FitParameters.muQ.value;
  ret.muS = m_FitParameters.muS.value;
  ret.muC = m_FitParameters.muC.value;
  ret.gq = m_FitParameters.gammaq.value;
  ret.gS = m_FitParameters.gammaS.value;
  ret.gC = m_FitParameters.gammaC.value;
  ret.VolumeR = m_FitParameters.R.value;
  ret.VolumeRSC = m_FitParameters.Rc.value;
  ret.B = spinB->value();
  ret.Q = spinQ->value();
  ret.S = spinS->value();
  ret.C = spinC->value();

  ret.Tkin = m_FitParameters.Tkin.value;

  ret.ComputeFluctations = false;
  ret.ResetMus = true;

  return ret;
}

ThermalModelFitParameters FitToExperimentTab::getFitParameters()
{
  return m_FitParameters;
}

int FitToExperimentTab::getCurrentRow()
{
    QModelIndexList selectedList = tableQuantities->selectionModel()->selectedRows();
    if (selectedList.count()==0) return -1;
    return selectedList.at(0).row();
}

void FitToExperimentTab::changedRow()
{
    QModelIndexList selectedList = tableQuantities->selectionModel()->selectedRows();
}

void FitToExperimentTab::performFit(const ThermalModelConfig & config, const ThermalModelFitParameters & params)
{
  ThermalModelBase *modelnew;

  if (config.ModelType == ThermalModelConfig::DiagonalEV) {
    modelnew = new ThermalModelEVDiagonal(model->TPS());
  }
  else if (config.ModelType == ThermalModelConfig::CrosstermsEV) {
    modelnew = new ThermalModelEVCrossterms(model->TPS());
  }
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
  model->SetCanonicalVolumeRadius(config.VolumeRSC);

  SetThermalModelConfiguration(model, config);

  timer.start();

  ThermalModelFit *fit = new ThermalModelFit(model);

  std::cout << "Tkin = " << params.Tkin.value << "\n";

  fit->SetParameters(params);

  fit->FixVcOverV(checkFixRc->isChecked());
  fit->SetVcOverV(spinVcV->value());

  fit->UseTkin(config.UsePCE);
  fit->UseSahaForNuclei(config.PCESahaForNuclei);
  fit->PCEFreezeLongLived(config.PCEFreezeLongLived);
  fit->SetPCEWidthCut(config.PCEWidthCut);


  SetThermalModelInteraction(model, config);

  if (config.ModelType == ThermalModelConfig::CE) {
    static_cast<ThermalModelCanonical*>(model)->CalculateQuantumNumbersRange(false);
  }

  for (int i = 0; i<quantities.size(); ++i) {
    fit->AddQuantity(quantities[i]);
  }

  if (fitcopy != NULL) { delete fitcopy; fitcopy = NULL; }
  fitcopy = fit;

  myModel->setModel(fitcopy);

  FitWorker *wrk = new FitWorker(fit);

  connect(wrk, SIGNAL(calculated()), this, SLOT(finalize()));
  connect(wrk, SIGNAL(finished()), wrk, SLOT(deleteLater()));


  buttonCalculate->setEnabled(false);

  wrk->start();

  calcTimer->start(100);

  lastconfig = config;
}

void FitToExperimentTab::calculate() {
  performFit(getConfig(), getFitParameters());
  return;
}

void FitToExperimentTab::showResults() {
    ResultDialog dialog(this, model);
    dialog.setWindowFlags(Qt::Window);
    dialog.exec();
}

void FitToExperimentTab::showChi2Profile()
{
  chi2ProfileDialog dialog(this, TPS, getConfig(), getFitParameters(), quantities, fitcopy);
  dialog.setWindowFlags(Qt::Window);
  dialog.setMinimumSize(QSize(800, 400));
  dialog.exec();
}

void FitToExperimentTab::setModel(ThermalModelBase *modelop) {
    *model = *modelop;
}

void FitToExperimentTab::removeQuantityFromFit() {
    QModelIndexList selectedList = tableQuantities->selectionModel()->selectedRows();
    for(unsigned int i = 0; i < selectedList.count(); ++i)
        myModel->removeQuantity(selectedList.at(i).row());
}

void FitToExperimentTab::quantityDoubleClick(const QModelIndex & index) {
    int row = index.row();
    if (row>=0) {
      labelHint->setVisible(false);
      QuantityDialog dialog(this, model, &quantities[row]);
      dialog.setWindowFlags(Qt::Window);
      dialog.exec();
    }
}

void FitToExperimentTab::addQuantity() {
    myModel->addQuantity();
    QuantityDialog dialog(this, model, &quantities[quantities.size()-1]);
    dialog.setWindowFlags(Qt::Window);
    if (dialog.exec()==QDialog::Rejected) myModel->removeQuantity(quantities.size()-1);
}

void FitToExperimentTab::loadFromFile() {
    QString path = QFileDialog::getOpenFileName(this, tr("Open file with experimental data to fit"), cpath);
    if (path.length()>0)
    {
        quantities = ThermalModelFit::loadExpDataFromFile(path.toStdString());
        if (fitcopy != NULL)
          fitcopy->ClearModelData();
        myModel->setQuantities(&quantities);
        tableQuantities->resizeColumnsToContents();
        cpath = path;
    }
}

void FitToExperimentTab::saveToFile()
{
  QString path = QFileDialog::getSaveFileName(this, tr("Save experimental data to file"), cpath);
  if (path.length()>0)
  {
    ThermalModelFit::saveExpDataToFile(quantities, path.toStdString());
    cpath = path;
  }
}

void FitToExperimentTab::modelChanged()
{
  if (configWidget->currentConfig.Ensemble != ThermalModelConfig::EnsembleGCE) {
    checkFixRc->setEnabled(true);
    spinVcV->setEnabled(checkFixRc->isChecked());
  }
  else {
    checkFixRc->setEnabled(false);
    spinVcV->setEnabled(false);
  }


  spinQ->setVisible(model->TPS()->hasCharged());
  spinS->setVisible(model->TPS()->hasStrange());
  spinC->setVisible(model->TPS()->hasCharmed());
  labelQ->setVisible(model->TPS()->hasCharged());
  labelS->setVisible(model->TPS()->hasStrange());
  labelC->setVisible(model->TPS()->hasCharmed());
  if (configWidget->currentConfig.Ensemble != ThermalModelConfig::EnsembleCE) {
    spinB->setEnabled(false);
    spinQ->setEnabled(false);
    spinS->setEnabled(false);
    spinC->setEnabled(false);
    labelB->setEnabled(false);
    labelQ->setEnabled(false);
    labelS->setEnabled(false);
    labelC->setEnabled(false);
  }
  else {
    spinB->setEnabled(configWidget->currentConfig.CanonicalB);
    spinQ->setEnabled(configWidget->currentConfig.CanonicalQ);
    spinS->setEnabled(configWidget->currentConfig.CanonicalS);
    spinC->setEnabled(configWidget->currentConfig.CanonicalC);
    labelB->setEnabled(configWidget->currentConfig.CanonicalB);
    labelQ->setEnabled(configWidget->currentConfig.CanonicalQ);
    labelS->setEnabled(configWidget->currentConfig.CanonicalS);
    labelC->setEnabled(configWidget->currentConfig.CanonicalC);
  }


  if (fitcopy != NULL && fitcopy->FittedQuantities().size() == fitcopy->ModelDataSize() && fitcopy->ModelDataSize() > 0) {
    buttonPlotDeviations->setEnabled(true);
    buttonPlotDataModel->setEnabled(true);
    buttonPlotDataVsModel->setEnabled(true);
  }
  else {
    buttonPlotDeviations->setEnabled(false);
    buttonPlotDataModel->setEnabled(false);
    buttonPlotDataVsModel->setEnabled(false);
  }

  // If GCE and only ratios fitted then volume drops out
  {
    int nYields = 0;
    for (int i = 0; i < quantities.size(); ++i)
      if (quantities[i].type == FittedQuantity::Multiplicity)
        nYields++;
    tableFitParameters->setRowHidden(m_FitParameters.IndexByName("R"), (nYields == 0));
  }

  // If GCE, or Vc fixed to V, then correlation volume drops out
  tableFitParameters->setRowHidden(m_FitParameters.IndexByName("Rc"),
    (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleGCE)
    || (checkFixRc->isChecked()));

  // If full CE, S/B fixed, or no baryons then muB drops out
  tableFitParameters->setRowHidden(m_FitParameters.IndexByName("muB"),
    (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleCE
      && configWidget->currentConfig.CanonicalB)
    || (configWidget->currentConfig.ConstrainMuB)
    || !(model->TPS()->hasBaryons()));

  // If full CE, Q/B fixed, or no charged particles then muQ drops out
  tableFitParameters->setRowHidden(m_FitParameters.IndexByName("muQ"),
    (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleCE
      && configWidget->currentConfig.CanonicalQ)
    || (configWidget->currentConfig.ConstrainMuQ)
    || !(model->TPS()->hasCharged()));

  // If full CE, SCE, or S fixed to zero, or no strange particles then muS drops out
  tableFitParameters->setRowHidden(m_FitParameters.IndexByName("muS"),
    (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleCE
      && configWidget->currentConfig.CanonicalS)
    || (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleSCE)
    || (configWidget->currentConfig.ConstrainMuS)
    || !(model->TPS()->hasStrange()));

  // If not GCE, or C fixed to zero, or no charm particles then muC drops out
  tableFitParameters->setRowHidden(m_FitParameters.IndexByName("muC"),
    (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleGCE
      && configWidget->currentConfig.CanonicalC)
    || (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleSCE)
    || (configWidget->currentConfig.Ensemble == ThermalModelConfig::EnsembleCCE)
    || (configWidget->currentConfig.ConstrainMuC)
    || !(model->TPS()->hasCharmed()));

  

  // If no strangeness, then gammaS drops out (but beware of neutral particles with strange quarks!)
  tableFitParameters->setRowHidden(m_FitParameters.IndexByName("gammaS"),
    !(model->TPS()->hasStrange()));

  // If no charm, then gammaC drops out (but beware of neutral particles with charm quarks!)
  tableFitParameters->setRowHidden(m_FitParameters.IndexByName("gammaC"),
    !(model->TPS()->hasCharmed()));

  // If PCE not used, no Tkin is fitted
  tableFitParameters->setRowHidden(m_FitParameters.IndexByName("Tkin"),
    !(configWidget->currentConfig.UsePCE));

  //m_FitParameters.Tkin.value = configWidget->currentConfig.Tkin;

  //tableParameters->hideRow(1);
  //tableParameters->setRowHidden(0, true);
  //tableParameters->resizeColumnsToContents();
}

void FitToExperimentTab::resetTPS() {
    model->ChangeTPS(model->TPS());
    myModel->updateAll();

    labelValid->setVisible(false);

    configWidget->setModel(model);
}

void FitToExperimentTab::writetofile() {
    if (model->IsCalculated()) {
        QString path = QFileDialog::getSaveFileName(this, tr("Save data to file"), QApplication::applicationDirPath() + "/output.out", tr("*.out"));
        if (path.length()>0)
        {
            std::string tpath = path.toStdString();
            tpath = tpath.substr(0, tpath.size()-4);
            if (fitcopy!=NULL) {
                fitcopy->PrintYieldsLatexAll(std::string(tpath + ".tex"), tpath);
                fitcopy->PrintYieldsTable(std::string(tpath + ".out"));
                std::string cmmnt = "Thermal fit to " + cpath.toStdString() + " within " + fitcopy->model()->TAG();
                fitcopy->PrintFitLog(std::string(tpath + ".txt"), cmmnt);
            }
        }
    }
}

void FitToExperimentTab::updateProgress() {
  //dbgstrm << "T\t= " << model->Parameters().T * 1.e3 << " MeV" << endl;
  dbgstrm << "T\t= " << fitcopy->Parameters().T.value * 1.e3 << " MeV" << endl;
  if (!(model->Ensemble() == ThermalModelBase::CE)) {
    dbgstrm << "muB\t= " << fitcopy->Parameters().muB.value * 1.e3 << " MeV" << endl;
    if (model->TPS()->hasCharged())
      dbgstrm << "muQ\t= " << fitcopy->Parameters().muQ.value * 1.e3 << " MeV" << endl;
    if (!(model->Ensemble() == ThermalModelBase::SCE)
      && model->TPS()->hasStrange())
      dbgstrm << "muS\t= " << fitcopy->Parameters().muS.value * 1.e3 << " MeV" << endl;
    if (!(model->Ensemble() == ThermalModelBase::SCE)
      && !(model->Ensemble() == ThermalModelBase::CCE)
      && model->TPS()->hasCharmed())
      dbgstrm << "muC\t= " << fitcopy->Parameters().muC.value * 1.e3 << " MeV" << endl;
  }
  else {
      dbgstrm << "B\t= " << fitcopy->BT() << endl;
      dbgstrm << "S\t= " << fitcopy->ST() << endl;
      dbgstrm << "Q\t= " << fitcopy->QT() << endl;
      dbgstrm << "C\t= " << fitcopy->CT() << endl;
  }
  if (fitcopy->Parameters().gammaq.toFit == true)
    dbgstrm << "gammaq\t= " << fitcopy->Parameters().gammaq.value << endl;
  if (fitcopy->Parameters().gammaS.toFit == true)
    dbgstrm << "gammaS\t= " << fitcopy->Parameters().gammaS.value << endl;
  if (fitcopy->Parameters().gammaC.toFit == true)
    dbgstrm << "gammaC\t= " << fitcopy->Parameters().gammaC.value << endl;
  dbgstrm << "V\t= " << (4./3.) * xMath::Pi() * pow(fitcopy->Parameters().R.value,3) << " fm^3" << endl;
  if (model->Ensemble() == ThermalModelBase::SCE || model->Ensemble() == ThermalModelBase::CE)
    dbgstrm << "Vc\t= " << (4. / 3.) * xMath::Pi() * pow(fitcopy->Parameters().Rc.value, 3) << " fm^3" << endl;
  if (fitcopy->UseTkin())
    dbgstrm << "Tkin\t= " << fitcopy->Parameters().Tkin.value * 1.e3 << " MeV" << endl;
  dbgstrm << endl;

  dbgstrm << "Iteration\t= " << fitcopy->Iters() << endl;
  dbgstrm << "chi2/Ndf\t= " << fitcopy->Chi2() << "/" << fitcopy->Ndf() << " = " << fitcopy->Chi2() / fitcopy->Ndf() << endl;

  teDebug->clear();
  teDebug->append(dbgstr);
  dbgstr.clear();

  myModel->updateAll();
}

void FitToExperimentTab::finalize() {
  calcTimer->stop();
  ThermalModelFitParameters result = fitcopy->Parameters();

  dbgstrm << "T\t= " << result.T.value * 1.e3 << " MeV" << endl;

  if (model->TPS()->hasBaryons() && !model->IsConservedChargeCanonical(ConservedCharge::BaryonCharge))
    dbgstrm << "muB\t= " << model->Parameters().muB * 1.e3 << " MeV" << endl;

  if (model->TPS()->hasStrange() && !model->IsConservedChargeCanonical(ConservedCharge::StrangenessCharge))
    dbgstrm << "muS\t= " << model->Parameters().muS * 1.e3 << " MeV" << endl;

  if (model->TPS()->hasStrange() && !model->IsConservedChargeCanonical(ConservedCharge::ElectricCharge))
    dbgstrm << "muQ\t= " << model->Parameters().muQ * 1.e3 << " MeV" << endl;

  if (model->TPS()->hasCharmed() && !model->IsConservedChargeCanonical(ConservedCharge::CharmCharge))
    dbgstrm << "muC\t= " << model->Parameters().muC * 1.e3 << " MeV" << endl;


  if (model->Ensemble() == ThermalModelBase::CE) {
    if (model->IsConservedChargeCanonical(ConservedCharge::BaryonCharge))
      dbgstrm << "B\t= " << model->CalculateBaryonDensity()      * model->Volume() << endl;
    if (model->IsConservedChargeCanonical(ConservedCharge::StrangenessCharge))
      dbgstrm << "S\t= " << model->CalculateStrangenessDensity() * model->Volume() << endl;
    if (model->IsConservedChargeCanonical(ConservedCharge::ElectricCharge))
      dbgstrm << "Q\t= " << model->CalculateChargeDensity()      * model->Volume() << endl;
    if (model->IsConservedChargeCanonical(ConservedCharge::CharmCharge))
      dbgstrm << "C\t= " << model->CalculateCharmDensity()       * model->Volume() << endl;
  }
  dbgstrm << "gammaq\t= " << model->Parameters().gammaq << endl;
  if (model->TPS()->hasStrange())
    dbgstrm << "gammaS\t= " << model->Parameters().gammaS << endl;
  if (model->TPS()->hasCharmed())
    dbgstrm << "gammaC\t= " << model->Parameters().gammaC << endl;
  //dbgstrm << "V\t= " << model->Parameters().V << " fm^3" << endl;
  dbgstrm << "V\t= " << (4./3.) * xMath::Pi() * pow(result.R.value, 3) << " fm^3" << endl;
  if (fitcopy->UseTkin()) {
    dbgstrm << "Tkin\t= " << result.Tkin.value * 1.e3 << " MeV" << endl;
    dbgstrm << "Vkin\t= " << model->Parameters().V << " fm^3" << endl;
  }
  dbgstrm << endl;

  dbgstrm << "Particle density\t= " << model->CalculateHadronDensity() << " fm^-3" << endl;

  double nb = model->CalculateBaryonDensity();
  dbgstrm << "Net baryon density\t= " << nb << " fm^-3" << endl;
  dbgstrm << "Net baryon number\t= " << nb * model->Parameters().V << endl;
  if (model->TPS()->hasCharged())
    dbgstrm << "Net electric charge\t= " << model->CalculateChargeDensity() * model->Parameters().V << endl;
  if (model->TPS()->hasStrange())
    dbgstrm << "Net strangeness\t= " << model->CalculateStrangenessDensity() * model->Parameters().V << endl;
  if (model->TPS()->hasCharmed())
    dbgstrm << "Net charm\t= " << model->CalculateCharmDensity() * model->Parameters().V << endl;
  dbgstrm << "E/N\t\t= " << model->CalculateEnergyDensity() / model->CalculateHadronDensity() << endl;
  if (fabs(nb) > 1.e-10)
    dbgstrm << "S/B\t\t= " << model->CalculateEntropyDensity() / nb << endl;
  if (fabs(nb) > 1.e-10)
    dbgstrm << "Q/B\t\t= " << model->CalculateChargeDensity() / model->CalculateBaryonDensity() << endl;
  if (model->TPS()->hasStrange())
    dbgstrm << "S/|S|\t\t= " << model->CalculateStrangenessDensity() / model->CalculateAbsoluteStrangenessDensity() << endl;
  if (model->TPS()->hasCharmed())
    dbgstrm << "C/|C|\t\t= " << model->CalculateCharmDensity() / model->CalculateAbsoluteCharmDensity() << endl;
  if (model->InteractionModel() == ThermalModelBase::DiagonalEV)
    dbgstrm << "EV/V\t\t= " << model->CalculateEigenvolumeFraction() << endl;
  dbgstrm << endl;
  dbgstrm << "chi2/ndf\t\t= " << result.chi2ndf * fitcopy->Ndf() << "/" << fitcopy->Ndf() << " = " << result.chi2ndf << endl;
  dbgstrm << endl;

  // Data description accuracy
  {
    std::pair<double, double> accuracy = fitcopy->ModelDescriptionAccuracy();
    dbgstrm << "Model accuracy = (" << QString::number(100. * accuracy.first, 'f', 2) 
      << QString::fromUtf8(" ± ") 
      << QString::number(100. * accuracy.second, 'f', 2)  << ") %" << endl;

    dbgstrm << endl;
  }

  qint64 elapsedTime = timer.elapsed();
  if (elapsedTime < 1000)
    dbgstrm << "Calculation time = " << timer.elapsed() << " ms" << endl;
  else if (elapsedTime < 10000)
    dbgstrm << "Calculation time = " << QString::number(timer.elapsed()/1000., 'f', 2) << " s" << endl;
  else
    dbgstrm << "Calculation time = " << QString::number(timer.elapsed()/1000., 'f', 1) << " s" << endl;

  dbgstrm << "----------------------------------------------------------" << endl;
  teDebug->clear();
  teDebug->append(dbgstr);
  dbgstr.clear();

  myModel->updateAll();

  tableQuantities->resizeColumnsToContents();

  tableParameters->clearContents();

  tableParameters->setColumnCount(3);
  int RowCount = 30;

  if (model->IsConservedChargeCanonical(ConservedCharge::BaryonCharge))
    RowCount--;

  if (model->IsConservedChargeCanonical(ConservedCharge::StrangenessCharge))
    RowCount--;

  if (model->IsConservedChargeCanonical(ConservedCharge::ElectricCharge))
    RowCount--;

  if (model->IsConservedChargeCanonical(ConservedCharge::CharmCharge))
    RowCount--;

  //if (model->Ensemble() == ThermalModelBase::CE)
  //  RowCount -= 4;

  //if (model->Ensemble() == ThermalModelBase::SCE)
  //  RowCount -= 2;

  RowCount -= 1;

  tableParameters->setRowCount(RowCount);
  tableParameters->verticalHeader()->hide();
  tableParameters->setHorizontalHeaderItem(0, new QTableWidgetItem(tr("Parameter")));
  tableParameters->setHorizontalHeaderItem(1, new QTableWidgetItem(tr("Value")));
  tableParameters->setHorizontalHeaderItem(2, new QTableWidgetItem(tr("Error")));

  int cRow = 0;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("T (MeV)"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.T.value*1.e3)));
  if (result.T.toFit==true) 
    tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.T.error*1.e3)));
  else
    tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));

  if (!model->IsConservedChargeCanonical(ConservedCharge::BaryonCharge)
    && model->TPS()->hasBaryons()) {
    cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("μB (MeV)"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.muB.value*1.e3)));
    if (result.muB.toFit==true) 
      tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.muB.error*1.e3)));
    else
      tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));
  }

  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("γq"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.gammaq.value)));
  if (result.gammaq.toFit==true) 
    tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.gammaq.error)));
  else
    tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));
  
  if (model->TPS()->hasStrange()) {
    cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("γs"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.gammaS.value)));
    if (result.gammaS.toFit == true)
      tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.gammaS.error)));
    else
      tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));
  }

  if (model->TPS()->hasCharmed()) {
    cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("γC"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.gammaC.value)));
    if (result.gammaC.toFit == true)
      tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.gammaC.error)));
    else
      tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));
  }

  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("R (fm)"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.R.value)));
  if (result.R.toFit==true) 
    tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.R.error)));
  else
    tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));

  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("V (fm^3)"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(4. * xMath::Pi() / 3. * pow(result.R.value, 3))));
  if (result.R.toFit == true)
    tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(4. * xMath::Pi() * pow(result.R.value, 2) * result.R.error)));
  else
    tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));

  if (model->Ensemble() != ThermalModelBase::GCE) {
    cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("Rc (fm)"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.Rc.value)));
    if (result.Rc.toFit==true) 
      tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.Rc.error)));
    else
      tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));

    cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("Vc (fm^3)"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(4. * xMath::Pi() / 3. * pow(result.Rc.value, 3))));
    if (result.R.toFit == true)
      tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(4. * xMath::Pi() * pow(result.Rc.value, 2) * result.Rc.error)));
    else
      tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));
  }

  if (model->UsePartialChemicalEquilibrium()) {
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("Tkin (MeV)"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.Tkin.value * 1.e3)));
    if (result.Tkin.toFit == true)
      tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.Tkin.error * 1.e3)));
    else
      tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));
  }

  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("chi2/dof"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.chi2) + "/" + QString::number(result.ndf)));
  tableParameters->setItem(cRow, 2, new QTableWidgetItem(""));

  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("chi2/dof"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.chi2ndf)));
  tableParameters->setItem(cRow, 2, new QTableWidgetItem(""));

  
  if (!model->IsConservedChargeCanonical(ConservedCharge::ElectricCharge)
    && model->TPS()->hasCharged()) {
      cRow++;
      tableParameters->setItem(cRow, 0, new QTableWidgetItem("μQ (MeV)"));
      tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.muQ.value*1.e3)));
      tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.muQ.error*1.e3)));
  }

  if (!model->IsConservedChargeCanonical(ConservedCharge::StrangenessCharge)
    && model->TPS()->hasStrange()) {
    cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("μS (MeV)"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.muS.value*1.e3)));
    tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.muS.error*1.e3)));
  }

  if (!model->IsConservedChargeCanonical(ConservedCharge::CharmCharge)
    && model->TPS()->hasCharmed()) {
    cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("μC (MeV)"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.muC.value*1.e3)));
    tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.muC.error*1.e3)));
  }

  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("nH (fm^-3)"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().nH.value)));
  tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().nH.error)));
    
  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("rhoB (fm^-3)"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().rhoB.value)));
  tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().rhoB.error)));
    
  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("rhoQ (fm^-3)"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().rhoQ.value)));
  tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().rhoQ.error)));
    
  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("e (MeV/fm^3)"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().en.value*1.e3)));
  tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().en.error*1.e3)));
    
  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("s (fm^-3)"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().entropy.value)));
  tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().entropy.error)));
    
  cRow++;
  tableParameters->setItem(cRow, 0, new QTableWidgetItem("p (MeV/fm^3)"));
  tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().pressure.value*1.e3)));
  tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().pressure.error*1.e3)));

  tableParameters->setRowCount(cRow + 1);

  tableParameters->resizeColumnsToContents();


  buttonCalculate->setEnabled(true);

  if (model->IsLastSolutionOK()) {
    labelValid->setText(tr("Calculation valid!"));
    labelValid->setStyleSheet("border : none; background-color : lightgreen;");
  }
  else {
    labelValid->setText(tr("Calculation NOT valid!"));
    labelValid->setStyleSheet("border : none; background-color : red;");
  }
  labelValid->setVisible(true);

  fitcopy->PrintYieldsLatexAll("Yield.dat", "p+p");

  modelChanged();
}

void FitToExperimentTab::showValidityCheckLog() {
  if (!model->IsLastSolutionOK()) {
    QMessageBox msgBox;
    msgBox.setText("There were some issues in calculation. See the log below:");
    msgBox.setDetailedText(model->ValidityCheckLog().c_str());
    msgBox.exec();
  }
}

void FitToExperimentTab::updateFontSizes() {
  QFont tmpf = QApplication::font();
  tmpf.setPointSize(tmpf.pointSize() - 1);
  labelHint->setFont(tmpf);
}

void FitToExperimentTab::plotYields()
{
  if (fitcopy != NULL) {
    fitcopy->SetQuantities(quantities);
    PlotDialog *dialog = new PlotDialog(this, new YieldsPlot(0, fitcopy));
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setWindowFlags(Qt::Window);
    dialog->setMinimumSize(QSize(800, 600));
    dialog->setModal(false);
    dialog->show();
    dialog->raise();
    dialog->activateWindow();
  }
}

void FitToExperimentTab::plotDeviations()
{
  if (fitcopy != NULL) {
    PlotDialog *dialog = new PlotDialog(this, new DeviationsPlot(0, fitcopy));
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setWindowFlags(Qt::Window);
    dialog->setMinimumSize(QSize(800, 300));
    dialog->setModal(false);
    dialog->show();
    dialog->raise();
    dialog->activateWindow();
  }
}

void FitToExperimentTab::plotDataModel()
{
  if (fitcopy != NULL) {
    PlotDialog *dialog = new PlotDialog(this, new DataModelPlot(0, fitcopy));
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setWindowFlags(Qt::Window);
    dialog->setMinimumSize(QSize(800, 300));
    dialog->setModal(false);
    dialog->show();
    dialog->raise();
    dialog->activateWindow();
  }
}

void FitToExperimentTab::plotDataVsModel()
{
  if (fitcopy != NULL) {
    PlotDialog *dialog = new PlotDialog(this, new DataVsModelPlot(0, fitcopy));
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setWindowFlags(Qt::Window);
    dialog->setMinimumSize(QSize(800, 600));
    dialog->setModal(false);
    dialog->show();
    dialog->raise();
    dialog->activateWindow();
  }
}

