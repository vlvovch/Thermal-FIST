/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "chi2ProfileDialog.h"

#include <QLayout>
#include <QLabel>
#include <QDebug>

#include <fstream>
#include <iomanip>

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

#include "qcustomplot.h"

using namespace thermalfist;

chi2ProfileDialog::chi2ProfileDialog(QWidget *parent, ThermalParticleSystem *inTPS, const ThermalModelConfig & inConfig, const ThermalModelFitParameters & inParams, const std::vector<FittedQuantity> & inQuantities, ThermalModelFit *inFit) :
    QDialog(parent), TPS(inTPS), config(inConfig), fitParams(inParams), quantities(inQuantities), modelFitInput(inFit)
{
  fitParams = inFit->Parameters();
  
  model = NULL;
  modelFit = NULL;

  lastFilePath = "";

  QHBoxLayout *layout = new QHBoxLayout();

  QHBoxLayout *layCombo = new QHBoxLayout();

  QLabel *labelCombo = new QLabel(tr("Fit parameter"));
  comboParameter = new QComboBox();
  std::vector<std::string> paramNames;
  paramNames.push_back("T");
  paramNames.push_back("R");
  paramNames.push_back("Rc");
  paramNames.push_back("muB");
  paramNames.push_back("muQ");
  paramNames.push_back("muS");
  paramNames.push_back("muC");
  paramNames.push_back("gammaq");
  paramNames.push_back("gammaS");
  paramNames.push_back("gammaC");

  int pind = 0;
  for (int i = 0; i < paramNames.size(); ++i) {
    if (fitParams.GetParameter(paramNames[i]).toFit == true) {
      if (config.ModelType == ThermalModelConfig::CE && 
        (paramNames[i] == "muB" || paramNames[i] == "muS" || paramNames[i] == "muQ" || paramNames[i] == "muC"))
        continue;

      if ((config.ModelType == ThermalModelConfig::SCE || config.ModelType == ThermalModelConfig::EVSCE || config.ModelType == ThermalModelConfig::VDWSCE) &&
        (paramNames[i] == "muS"))
        continue;

      comboParameter->addItem(QString(paramNames[i].c_str()));
      paramNamesMap[pind++] = paramNames[i];

      double mn = 1.;
      std::string ParameterName = paramNames[i];
      if (ParameterName == "T" || ParameterName == "muB" ||
        ParameterName == "muQ" || ParameterName == "muS" ||
        ParameterName == "muC")
        mn = 1.e3;

      const FitParameter& param = modelFitInput->Parameters().GetParameter(paramNames[i]);

      
      vecAleft.push_back((param.value - 4. * param.errm) * mn);
      vecAright.push_back((param.value + 4. * param.errp) * mn);

      //vecAleft.push_back(fitParams.GetParameter(paramNames[i]).xmin * mn);
      //vecAright.push_back(fitParams.GetParameter(paramNames[i]).xmax * mn);

      vecAvalues.push_back(std::vector<double>(0));
      vecParams.push_back(std::vector<double>(0));
      vecCurrentSize.push_back(0);
      vecTotalSize.push_back(0);
    }
  }

  connect(comboParameter, SIGNAL(currentIndexChanged(int)), this, SLOT(parameterChanged(int)));

  layCombo->addWidget(labelCombo);
  layCombo->addWidget(comboParameter);

  plot = new QCustomPlot;
  plot->xAxis->setLabel("T");
  plot->yAxis->setLabel("χ2");
  plot->xAxis->setRange(50, 700);
  plot->yAxis->setRange(0,  100);

  if (comboParameter->count() > 0) {
    plot->xAxis->setLabel(paramNamesMap[0].c_str());
    plot->xAxis->setRange(fitParams.GetParameter(paramNamesMap[0]).xmin, fitParams.GetParameter(paramNamesMap[0]).xmax);
  }

  if (modelFitInput != NULL) {
    plot->yAxis->setRange(0, 3. * modelFitInput->Chi2());
  }

  plot->addGraph();
  plot->graph(0)->setName("Calculated");
  plot->graph(0)->setPen(QPen(Qt::black, 2, Qt::SolidLine));
  plot->graph(0)->setLineStyle(QCPGraph::lsLine);
  plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 9.));

  if (modelFitInput != NULL) {
    plot->addGraph();
    plot->graph(1)->setName("Parabolic");
    plot->graph(1)->setPen(QPen(Qt::red, 2, Qt::DashLine));

    plot->addGraph();
    plot->graph(2)->setPen(QPen(Qt::red, 2));
    plot->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 12.));

    plot->legend->removeAt(2);
  }

  
  plot->legend->setVisible(true);

  QVBoxLayout *layInterface = new QVBoxLayout();
  layInterface->setAlignment(Qt::AlignTop);

  QGridLayout *layParams = new QGridLayout();

  labelAmin = new QLabel(comboParameter->currentText() + "<sub>min</sub>:");
  spinAmin = new QDoubleSpinBox();
  spinAmin->setDecimals(5);
  spinAmin->setMinimum(-2000.);
  spinAmin->setMaximum(2000.);
  spinAmin->setValue(50.);
  labelAmax = new QLabel(comboParameter->currentText() + "<sub>max</sub>:");
  spinAmax = new QDoubleSpinBox();
  spinAmax->setDecimals(5);
  spinAmax->setMinimum(-2000.);
  spinAmax->setMaximum(2000.);
  spinAmax->setValue(180.);
  labelAiter = new QLabel(comboParameter->currentText() + "<sub>iters</sub>:");
  spinAiter = new QSpinBox();
  spinAiter->setMinimum(1);
  spinAiter->setMaximum(100000);
  spinAiter->setValue(10);

  if (comboParameter->count() > 0) {
    spinAmin->setValue(vecAleft[0]);
    spinAmax->setValue(vecAright[0]);
  }

  connect(spinAmin, SIGNAL(valueChanged(double)), this, SLOT(limitsChanged()));
  connect(spinAmax, SIGNAL(valueChanged(double)), this, SLOT(limitsChanged()));

  


  layParams->addWidget(labelAmin, 0, 0);
  layParams->addWidget(spinAmin, 0, 1);
  layParams->addWidget(labelAmax, 0, 2);
  layParams->addWidget(spinAmax, 0, 3);
  layParams->addWidget(labelAiter, 0, 4);
  layParams->addWidget(spinAiter, 0, 5);

  QHBoxLayout *laychi2Range = new QHBoxLayout();
  laychi2Range->setAlignment(Qt::AlignLeft);
  QLabel *labelChi2min = new QLabel(tr("χ<sup>2</sup><sub>min</sub>:"));
  spinchi2min = new QDoubleSpinBox();
  spinchi2min->setMinimum(0.);
  spinchi2min->setMaximum(100000.);
  spinchi2min->setValue(0.);
  QLabel *labelChi2max = new QLabel(tr("χ<sup>2</sup><sub>max</sub>:"));
  spinchi2max = new QDoubleSpinBox();
  spinchi2max->setMinimum(0.);
  spinchi2max->setMaximum(100000.);
  spinchi2max->setValue(50.);

  if (modelFitInput != NULL) {
    spinchi2max->setValue(3.*modelFitInput->Chi2());
  }

  laychi2Range->addWidget(labelChi2min);
  laychi2Range->addWidget(spinchi2min);
  laychi2Range->addWidget(labelChi2max);
  laychi2Range->addWidget(spinchi2max);


  QHBoxLayout *layButtons = new QHBoxLayout();
  layButtons->setAlignment(Qt::AlignLeft);
  buttonCalculate = new QPushButton(tr("Calculate χ2 profile"));
  connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));
  buttonReplot = new QPushButton(tr("Replot"));
  connect(buttonReplot, SIGNAL(clicked()), this, SLOT(replot()));
  buttonToFile = new QPushButton(tr("Write to file..."));
  connect(buttonToFile, SIGNAL(clicked()), this, SLOT(writetoFile()));
  layButtons->addWidget(buttonCalculate);
  layButtons->addWidget(buttonReplot);
  layButtons->addWidget(buttonToFile);

  progBar = new QProgressBar();
  progBar->setVisible(false);

  layInterface->addLayout(layCombo);
  layInterface->addLayout(layParams);
  layInterface->addLayout(laychi2Range);
  layInterface->addLayout(layButtons);
  layInterface->addWidget(progBar);


  layout->addWidget(plot, 1);
  layout->addLayout(layInterface);

  parameterChanged(comboParameter->currentIndex());

  setLayout(layout);

  fRunning = false;

  calcTimer = new QTimer(this);
  connect(calcTimer, SIGNAL(timeout()), this, SLOT(updateProgress()));

  fStop = 0;

  setWindowTitle(tr("chi2 profile"));

  replot();
}

void chi2ProfileDialog::setModel()
{
  if (model != NULL)
    delete model;

  if (config.ModelType == ThermalModelConfig::DiagonalEV)
    model = new ThermalModelEVDiagonal(TPS);
  else if (config.ModelType == ThermalModelConfig::CrosstermsEV)
    model = new ThermalModelEVCrossterms(TPS);
  else if (config.ModelType == ThermalModelConfig::CE) {
    model = new ThermalModelCanonical(TPS);
  }
  else if (config.ModelType == ThermalModelConfig::SCE)
    model = new ThermalModelCanonicalStrangeness(TPS);
  else if (config.ModelType == ThermalModelConfig::EVSCE)
    model = new ThermalModelEVCanonicalStrangeness(TPS);
  else if (config.ModelType == ThermalModelConfig::VDWSCE)
    model = new ThermalModelVDWCanonicalStrangeness(TPS);
  else if (config.ModelType == ThermalModelConfig::QvdW)
    model = new ThermalModelVDW(TPS);
  else if (config.ModelType == ThermalModelConfig::CCE)
    model = new ThermalModelCanonicalCharm(TPS);
  else
    model = new ThermalModelIdeal(TPS);

  if (model != NULL)
    model->SetNormBratio(config.RenormalizeBR);

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

  model->SetUseWidth(config.FiniteWidth != 0);
  if (config.FiniteWidth == 1)
    model->SetResonanceWidthIntegrationType(ThermalParticle::BWTwoGamma);
  else if (config.FiniteWidth == 2)
    model->SetResonanceWidthIntegrationType(ThermalParticle::eBW);

  //model->SetQoverB(config.QoverB);
  //model->ConstrainMuQ(config.ConstrainMuQ);
  //model->ConstrainMuS(config.ConstrainMuS);
  //model->ConstrainMuC(config.ConstrainMuC);

  //model->SetNormBratio(config.RenormalizeBR);

  //if (config.QuantumStatisticsType)
  //  model->SetCalculationType(IdealGasFunctions::Quadratures);
  //else
  //  model->SetCalculationType(IdealGasFunctions::ClusterExpansion);

  //model->SetStatistics(config.QuantumStatistics);

  //if (config.QuantumStatisticsInclude == 1 || config.QuantumStatisticsInclude == 2) {
  //  for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
  //    ThermalParticle &part = model->TPS()->Particle(i);
  //    if (config.QuantumStatisticsInclude == 2) {
  //      if (part.PdgId() != 211 && part.PdgId() != 111 && part.PdgId() != -211) {
  //        part.UseStatistics(false);
  //      }
  //    }
  //    else if (config.QuantumStatisticsInclude == 1) {
  //      if (part.BaryonCharge() != 0) {
  //        part.UseStatistics(false);
  //      }
  //    }
  //  }
  //}

  SetThermalModelConfiguration(model, config);

  SetThermalModelInteraction(model, config);

  if (config.ModelType == ThermalModelConfig::CE) {
    static_cast<ThermalModelCanonical*>(model)->CalculateQuantumNumbersRange(false);
  }
}

void chi2ProfileDialog::calculate() {
  if (!fRunning) {
      plot->xAxis->setLabel(comboParameter->currentText());
      plot->xAxis->setRange(spinAmin->value(),spinAmax->value());

      int cindex = comboParameter->currentIndex();
      if (cindex<0 || cindex>=vecParams.size())
        return;

      int &fCurrentSize = vecCurrentSize[cindex];
      int &fTotalSize   = vecTotalSize[cindex];
      std::vector<double> &Avalues = vecAvalues[cindex];
      std::vector<double> &params  = vecParams[cindex];

      Avalues.resize(0);
      double da = (spinAmax->value() - spinAmin->value()) / spinAiter->value();
      for (double a = spinAmin->value(); a <= spinAmax->value() + 1.e-8; a += da) {
              Avalues.push_back(a);
              params.push_back(0.);
          }

      fCurrentSize = 0;
      //fPreviousSize = 0;
      fTotalSize = spinAiter->value();
      fStop = 0;

      replot();

      setModel();

      modelFit = new ThermalModelFit(model);
      modelFit->SetParameters(fitParams);

      modelFit->FixVcOverV(modelFitInput->FixVcOverV());
      modelFit->SetVcOverV(modelFitInput->VcOverV());

      for (int i = 0; i<quantities.size(); ++i) {
        modelFit->AddQuantity(quantities[i]);
      }

      chi2ProfileWorker *wrk = new chi2ProfileWorker(modelFit,comboParameter->currentText().toStdString(), &Avalues, &params,
                                                      &fCurrentSize, &fStop, this);

      connect(wrk, SIGNAL(calculated()), this, SLOT(finalize()));
      connect(wrk, SIGNAL(finished()), wrk, SLOT(deleteLater()));

      wrk->start();

      comboParameter->setEnabled(false);
      buttonCalculate->setText(tr("Stop"));
      buttonReplot->setEnabled(false);
      buttonToFile->setEnabled(false);
      fRunning = true;
      progBar->setVisible(true);

      calcTimer->start(10);
  }
  else {
      fStop = 1;
  }
}

void chi2ProfileDialog::replot() {
  plot->xAxis->setRange(spinAmin->value(),    spinAmax->value());
  plot->yAxis->setRange(spinchi2min->value(), spinchi2max->value());

  int cindex = comboParameter->currentIndex();
  if (cindex<0 || cindex>vecParams.size())
    return;

  plot->graph(0)->data()->clear();

  int fCurrentSize = vecCurrentSize[cindex];
  std::vector<double> &Avalues = vecAvalues[cindex];
  std::vector<double> &params  = vecParams[cindex];

  for (int i = 0; i<fCurrentSize; ++i) {
    plot->graph(0)->addData(Avalues[i], params[i]);
  }

  if (modelFitInput != NULL) {
    std::string ParameterName = comboParameter->currentText().toStdString();

    double mn = 1.;
    if (ParameterName == "T" || ParameterName == "muB" ||
      ParameterName == "muQ" || ParameterName == "muS" ||
      ParameterName == "muC")
      mn = 1.e-3; // MeV to GeV

    double chi2min = modelFitInput->Chi2();
    double amin = modelFitInput->Parameters().GetParameter(comboParameter->currentText().toStdString()).value;
    double Aco = 1 / pow(modelFitInput->Parameters().GetParameter(comboParameter->currentText().toStdString()).error, 2);

    


    plot->graph(1)->data()->clear();

    if (!(chi2min != chi2min || Aco != Aco)) {

      int iters = 200;
      double da = (spinAmax->value() - spinAmin->value()) / iters;
      for (double a = spinAmin->value(); a <= spinAmax->value() + 1.e-8; a += da) {
        double ta = a * mn;
        plot->graph(1)->addData(a, chi2min + Aco * (ta - amin) * (ta - amin));
      }
    }

    plot->graph(2)->data()->clear();
    if (!(chi2min != chi2min))
      plot->graph(2)->addData(amin / mn, chi2min);
  }
  plot->replot();
}

void chi2ProfileDialog::updateProgress() {
  int cindex = comboParameter->currentIndex();
  if (cindex<0 || cindex>vecParams.size())
    return;
  progBar->setRange(0, vecTotalSize[cindex]);
  progBar->setValue(vecCurrentSize[cindex]);
  replot();
}

void chi2ProfileDialog::parameterChanged(int index)
{
  labelAmin->setText(comboParameter->currentText() + "<sub>min</sub>:");
  labelAmax->setText(comboParameter->currentText() + "<sub>max</sub>:");
  labelAiter->setText(comboParameter->currentText() + "<sub>iters</sub>:");
  if (comboParameter->count() > 0 && index != -1) {
    double mn = 1.;
    std::string ParameterName = comboParameter->currentText().toStdString();
    if (ParameterName == "T" || ParameterName == "muB" ||
      ParameterName == "muQ" || ParameterName == "muS" ||
      ParameterName == "muC")
      mn = 1.e3;
    
    spinAmin->blockSignals(true);
    spinAmax->blockSignals(true);
    spinAmin->setValue(vecAleft[index]);
    spinAmax->setValue(vecAright[index]);
    spinAmin->blockSignals(false);
    spinAmax->blockSignals(false);

    plot->xAxis->setLabel(comboParameter->currentText());
    replot();
  }
}

void chi2ProfileDialog::limitsChanged()
{
  int cindex = comboParameter->currentIndex();
  if (cindex<0 || cindex>vecParams.size())
    return;

  vecAleft[cindex]  = spinAmin->value();
  vecAright[cindex] = spinAmax->value();
}


void chi2ProfileDialog::finalize() {
    calcTimer->stop();

    if (modelFit != NULL)
      delete modelFit;
    modelFit = NULL;

    if (model != NULL)
      delete model;
    model = NULL;

    buttonCalculate->setText(tr("Calculate"));
    comboParameter->setEnabled(true);
    buttonReplot->setEnabled(true);
    buttonToFile->setEnabled(true);
    fRunning = false;
    progBar->setVisible(false);
    replot();
}

void chi2ProfileDialog::writetoFile() {
  if (!fRunning) {

    int cindex = comboParameter->currentIndex();
    if (cindex<0 || cindex>vecParams.size())
      return;

    if (lastFilePath == "")
      lastFilePath = QApplication::applicationDirPath() + "/chi2profile.out";
    QString path = QFileDialog::getSaveFileName(this, tr("Save chi2 profile to file"), lastFilePath, tr("*.out"));
    if (path.length()>0)
    {
      lastFilePath = path;
      std::ofstream fout(path.toStdString());

      fout << std::setw(20) << comboParameter->currentText().toStdString()
           << std::setw(20) << "chi2" << std::endl;

      int fCurrentSize = vecCurrentSize[cindex];
      std::vector<double> &Avalues = vecAvalues[cindex];
      std::vector<double> &params = vecParams[cindex];
      for (int i = 0; i<fCurrentSize; ++i) {
        fout << std::setw(20) << Avalues[i]
             << std::setw(20) << params[i] << std::endl;
      }

      fout.close();
    }
  }
}

