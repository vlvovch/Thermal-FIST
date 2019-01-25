/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "SpectralFunctionDialog.h"

#include <QLayout>
#include <QLabel>
#include <QDebug>

#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "ThermalFISTConfig.h"
#include "HRGBase/NumericalIntegration.h"

#include "qcustomplot.h"

using namespace thermalfist;

SpectralFunctionDialog::SpectralFunctionDialog(QWidget *parent, ThermalParticle *part, double Temp, double Chem, int scheme) :
    QDialog(parent), particle(part), T(Temp), Mu(Chem), WidthScheme(scheme)
{
  QHBoxLayout *layout = new QHBoxLayout();

  QVBoxLayout *layoutL = new QVBoxLayout();

  QHBoxLayout *layCombo = new QHBoxLayout();

  QLabel *labelCombo = new QLabel(tr("Scheme:"));
  comboScheme = new QComboBox();
  comboScheme->addItem(QString("Breit-Wigner"));
  comboScheme->addItem(QString("eBW"));

  connect(comboScheme, SIGNAL(currentIndexChanged(int)), this, SLOT(calculate()));

  layCombo->addWidget(labelCombo);
  layCombo->addWidget(comboScheme);

  comboView = new QComboBox();
  comboView->addItem(QString("Spectral function"));
  comboView->addItem(QString("Branching ratios"));
  comboView->addItem(QString("Partial widths"));
  comboView->setCurrentIndex(0);
  connect(comboView, SIGNAL(currentIndexChanged(int)), this, SLOT(replot()));

  QHBoxLayout *layChannel = new QHBoxLayout();
  layChannel->setAlignment(Qt::AlignLeft);

  QLabel *labelChannel = new QLabel(tr("Channel:"));
  comboChannel = new QComboBox();
  comboChannel->addItem(QString("Full"));
  for (int i = 0; i < part->Decays().size(); ++i) {
    comboChannel->addItem(QString(part->Decays()[i].mChannelName.c_str()));
  }
  comboChannel->setCurrentIndex(0);
  connect(comboChannel, SIGNAL(currentIndexChanged(int)), this, SLOT(replot()));

  layChannel->addWidget(labelChannel);
  layChannel->addWidget(comboChannel);


  plot = new QCustomPlot;
  plot->axisRect()->setupFullAxesBox();

  layoutL->addWidget(comboView, 0, Qt::AlignLeft);
  layoutL->addLayout(layChannel);
  layoutL->addWidget(plot, 1);

  QVBoxLayout *layInterface = new QVBoxLayout();
  layInterface->setAlignment(Qt::AlignTop);

  QGridLayout *layParams = new QGridLayout();

  QLabel *labelT = new QLabel("T [MeV]:");
  spinT = new QDoubleSpinBox();
  spinT->setDecimals(5);
  spinT->setMinimum(20.);
  spinT->setMaximum(2000.);
  spinT->setValue(T * 1.e3);
  connect(spinT, SIGNAL(editingFinished()), this, SLOT(calculate()));

  QLabel *labelMu = new QLabel("Mu [MeV]:");
  spinMu = new QDoubleSpinBox();
  spinMu->setDecimals(5);
  spinMu->setMinimum(-2000.);
  spinMu->setMaximum(2000.);
  spinMu->setValue(Mu * 1.e3);
  connect(spinMu, SIGNAL(editingFinished()), this, SLOT(calculate()));


  layParams->addWidget(labelT, 0, 0);
  layParams->addWidget(spinT, 0, 1);
  
  layParams->addWidget(labelMu, 1, 0);
  layParams->addWidget(spinMu, 1, 1);
  
  layInterface->addLayout(layCombo);
  layInterface->addLayout(layParams);


  layout->addLayout(layoutL, 1);
  layout->addLayout(layInterface);

  setLayout(layout);

  setWindowTitle(QString(particle->Name().c_str()) + tr(" spectral function"));

  calculate();
}

void SpectralFunctionDialog::replot()
{
  if (comboView->currentIndex() == 0)
    comboChannel->setEnabled(true);
  else
    comboChannel->setEnabled(false);
  
  if (comboView->currentIndex() == 0) replotSpectralFunctions();
  else if (comboView->currentIndex() == 1)  replotBranchingRatios();
  else replotPartialWidths();
}

void SpectralFunctionDialog::replotSpectralFunctions()
{
  plot->xAxis->setLabel("M [GeV]");
  plot->xAxis->setLabelFont(QFont("Arial", font().pointSize() + 4));
  plot->yAxis->setLabel("rho(M) [GeV-1]");
  plot->yAxis->setLabelFont(QFont("Arial", font().pointSize() + 4));



  plot->yAxis->setScaleType(QCPAxis::stLinear);
  plot->yAxis->setTicker(QSharedPointer<QCPAxisTicker>(new QCPAxisTicker));

  plot->clearGraphs();
  plot->clearItems();

  QCPItemStraightLine *item = new QCPItemStraightLine(plot);
  item->setPen(QPen(Qt::black, 1, Qt::DashLine));
  item->point1->setCoords(particle->Mass(), 0.);
  item->point2->setCoords(particle->Mass(), 1.);
  //plot->addItem(item);

  plot->addGraph();
  plot->graph(0)->setName("Vacuum");
  plot->graph(0)->setPen(QPen(Qt::black, 2, Qt::DashLine));
  plot->graph(0)->setLineStyle(QCPGraph::lsLine);

  plot->addGraph();
  plot->graph(1)->setName("Thermal");
  plot->graph(1)->setPen(QPen(Qt::red, 2, Qt::SolidLine));
  plot->graph(1)->setBrush(QBrush(QColor(255, 0, 0, 10)));
  plot->graph(1)->setLineStyle(QCPGraph::lsLine);

  plot->legend->setFont(QFont("Arial", font().pointSize() + 2));
  plot->legend->setVisible(true);
  
  double maxval = 0.;
  for (int ind = 0; ind < BWm.size(); ++ind) {
    maxval = std::max(maxval, BWval[ind]);
    maxval = std::max(maxval, BWTHval[ind]);
  }

  plot->xAxis->setRange(xleft, xright);
  plot->yAxis->setRange(0., 1.1*maxval);

  //// Fill static and thermal
  plot->graph(0)->data()->clear();
  plot->graph(1)->data()->clear();
  for (int i = 0; i<BWm.size(); ++i) {
    double mnozh = 1.;
    if (comboChannel->currentIndex() != 0) {
      if (WidthScheme == 0)
        mnozh = particle->Decays()[comboChannel->currentIndex() - 1].mBratio;
      else
        mnozh = particle->Decays()[comboChannel->currentIndex() - 1].ModifiedWidth(BWm[i]) * particle->ResonanceWidth() / particle->TotalWidtheBW(BWm[i]);
    }

    plot->graph(0)->addData(BWm[i], mnozh * BWval[i]);
    plot->graph(1)->addData(BWm[i], mnozh * BWTHval[i]);
  }

  

  plot->replot();
}

void SpectralFunctionDialog::replotBranchingRatios()
{
  plot->xAxis->setLabel("M [GeV]");
  plot->yAxis->setLabel("BR");

  plot->yAxis->setScaleType(QCPAxis::stLinear);
  plot->yAxis->setTicker(QSharedPointer<QCPAxisTicker>(new QCPAxisTicker));

  plot->xAxis->setRange(xleft, xright);
  plot->yAxis->setRange(0., 1.);

  plot->clearGraphs();
  plot->clearItems();

  QCPItemStraightLine *item = new QCPItemStraightLine(plot);
  item->setPen(QPen(Qt::black, 1, Qt::DashLine));
  item->point1->setCoords(particle->Mass(), 0.);
  item->point2->setCoords(particle->Mass(), 1.);
  //plot->addItem(item);

  QVector<QColor> colors;
  colors.push_back(QColor(Qt::black));
  colors.push_back(QColor(Qt::red));
  colors.push_back(QColor(Qt::blue));
  colors.push_back(QColor(Qt::green));

  QVector<Qt::PenStyle> styles;
  styles.push_back(Qt::SolidLine);
  styles.push_back(Qt::DashLine);
  styles.push_back(Qt::DashDotLine);
  

  for (int i = 0; i < particle->Decays().size(); ++i) {
    plot->addGraph();

    QColor tcolor = colors[i % 4];
    Qt::PenStyle penstyle = styles[ (i/4)%3 ];

    plot->graph(i)->setName(QString(particle->Decays()[i].mChannelName.c_str()));
    plot->graph(i)->setPen(QPen(tcolor, 2, penstyle));
    plot->graph(i)->setLineStyle(QCPGraph::lsLine);
  }

  for (int ind = 0; ind < particle->Decays().size(); ++ind) {
    plot->graph(ind)->data()->clear();
    for (int i = 0; i < BWm.size(); ++i) {
      double BR = 1.;
      if (WidthScheme == 0)
        BR = particle->Decays()[ind].mBratio;
      else
        BR = particle->Decays()[ind].ModifiedWidth(BWm[i]) * particle->ResonanceWidth() / particle->TotalWidtheBW(BWm[i]);

      plot->graph(ind)->addData(BWm[i], BR);
    }
  }

  plot->replot();
}

void SpectralFunctionDialog::replotPartialWidths()
{
  plot->xAxis->setLabel("M [GeV]");
  plot->yAxis->setLabel("Gamma_i [GeV]");

  plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
  plot->yAxis->setTicker(QSharedPointer<QCPAxisTickerLog>(new QCPAxisTickerLog));

  plot->xAxis->setRange(xleft, xright);

  double MaxWidth = 1.1 * particle->ResonanceWidth();
  if (WidthScheme != 0)
    MaxWidth = 1.1 * particle->TotalWidtheBW(BWm[BWm.size() - 1]);

  plot->yAxis->setRange(1.e-3, MaxWidth);

  plot->clearGraphs();
  plot->clearItems();

  QCPItemStraightLine *item = new QCPItemStraightLine(plot);
  item->setPen(QPen(Qt::black, 1, Qt::DashLine));
  item->point1->setCoords(particle->Mass(), 0.);
  item->point2->setCoords(particle->Mass(), 1.);
  //plot->addItem(item);

  QCPItemStraightLine *itemGamma = new QCPItemStraightLine(plot);
  itemGamma->setPen(QPen(Qt::black, 1, Qt::DashLine));
  itemGamma->point1->setCoords(0., particle->ResonanceWidth());
  itemGamma->point2->setCoords(1., particle->ResonanceWidth());
  //plot->addItem(itemGamma);

  QVector<QColor> colors;
  colors.push_back(QColor(Qt::black));
  colors.push_back(QColor(Qt::red));
  colors.push_back(QColor(Qt::blue));
  colors.push_back(QColor(Qt::green));

  QVector<Qt::PenStyle> styles;
  styles.push_back(Qt::SolidLine);
  styles.push_back(Qt::DashLine);
  styles.push_back(Qt::DashDotLine);


  plot->addGraph();
  plot->graph(0)->setName("Total");
  plot->graph(0)->setPen(QPen(Qt::black, 4, Qt::SolidLine));
  plot->graph(0)->setLineStyle(QCPGraph::lsLine);

  for (int i = 0; i < particle->Decays().size(); ++i) {
    plot->addGraph();

    QColor tcolor = colors[i % 4];
    Qt::PenStyle penstyle = styles[(i / 4) % 3];

    plot->graph(i+1)->setName(QString(particle->Decays()[i].mChannelName.c_str()));
    plot->graph(i+1)->setPen(QPen(tcolor, 2, penstyle));
    plot->graph(i+1)->setLineStyle(QCPGraph::lsLine);
  }


  for (int i = 0; i < BWm.size(); ++i) {
    double BR = 1.;
    if (WidthScheme == 0)
      BR = particle->ResonanceWidth();
    else
      BR = particle->TotalWidtheBW(BWm[i]);

    plot->graph(0)->addData(BWm[i], BR);
  }


  for (int ind = 0; ind < particle->Decays().size(); ++ind) {
    plot->graph(ind + 1)->data()->clear();
    for (int i = 0; i < BWm.size(); ++i) {
      double BR = 1.;
      if (WidthScheme == 0)
        BR = particle->Decays()[ind].mBratio * particle->ResonanceWidth();
      else
        BR = particle->Decays()[ind].ModifiedWidth(BWm[i]) * particle->ResonanceWidth();

      plot->graph(ind + 1)->addData(BWm[i], BR);
    }
  }

  plot->replot();
}

void SpectralFunctionDialog::calculate()
{
  WidthScheme = comboScheme->currentIndex();
  T = spinT->value() * 1.e-3;
  Mu = spinMu->value() * 1.e-3;

  if (WidthScheme == 0) // BWTwoGamma scheme
  {
    double Threshold = particle->DecayThresholdMass();
    double Width     = particle->ResonanceWidth();
    double Mass      = particle->Mass();

    double a = std::max(Threshold, Mass - 2.*Width);
    double b = Mass + 2.*Width;


    // Normalization
    std::vector<double> xlegBW, wlegBW;

    NumericalIntegration::GetCoefsIntegrateLegendre10(a, b, &xlegBW, &wlegBW);

    norm = 0.; 
    normTH = 0.;

    for (int i = 0; i < xlegBW.size(); ++i) {
      double M = xlegBW[i];
      norm += wlegBW[i] * particle->MassDistribution(M, Width);
      normTH += wlegBW[i] * particle->ThermalMassDistribution(M, T, Mu, Width);
    }


    BWm.resize(0);
    BWval.resize(0);
    BWTHval.resize(0);
    int iters = 200;
    double dM = (b - a) / (iters - 1);
    for (int ind = 0; ind < iters; ++ind) {
      double M = a + dM * ind;

      double ty = particle->MassDistribution(M, Width) / norm;
      BWm.push_back(M);
      BWval.push_back(ty);

      ty = particle->ThermalMassDistribution(M, T, Mu, Width) / normTH;
      BWTHval.push_back(ty);
    }

    // Set ranges
    xleft = std::min(a, particle->DecayThresholdMassDynamical());
    xright = b + 0.2 * (b - xleft);
  }
  else // eBW
  {
    double Threshold = particle->DecayThresholdMass();
    double Width = particle->ResonanceWidth();
    double Mass = particle->Mass();

    double a = std::max(Threshold, Mass - 2.*Width);
    double b = Mass + 2.*Width;


    // Normalization
    std::vector<double> xlegpdyn, wlegpdyn; // [Thr, a]
    std::vector<double> xlegdyn, wlegdyn;   // [a, b]
    std::vector<double> xlagdyn, wlagdyn;   // [b, inf]
    std::vector<double> xall, wall;         // All

    if (particle->Decays().size() == 0)
      a = Mass - 2.*Width + 1.e-6;
    else
      a = particle->DecayThresholdMassDynamical();


    b = Mass + 2.*Width;
    if (a >= Mass - 2.*Width) {
      xlegpdyn.resize(0);
      if (a >= Mass + 2.*Width)
        xlegdyn.resize(0);
      else
        NumericalIntegration::GetCoefsIntegrateLegendre32(a, b, &xlegdyn, &wlegdyn);
    }
    else {
      NumericalIntegration::GetCoefsIntegrateLegendre32(a, Mass - 2.*Width, &xlegpdyn, &wlegpdyn);
      NumericalIntegration::GetCoefsIntegrateLegendre32(Mass - 2.*Width, b, &xlegdyn, &wlegdyn);
    }
    NumericalIntegration::GetCoefsIntegrateLaguerre32(&xlagdyn, &wlagdyn);

    for (int i = 0; i < xlegpdyn.size(); ++i) {
      xall.push_back(xlegpdyn[i]);
      wall.push_back(wlegpdyn[i]);
    }

    for (int i = 0; i < xlegdyn.size(); ++i) {
      xall.push_back(xlegdyn[i]);
      wall.push_back(wlegdyn[i]);
    }

    for (int i = 0; i < xlagdyn.size(); ++i) {
      xall.push_back(Mass + 2. * Width + xlagdyn[i] * Width);
      wall.push_back(wlagdyn[i] * Width);
    }

    norm = 0.;
    normTH = 0.;

    for (int i = 0; i < xall.size(); ++i) {
      double M = xall[i];
      double MWidth = particle->TotalWidtheBW(M);
      norm += wall[i] * particle->MassDistribution(M, MWidth);
      normTH += wall[i] * particle->ThermalMassDistribution(M, T, Mu, MWidth);
    }

    xleft = xall[0];
    xright = b + 0.2 * (b - xleft);

    a = xleft;
    b = xright;

    BWm.resize(0);
    BWval.resize(0);
    BWTHval.resize(0);
    // Compute all values
    int iters = 400;
    double dM = (b - a) / (iters - 1);
    for (int ind = 0; ind < iters; ++ind) {
      double M = a + dM * ind;
      double MWidth = particle->TotalWidtheBW(M);

      double ty = particle->MassDistribution(M, MWidth) / norm;
      BWm.push_back(M);
      BWval.push_back(ty);

      ty = particle->ThermalMassDistribution(M, T, Mu, MWidth) / normTH;
      BWTHval.push_back(ty);
    }
  }

  replot();
}


