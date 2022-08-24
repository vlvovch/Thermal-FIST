/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "eventgeneratortab.h"

#include <algorithm>

#include <QLayout>
//#include <QFileDialog>
#include <QLabel>
#include <QHeaderView>
#include <QGroupBox>
#include <QApplication>
#include <QMessageBox>
#include <QElapsedTimer>
#include <QDebug>
#include <QFileDialog>
#include <QTimer>
#include <QSplitter>
#include <QStackedWidget>

#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGEV/ThermalModelEVCrossterms.h"
#include "HRGEV/ExcludedVolumeHelper.h"
#include "HRGVDW/ThermalModelVDW.h"
#include "HRGBase/ThermalModelCanonical.h"
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonicalCharm.h"
#include "HRGEventGenerator/SphericalBlastWaveEventGenerator.h"
#include "HRGEventGenerator/CylindricalBlastWaveEventGenerator.h"
#include "HRGEventGenerator/CracowFreezeoutEventGenerator.h"

//#include "EventGeneratorExtensions.h"

#include "DebugText.h"
#include "spectramodel.h"
#include "particlespectra.h"
#include "qcustomplot.h"

using namespace thermalfist;

EventGeneratorWorker::EventGeneratorWorker(
  thermalfist::EventGeneratorBase* gen,
  ParticlesSpectra* spec,
  QMutex* mut,
  int totalEvents,
  int* evproc,
  int* stopo,
  double* nEp,
  bool pDecays,
  std::string fileout,
  QObject* parent):
      QThread(parent), generator(gen), spectra(spec), mutex(mut),
          events(totalEvents), eventsProcessed(evproc), stop(stopo), nE(nEp), performDecays(pDecays)/*,
          hepmcout(fileout + ".hepmc3")*/
  {
      wsum = w2sum = 0.;
      fout.clear();
      if (fileout != "") 
        fout.open(fileout.c_str());

  }

void EventGeneratorWorker::run()
{
     if (mutex!=NULL) {
        SimpleEvent::EventOutputConfig outconfig;
        //outconfig.printEnergy = true;
        //outconfig.printDecayEpoch = true;
        outconfig.printMotherPdg = true;
        outconfig.printPhotonsLeptons = true;

        for(int i=0;i<events && !(*stop);++i) {
            SimpleEvent ev = generator->GetEvent(performDecays);
            mutex->lock();
            wsum += ev.weight;
            w2sum += ev.weight*ev.weight;
            spectra->ProcessEvent(ev);
            (*eventsProcessed)++;

            if (fout.is_open())
              ev.writeToFile(fout, outconfig, (*eventsProcessed));

            //hepmcout.WriteEvent(ev);

            mutex->unlock();
        }
     }
     if (fout.is_open())
        fout.close();
     *nE = wsum * wsum / w2sum;
     emit calculated();
 }

BinningDialog::BinningDialog(BinningConfig* config, QWidget* parent) : QDialog(parent), m_config(config)
{
  if (config == 0) {
    QMessageBox msgBox;
    msgBox.setText(tr("Empty binning config!"));
    msgBox.exec();
    QDialog::reject();
  }
  
  QVBoxLayout* layout = new QVBoxLayout();

  QGroupBox* gr1D = new QGroupBox(tr("1D plots"));
  QHBoxLayout* grLayout = new QHBoxLayout();
  grLayout->addWidget(new QLabel(tr("Number of bins:")));
  spinBins = new QSpinBox();
  spinBins->setMinimum(1);
  spinBins->setMaximum(10000);
  spinBins->setValue(m_config->bins1D);
  grLayout->addWidget(spinBins);
  gr1D->setLayout(grLayout);

  QGroupBox* gr2D = new QGroupBox(tr("2D plots"));
  QGridLayout* gr2Layout = new QGridLayout();
  gr2Layout->addWidget(new QLabel(tr("Number of X bins:")), 0, 0);
  spinBinsX = new QSpinBox();
  spinBinsX->setMinimum(1);
  spinBinsX->setMaximum(10000);
  spinBinsX->setValue(m_config->bins2D_x);
  gr2Layout->addWidget(spinBinsX, 0, 1);
  gr2Layout->addWidget(new QLabel(tr("Number of Y bins:")), 1, 0);
  spinBinsY = new QSpinBox();
  spinBinsY->setMinimum(1);
  spinBinsY->setMaximum(10000);
  spinBinsY->setValue(m_config->bins2D_y);
  gr2Layout->addWidget(spinBinsY, 1, 1);
  gr2D->setLayout(gr2Layout);

  std::cout << m_config->bins1D << " " << m_config->bins2D_x << " " << m_config->bins2D_y << std::endl;

  layout->addWidget(gr1D);
  layout->addWidget(gr2D);
  setLayout(layout);

  QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok
    | QDialogButtonBox::Cancel);
  connect(buttonBox, &QDialogButtonBox::accepted, this, &BinningDialog::OK);
  connect(buttonBox, &QDialogButtonBox::rejected, this, &BinningDialog::Discard);

  layout->addWidget(buttonBox);

  setWindowTitle(tr("Binning configuration"));
}

void BinningDialog::OK()
{
  m_config->bins1D = spinBins->value();
  m_config->bins2D_x = spinBinsX->value();
  m_config->bins2D_y = spinBinsY->value();
  QDialog::accept();
}



ParticlesAnalyzeDialog::ParticlesAnalyzeDialog(ParticlesAnalyzeConfig* config, QWidget* parent) : QDialog(parent), m_config(config)
{
  if (config == 0) {
    QMessageBox msgBox;
    msgBox.setText(tr("Empty particle list config!"));
    msgBox.exec();
    QDialog::reject();
  }

  QVBoxLayout* layout = new QVBoxLayout();

  RBAll = new QRadioButton(tr("All hadrons and resonances"));
  RBStable = new QRadioButton(tr("Only stable"));
  RBStablePlus = new QRadioButton(tr("Only stable plus specified (enter comma-separated list of pdg codes)"));
  if (m_config->type == 0)
    RBAll->setChecked(true);
  else if (m_config->type == 1)
    RBStable->setChecked(true);
  else
    RBStablePlus->setChecked(true);

  lePdgs = new QLineEdit();
  {
    std::set<long long>::iterator it = m_config->pdgCodes.begin();
    while (it != m_config->pdgCodes.end())
    {
      if (it != m_config->pdgCodes.begin()) {
        lePdgs->setText(lePdgs->text() + ",");
      }
      lePdgs->setText(lePdgs->text() + QString::number(*it));
      it++;
    }
  }

  layout->addWidget(RBAll);
  layout->addWidget(RBStable);
  layout->addWidget(RBStablePlus);
  layout->addWidget(lePdgs);

  QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok
    | QDialogButtonBox::Cancel);
  connect(buttonBox, &QDialogButtonBox::accepted, this, &ParticlesAnalyzeDialog::OK);
  connect(buttonBox, &QDialogButtonBox::rejected, this, &ParticlesAnalyzeDialog::Discard);

  layout->addWidget(buttonBox);

  setLayout(layout);

  setWindowTitle(tr("Particle list configuration"));
}

void ParticlesAnalyzeDialog::OK()
{
  if (RBAll->isChecked())
    m_config->type = 0;
  else if (RBStable->isChecked())
    m_config->type = 1;
  else
    m_config->type = 2;

  std::vector<std::string> pdgs = CuteHRGHelper::split(lePdgs->text().toStdString(), ',');
  m_config->pdgCodes.clear();
  for (int i = 0; i < pdgs.size(); ++i)
    m_config->pdgCodes.insert(std::stoll(pdgs[i]));
  QDialog::accept();
}


EventGeneratorTab::EventGeneratorTab(QWidget *parent, ThermalModelBase *modelop) :
    QWidget(parent)
{
    int index = 0;
    QString tname;

    tname = "dN/dp";
    paramnames.push_back(tname);
    paramnamesx.push_back("p (GeV)");
    parammap[tname] = index;
    index++;

    tname = "dN/dy";
    paramnames.push_back(tname);
    paramnamesx.push_back("y");
    parammap[tname] = index;
    index++;

    tname = "dN/mtdmt";
    paramnames.push_back(tname);
    paramnamesx.push_back("mt-m0 (GeV)");
    parammap[tname] = index;
    index++;

    tname = "d2N/dydpt";
    paramnames.push_back(tname);
    paramnamesx.push_back("y");
    parammap[tname] = index;
    index++;

    tname = "dN/dpt";
    paramnames.push_back(tname);
    paramnamesx.push_back("pT (GeV)");
    parammap[tname] = index;
    index++;

    //fitcopy = NULL;
    dbgstrm.setString(&dbgstr);

    TPS = modelop->TPS();
    model = new ThermalModelIdeal(modelop->TPS());
    *model = *modelop;
    spectra = new ParticlesSpectra(model);


    QHBoxLayout *DataEditLay = new QHBoxLayout();

    QVBoxLayout *dataLayv = new QVBoxLayout();
    QVBoxLayout *dataLayv2 = new QVBoxLayout();

    QHBoxLayout* layTop = new QHBoxLayout();
    QLabel *labelParticleList = new QLabel(tr("Particle list:"));
    //checkStableOnly = new QCheckBox(tr("Analyze only stable hadrons"));
    //checkStableOnly->setChecked(true);
    //layTop->addLayout(labelParticleList, 1);
    QPushButton* buttonPartsConfig = new QPushButton(tr("Edit particle list for analysis..."));
    connect(buttonPartsConfig, SIGNAL(clicked()), this, SLOT(changeParticles()));
    layTop->addWidget(labelParticleList, 1);
    layTop->addWidget(buttonPartsConfig);


    myModel = new SpectraModel(this);
    myModel->setSpectra(spectra);
    tableSpectra = new QTableView();
    tableSpectra->setModel(myModel);
    tableSpectra->setSelectionBehavior(QAbstractItemView::SelectRows);
    tableSpectra->setSelectionMode(QAbstractItemView::SingleSelection);
    tableSpectra->resizeColumnsToContents();
    //connect(tableSpectra, SIGNAL(doubleClicked(QModelIndex)), this, SLOT(quantityDoubleClick(QModelIndex)));

    connect(tableSpectra->selectionModel(), SIGNAL(selectionChanged(QItemSelection,QItemSelection)), this, SLOT(changedRow()));

    QHBoxLayout *selLay = new QHBoxLayout();
    selLay->setAlignment(Qt::AlignLeft);
    QLabel *labelSel = new QLabel(tr("Distribution:"));
    comboDistr = new QComboBox();
    for(int i=0;i<index;++i) comboDistr->addItem(paramnames[i]);
    comboDistr->setCurrentIndex(0);
    connect(comboDistr, SIGNAL(currentIndexChanged(int)), this, SLOT(changePlot()));

    QPushButton* buttonBinning = new QPushButton(tr("Binning..."));
    connect(buttonBinning, SIGNAL(clicked()), this, SLOT(changeBinning()));

    selLay->addWidget(labelSel);
    selLay->addWidget(comboDistr);
    selLay->addWidget(buttonBinning, 1, Qt::AlignRight);

    plotDistr = new QCustomPlot();
    plotDistr->xAxis->setLabel("p (GeV/c)");
    plotDistr->yAxis->setLabel(comboDistr->currentText());

    plotDistr->addGraph();
    //plotDistr->graph(0)->setPen(QPen(Qt::blue, 2, Qt::SolidLine));
    plotDistr->graph(0)->setLineStyle(QCPGraph::lsNone);
    plotDistr->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    plotDistr->graph(0)->setName("Monte Carlo");
    //plotDistr->graph(0)->setErrorType(QCPGraph::etValue);
    plotDistr->addGraph();
    plotDistr->graph(1)->setName(comboDistr->currentText());
    plotDistr->graph(1)->setPen(QPen(Qt::blue, 2, Qt::DashLine));

    plotDistr->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(plotDistr, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequestPlot1D(QPoint)));

    errorBars = new QCPErrorBars(plotDistr->xAxis, plotDistr->yAxis);
    errorBars->removeFromLegend();
    errorBars->setErrorType(QCPErrorBars::etValueError);
    errorBars->setPen(QPen(Qt::blue));
    errorBars->setSymbolGap(5);

    plot2D = new QCustomPlot();
    plot2D->xAxis->setLabel(tr("y"));
    plot2D->yAxis->setLabel(tr("pT (GeV/c)"));

    plot2D->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(plot2D, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequestPlot2D(QPoint)));

    colormap = new QCPColorMap(plot2D->xAxis, plot2D->yAxis);

    //plot2D->addPlottable(colormap);

    colorScale = new QCPColorScale(plot2D);
    plot2D->plotLayout()->addElement(0, 1, colorScale);

    plot = new QStackedWidget();
    plot->addWidget(plotDistr);

    plot->addWidget(plot2D);
    if (comboDistr->currentIndex()==3) plot->setCurrentIndex(1);
    else plot->setCurrentIndex(0);

    //dataLayv->addWidget(labelParticleList);
    dataLayv->addLayout(layTop, 1);
    dataLayv->addWidget(tableSpectra,1);
    dataLayv2->addLayout(selLay);
    dataLayv2->addWidget(plot, 1);

    QWidget *container = new QWidget;
    container->setLayout(dataLayv);
    QWidget *container2 = new QWidget;
    container2->setLayout(dataLayv2);

    QSplitter *splitter = new QSplitter();
    splitter->setOrientation(Qt::Vertical);
    splitter->addWidget(container);
    splitter->addWidget(container2);


    QVBoxLayout *editorLay = new QVBoxLayout();
    editorLay->setContentsMargins(15, 0, 0, 0);
    editorLay->setAlignment(Qt::AlignTop);

    QGroupBox *grModelConfig = new QGroupBox(tr("HRG model configuration:"));

    configWidget = new ModelConfigWidget(NULL, model, true);
    connect(configWidget, SIGNAL(changed()), this, SLOT(modelChanged()));

    QHBoxLayout *layModelConfig = new QHBoxLayout();
    layModelConfig->setAlignment(Qt::AlignLeft);
    layModelConfig->addWidget(configWidget);

    grModelConfig->setLayout(layModelConfig);


    QGroupBox *grParameters = new QGroupBox(tr("Chemical freeze-out parameters:"));

    QGridLayout *layParameters = new QGridLayout();
    layParameters->setAlignment(Qt::AlignLeft);
    //layQuant->setHorizontalSpacing(15);
    QLabel *labelTemperature = new QLabel(tr("T<sub>ch</sub> (MeV):"));
    spinTemperature = new QDoubleSpinBox();
    spinTemperature->setMinimum(1.);
    spinTemperature->setMaximum(10000.);
    spinTemperature->setValue(model->Parameters().T * 1e3);
    spinTemperature->setToolTip(tr("Chemical freeze-out temperature (fixes the yields)"));
    connect(spinTemperature, SIGNAL(valueChanged(double)), this, SLOT(changeTkin(double)));
    QLabel *labelmuB = new QLabel(tr("μ<sub>B</sub> (MeV):"));
    spinmuB = new QDoubleSpinBox();
    spinmuB->setMinimum(-1000.);
    spinmuB->setMaximum(1000.);
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

    QLabel *labelVolumeRSC = new QLabel(tr("R<sub>SC</sub> (fm):"));
    spinVolumeRSC = new QDoubleSpinBox();
    spinVolumeRSC->setDecimals(4);
    spinVolumeRSC->setMinimum(0.);
    spinVolumeRSC->setMaximum(25.);
    spinVolumeRSC->setValue(spinVolumeR->value());
    spinVolumeRSC->setEnabled(false);
    spinVolumeRSC->setToolTip(tr("Correlation radius: the (canonical) correlation volume is a sphere of this radius"));


    QLabel* labelVolume    = new QLabel(tr("V (fm<sup>3</sup>):"));
    labelVolumeVal = new QLabel("4000");
    QLabel* labelVolumeSC = new QLabel(tr("V<sub>SC</sub> (fm<sup>3</sup>):"));
    labelVolumeSCVal = new QLabel("4000");

    changeVolumeRSC(spinVolumeRSC->value());

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
    //layParameters->addWidget(labelVolumeSC, 2, 6, 1, 1, Qt::AlignRight);
    //layParameters->addWidget(labelVolumeSCVal, 2, 7);

    layParameters->addWidget(labelB, 3, 0, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinB, 3, 1);
    layParameters->addWidget(labelS, 3, 4, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinS, 3, 5);
    layParameters->addWidget(labelQ, 3, 2, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinQ, 3, 3);
    layParameters->addWidget(labelC, 3, 6, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinC, 3, 7);

    grParameters->setLayout(layParameters);


    QGroupBox *grParametersMom = new QGroupBox(tr("Blast-wave momentum spectrum:"));
    QVBoxLayout *layMom = new QVBoxLayout();

    QHBoxLayout *layMom1 = new QHBoxLayout();
    layMom1->setAlignment(Qt::AlignLeft);
    radioSR = new QRadioButton(tr("Spherically symmetric"));
    radioSR->setToolTip(tr("Siemens-Rasmussen momentum distribution"));
    radioSSH = new QRadioButton(tr("Cylindrically symmetric"));
    radioSSH->setToolTip(tr("Schnedermann-Sollfrank-Heinz prescription"));
    radioCracow = new QRadioButton(tr("Cracow model"));
    radioCracow->setToolTip(tr("Cracow model"));
    radioSR->setChecked(true);
    connect(radioSR, SIGNAL(clicked()), this, SLOT(modelChanged()));
    connect(radioSSH, SIGNAL(clicked()), this, SLOT(modelChanged()));
    connect(radioCracow, SIGNAL(clicked()), this, SLOT(modelChanged()));


    layMom1->addWidget(radioSR);
    layMom1->addWidget(radioSSH);
    layMom1->addWidget(radioCracow);

    QHBoxLayout *layMom2 = new QHBoxLayout();
    layMom2->setAlignment(Qt::AlignLeft);

    QLabel *labelTkin = new QLabel(tr("T<sub>BW</sub> (MeV):"));
    spinTkin = new QDoubleSpinBox();
    spinTkin->setMinimum(1.);
    spinTkin->setMaximum(10000.);
    spinTkin->setValue(model->Parameters().T * 1e3);
    spinTkin->setToolTip(tr("Blast-wave temperature which fixes the momentum spectrum. Set it after T<sub>ch</sub> ans T<sub>kin</sub> (PCE) is fixed."));
    labelBeta = new QLabel(tr("β:"));
    spinBeta = new QDoubleSpinBox();
    spinBeta->setMinimum(0.);
    spinBeta->setMaximum(1.);
    spinBeta->setDecimals(3);
    prevBeta = 0.5;
    prevRmax = 9.0;
    spinBeta->setValue(prevBeta);
    spinBeta->setToolTip(tr("Radial flow velocity"));
    labelBetat = new QLabel(tr("〈β〉<sub>T</sub>:"));
    spinBetat = new QDoubleSpinBox();
    spinBetat->setMinimum(0.);
    spinBetat->setMaximum(1.);
    spinBetat->setDecimals(3);
    spinBetat->setValue(0.5);
    spinBetat->setToolTip(tr("Transverse flow parameter"));
    QLabel *labelEtaMax = new QLabel(tr("η<sub>max</sub>:"));
    spinEtaMax = new QDoubleSpinBox();
    spinEtaMax->setMinimum(0.);
    spinEtaMax->setMaximum(100.);
    spinEtaMax->setDecimals(3);
    spinEtaMax->setValue(1.);
    spinEtaMax->setToolTip(tr("Longitudinal space-time rapidity cut-off |η|<η<sub>max</sub>"));
    QLabel *labeln = new QLabel(tr("n:"));
    spinn = new QDoubleSpinBox();
    spinn->setMinimum(0.);
    spinn->setMaximum(10.);
    spinn->setDecimals(3);
    spinn->setValue(1.);
    spinn->setToolTip(tr("Transverse radial flow profile"));

    layMom2->addWidget(labelTkin);
    layMom2->addWidget(spinTkin);
    layMom2->addSpacing(20);
    layMom2->addWidget(labelBeta);
    layMom2->addWidget(spinBeta);
    layMom2->addSpacing(20);
    layMom2->addWidget(labelBetat);
    layMom2->addWidget(spinBetat);
    layMom2->addSpacing(5);
    layMom2->addWidget(labeln);
    layMom2->addWidget(spinn);
    layMom2->addSpacing(5);
    layMom2->addWidget(labelEtaMax);
    layMom2->addWidget(spinEtaMax);

    layMom->addLayout(layMom1);
    layMom->addLayout(layMom2);

    grParametersMom->setLayout(layMom);

    QHBoxLayout *layFlags = new QHBoxLayout();
    layFlags->setAlignment(Qt::AlignLeft);

    checkDecays = new QCheckBox(tr("Perform decays"));
    checkDecays->setChecked(false);
    checkDecays->setToolTip(tr("Perform chain of decays of all unstable particles"));

    layFlags->addWidget(checkDecays);
    
    
    


    QHBoxLayout *layEvents = new QHBoxLayout();
    layEvents->setAlignment(Qt::AlignLeft);

    checkFile = new QCheckBox(tr("Write events to file"));
    checkFile->setChecked(false);
    checkFile->setToolTip(tr("Writes generated events to file"));
    connect(checkFile, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));

    leFilePath = new QLineEdit("");
    leFilePath->setReadOnly(true);
    leFilePath->setText(QDir::currentPath() + "/events.dat");
    leFilePath->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    buttonChooseFile = new QPushButton(tr("Choose path..."));
    connect(buttonChooseFile, SIGNAL(clicked()), this, SLOT(chooseOutputFile()));
    

    QLabel *labelEvents = new QLabel(tr("Events:"));
    spinEvents = new QSpinBox();
    spinEvents->setMinimum(0);
    spinEvents->setMaximum(1000000000);
    spinEvents->setValue(10000);
    spinEvents->setToolTip(tr("Number of events to generate"));

    layEvents->addWidget(labelEvents);
    layEvents->addWidget(spinEvents);
    layEvents->addStretch(1);
    layEvents->addWidget(checkFile);
    layEvents->addWidget(leFilePath);
    layEvents->addWidget(buttonChooseFile);

    buttonCalculate = new QPushButton(tr("Generate"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));


    progBar = new QProgressBar();
    progBar->setVisible(false);
    progBar->setSizePolicy(QSizePolicy::Preferred,QSizePolicy::Fixed);
    progBar->setAlignment(Qt::AlignCenter);


    teDebug = new QTextEdit;
    //teDebug->setMaximumWidth(620);
    //teDebug->setFixedHeight(200);
    teDebug->setReadOnly(true);

    //editorLay->addLayout(layModelEnsemble);
    //editorLay->addWidget(grEV);
    editorLay->addWidget(grModelConfig);
    editorLay->addWidget(grParameters);
    editorLay->addWidget(grParametersMom);
    editorLay->addLayout(layFlags);
    editorLay->addLayout(layEvents);
    editorLay->addWidget(buttonCalculate, 0, Qt::AlignLeft);
    editorLay->addWidget(progBar);
    editorLay->addWidget(teDebug);

    DataEditLay->addWidget(splitter, 1);
    DataEditLay->addLayout(editorLay, 0);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(DataEditLay);
    setLayout(mainLayout);

    fRunning = false;



    calcTimer = new QTimer(this);
    connect(calcTimer, SIGNAL(timeout()), this, SLOT(updateProgress()));

    modelChanged();

    spinTemperature->setValue(155.0);
    spinmuB->setValue(0.0);
    spinmuQ->setValue(0.0);
    spinmuS->setValue(0.0);
    spinmuC->setValue(0.0);
    spingammaS->setValue(1.0);
    spingammaq->setValue(1.0);
    spinVolumeR->setValue(8.0);

    cpath = QApplication::applicationDirPath();
}



EventGeneratorTab::~EventGeneratorTab()
{
    delete spectra;
    delete model;
}

ThermalModelConfig EventGeneratorTab::getConfig()
{
  ThermalModelConfig ret = configWidget->updatedConfig();

  //ret.QuantumStatistics = 0;
  //ret.QuantumStatisticsType = 0;
  //ret.QuantumStatisticsInclude = 0;


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

  ret.ComputeFluctations = false;

  ret.CanonicalB &= (ret.Ensemble == ThermalModelConfig::EnsembleCE);
  ret.CanonicalQ &= (ret.Ensemble == ThermalModelConfig::EnsembleCE);
  ret.CanonicalS &= (ret.Ensemble == ThermalModelConfig::EnsembleCE || ret.Ensemble == ThermalModelConfig::EnsembleSCE);
  ret.CanonicalC &= (ret.Ensemble == ThermalModelConfig::EnsembleCE || ret.Ensemble == ThermalModelConfig::EnsembleCCE);

  return ret;
}


int EventGeneratorTab::getCurrentRow()
{
    QModelIndexList selectedList = tableSpectra->selectionModel()->selectedRows();
    if (selectedList.count()==0) return -1;
    return selectedList.at(0).row();
}

void EventGeneratorTab::changedRow()
{
    QModelIndexList selectedList = tableSpectra->selectionModel()->selectedRows();
    int id = getCurrentRow();
    if (id<0) id = 0;
    mutex.lock();
    if (id<spectra->fParticles.size()) {
        int ttype = comboDistr->currentIndex();
        if ((ttype>=0 && ttype<=2) || ttype == 4) {
            double tl = 0.;
            if (ttype==1) tl = -3.;
            double tr = 2.;
            if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1)
              tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;
            if (ttype==1) tr = 3.;
            if (ttype==1 && spectra->fDistributionType==1) {
                tl = -3. - spectra->fEtaMax;
                tr = 3. + spectra->fEtaMax;
            }
            spectra->fParticles[id].FillModelDistribution(tl, tr, 100, ttype);
        }
    }
    mutex.unlock();
    if (!fRunning) replot();
}

void EventGeneratorTab::replot() {
    int index = comboDistr->currentIndex();
    int id = getCurrentRow();
    if (id<0) id = 0;
    if (((index>=0 && index<3) || index == 4) && id<spectra->fParticles.size()) {

        double tmin = 1.e5, tmax = 0.;
        double tl = 0., tr = 2.;
        if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1)
          tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;

        plotDistr->graph(0)->data()->clear();
        plotDistr->graph(0)->setName(paramnames[index]);
        plotDistr->graph(1)->data()->clear();
        plotDistr->graph(1)->setName(paramnames[index]);

        int ttype = comboDistr->currentIndex();
        if ((ttype>=0 && ttype<=2) || ttype==4) {
            QVector<double> x1,y1,y1err;
            std::vector<double> tvec = spectra->fParticles[id].GetXVector(ttype);
            x1 = QVector<double>::fromStdVector(tvec);
            //x1    = QVector<double>::fromStdVector (spectra->fParticles[id].GetXVector(ttype) );
            y1    = QVector<double>::fromStdVector (spectra->fParticles[id].GetYVector(ttype) );
            y1err = QVector<double>::fromStdVector (spectra->fParticles[id].GetYErrorVector(ttype) );
            for(int i=0;i<x1.size();++i) {
              if (y1[i] != 0.)
                tmin = std::min(tmin, y1[i]);
              tmax = std::max(tmax, y1[i]);
            }
            plotDistr->graph(0)->setData(x1, y1);
            //plotDistr->graph(0)->setDataValueError(x1,y1,y1err);
            errorBars->data()->clear();
            errorBars->setData(y1err);
            errorBars->setDataPlottable(plotDistr->graph(0));
        }


        plotDistr->xAxis->setLabel(paramnamesx[index]);
        plotDistr->yAxis->setLabel(paramnames[index]);

        QVector<double> x2,y2;
        x2 = QVector<double>::fromStdVector (spectra->fParticles[id].GetModelX() );
        y2 = QVector<double>::fromStdVector (spectra->fParticles[id].GetModelY() );
        for(int i=0;i<x2.size();++i) {
          if (y2[i] != 0.)
            tmin = std::min(tmin, y2[i]);
          tmax = std::max(tmax, y2[i]);
        }
        plotDistr->graph(1)->setData(x2,y2);

        plotDistr->xAxis->setRange(tl, tr);
        if (x2.size()>1) plotDistr->xAxis->setRange(x2[0]-0.5*(x2[1]-x2[0]), x2[x2.size()-1]+0.5*(x2[1]-x2[0]));

        if (comboDistr->currentIndex()==2) {
            //plotDistr->yAxis->setScaleLogBase(100);
            plotDistr->yAxis->setScaleType(QCPAxis::stLogarithmic);
            QSharedPointer<QCPAxisTickerLog> ticker(new QCPAxisTickerLog);
            ticker->setSubTickCount(10);
            ticker->setLogBase(100);
            plotDistr->yAxis->setTicker(ticker);
            plotDistr->yAxis->setRange(1.e-1*tmin, 1.e1*tmax);
            plotDistr->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
            plotDistr->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
        }
        else {
            plotDistr->yAxis->setScaleType(QCPAxis::stLinear);
            QSharedPointer<QCPAxisTicker> ticker(new QCPAxisTicker);
            plotDistr->yAxis->setTicker(ticker);
            plotDistr->yAxis->setRange(0., 1.1*tmax);
            plotDistr->yAxis->setNumberFormat("g");
            plotDistr->yAxis->setNumberPrecision(5);
        }

        //plotDistr->xAxis->setRange(tl, tr);
        plot->setCurrentIndex(0);
        plotDistr->replot();
    }

    if (index==3 && id<spectra->fParticles.size()) {
        if (spectra->fDistributionType == 0) plot2D->xAxis->setRange(-3., 3.);
        else plot2D->xAxis->setRange(-3. - spectra->fEtaMax, 3. + spectra->fEtaMax);
        plot2D->yAxis->setRange( 0., 2.);

        QCPColorMapData *data = new QCPColorMapData(spectra->fParticles[id].d2ndptdy.szx,
                                                    spectra->fParticles[id].d2ndptdy.szy,
                                                    plot2D->xAxis->range(),
                                                    plot2D->yAxis->range());

        fXv = spectra->fParticles[id].GetXVector(3);
        fYv = spectra->fParticles[id].GetYVector(3);
        fZv = spectra->fParticles[id].GetZVector(3);
        fZvErr = spectra->fParticles[id].GetZErrorVector(3);
        for(int i=0;i<fXv.size();++i) {
           data->setData(fXv[i], fYv[i], fZv[i]);
        }

        colormap->setGradient(QCPColorGradient::gpPolar);
        colormap->setData(data, true);
        colormap->rescaleDataRange(true);

        colorScale->setGradient(colormap->gradient());
        colorScale->setDataRange(colormap->dataRange());

        delete data;

        plot2D->rescaleAxes(true);

        plot->setCurrentIndex(1);
        plot2D->replot();
    }
}

void EventGeneratorTab::replot(const QVector<double> &x1, const QVector<double> &y1, const QVector<double> &y1err,
                               const QVector<double> &x2, const QVector<double> &y2, int index, double rightlimit)
{
        double tmin = 1.e5, tmax = 0.;
        double tl = 0., tr = rightlimit;
        //double tl = 0., tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;

        plotDistr->graph(0)->data()->clear();
        plotDistr->graph(0)->setName(paramnames[index]);
        plotDistr->graph(1)->data()->clear();
        plotDistr->graph(1)->setName(paramnames[index]);

        for(int i=0;i<x1.size();++i) {
          if (y1[i] != 0.)
            tmin = std::min(tmin, y1[i]);
          tmax = std::max(tmax, y1[i]);
        }
        plotDistr->graph(0)->setData(x1, y1);
        errorBars->data()->clear();
        errorBars->setData(y1err);
        errorBars->setDataPlottable(plotDistr->graph(0));

        plotDistr->xAxis->setLabel(paramnamesx[index]);
        plotDistr->yAxis->setLabel(paramnames[index]);
        for(int i=0;i<x2.size();++i) {
          if (y2[i] != 0.)
            tmin = std::min(tmin, y2[i]);
           tmax = std::max(tmax, y2[i]);
        }
        plotDistr->graph(1)->setData(x2,y2);

        if (comboDistr->currentIndex()==2) {
            plotDistr->yAxis->setScaleType(QCPAxis::stLogarithmic);
            QSharedPointer<QCPAxisTickerLog> ticker(new QCPAxisTickerLog);
            ticker->setSubTickCount(10);
            ticker->setLogBase(100);
            plotDistr->yAxis->setTicker(ticker);
            //plotDistr->yAxis->setScaleLogBase(100);
            plotDistr->yAxis->setRange(1.e-1*tmin, 1.e1*tmax);
            plotDistr->yAxis->setNumberPrecision(0);
            plotDistr->yAxis->setNumberFormat("ebc"); // e = exponential, b = beautiful decimal powers
        }
        else {
            plotDistr->yAxis->setScaleType(QCPAxis::stLinear);
            QSharedPointer<QCPAxisTicker> ticker(new QCPAxisTicker);
            plotDistr->yAxis->setTicker(ticker);
            plotDistr->yAxis->setRange(0., 1.1*tmax);
            plotDistr->yAxis->setNumberFormat("g");
            plotDistr->yAxis->setNumberPrecision(5); 
            //plotDistr->yAxis->setSubTickCount(5);
        }

        plotDistr->xAxis->setRange(tl, tr);
        if (x2.size()>1) plotDistr->xAxis->setRange(x2[0] - 0.5*(x2[1] - x2[0]), x2[x2.size() - 1] + 0.5*(x2[1] - x2[0]));
        plot->setCurrentIndex(0);
        plotDistr->replot();
}

void EventGeneratorTab::replot2D(const QVector<double> &xv, const QVector<double> &yv, const QVector<double> &zv, int index, double rightlimit)
{
    if (spectra->fDistributionType == 0) plot2D->xAxis->setRange(-3., 3.);
    else plot2D->xAxis->setRange(-3. - spectra->fEtaMax, 3. + spectra->fEtaMax);
    plot2D->yAxis->setRange( 0., rightlimit);
    
    int id = getCurrentRow();
    if (id<0) id = 0;

    plot2D->yAxis->setRange(0., TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.);


    QCPColorMapData *data = new QCPColorMapData(spectra->fParticles[id].d2ndptdy.szx,
                                                spectra->fParticles[id].d2ndptdy.szy,
                                                plot2D->xAxis->range(),
                                                plot2D->yAxis->range());

    for(int i=0;i<xv.size();++i) {
       data->setData(xv[i], yv[i], zv[i]);
    }

    colormap->setGradient(QCPColorGradient::gpPolar);
    colormap->setData(data, true);
    colormap->rescaleDataRange(true);

    colorScale->setGradient(colormap->gradient());
    colorScale->setDataRange(colormap->dataRange());

    delete data;

    plot2D->rescaleAxes(true);

    plot->setCurrentIndex(1);
    plot2D->replot();
}

void EventGeneratorTab::calculate() {
    if (!fRunning) {
      if (radioSSH->isChecked()) {
        double betaT = spinBetat->value();
        double npow = spinn->value();
        double betaS = (2. + npow) / 2. * betaT;
        if (betaS >= 1.0) {
          QMessageBox msgBox;
          msgBox.setText(QString("Too high transverse flow! The flow velocity at the surface exceeds the speed of light, β<sub>s</sub> = %1").arg(betaS));
          msgBox.exec();
          return;
        }
      }

      generateEvents(getConfig());
    }
    else {
        fStop = 1;
    }
}


void EventGeneratorTab::finalize() {
    calcTimer->stop();
    buttonCalculate->setText(tr("Generate"));
    fRunning = false;
    progBar->setVisible(false);

    for(int i=0;i<spectra->fParticles.size();++i) spectra->fParticles[i].CalculateAverages();
    for(int i=0;i<spectra->fNetParticles.size();++i) spectra->fNetParticles[i].CalculateCentralMoments();
    for(int i=0;i<spectra->fNetCharges.size();++i) spectra->fNetCharges[i].CalculateCentralMoments();
    for(int i=0;i<spectra->fTotalCharges.size();++i) spectra->fTotalCharges[i].CalculateCentralMoments();
    for(int i=0;i<spectra->fPositiveCharges.size();++i) spectra->fPositiveCharges[i].CalculateCentralMoments();
    for(int i=0;i<spectra->fNegativeCharges.size();++i) spectra->fNegativeCharges[i].CalculateCentralMoments();

    dbgstrm << "Generated " << fCurrentSize << " events" << endl;
    dbgstrm << "Effective event number = " << nE << endl;
    dbgstrm << "CE acceptance rate: " << EventGeneratorBase::fCEAccepted / (double)(EventGeneratorBase::fCETotal) << endl;
    dbgstrm << "Calculation time = " << timer.elapsed() << " ms" << endl;
    dbgstrm << "Per event = " << timer.elapsed()/(double)(fCurrentSize) << " ms" << endl;
    dbgstrm << "----------------------------------------------------------" << endl;
    teDebug->append(dbgstr);
    dbgstr.clear();
    teDebug->verticalScrollBar()->setValue(teDebug->verticalScrollBar()->maximum());
    replot();
    myModel->updateAll();

    delete generator;
    fRunning = false;
    tableSpectra->resizeColumnsToContents();
}

void EventGeneratorTab::updateProgress() {
    progBar->setRange(0, fTotalSize);
    progBar->setValue(fCurrentSize);
    progBar->setFormat(QString::number(fCurrentSize) + " events (" + QString::number((fCurrentSize/(double)(fTotalSize)*100.),'f',2)+"%" ")");

    mutex.lock();
    for(int i=0;i<spectra->fParticles.size();++i) spectra->fParticles[i].CalculateAverages();
    for(int i=0;i<spectra->fNetParticles.size();++i) spectra->fNetParticles[i].CalculateCentralMoments();
    for(int i=0;i<spectra->fNetCharges.size();++i) spectra->fNetCharges[i].CalculateCentralMoments();
    for(int i=0;i<spectra->fTotalCharges.size();++i) spectra->fTotalCharges[i].CalculateCentralMoments();
    for(int i=0;i<spectra->fPositiveCharges.size();++i) spectra->fPositiveCharges[i].CalculateCentralMoments();
    for(int i=0;i<spectra->fNegativeCharges.size();++i) spectra->fNegativeCharges[i].CalculateCentralMoments();
    //qDebug() << "mlock";
    myModel->updateAll();
    tableSpectra->resizeColumnsToContents();
    int index = comboDistr->currentIndex();
    int id = getCurrentRow();
    if (id<0) id = 0;
    QVector<double> x1, y1, y1err, z1;
    QVector<double> x2, y2;
    if (index>=0 && index<paramnames.size() && id<spectra->fParticles.size()) {
        int ttype = comboDistr->currentIndex();
        if (ttype>=0 && ttype<=4) {
            x1 = QVector<double>::fromStdVector (spectra->fParticles[id].GetXVector(ttype) );
            y1 = QVector<double>::fromStdVector (spectra->fParticles[id].GetYVector(ttype) );
            if (ttype<=2 || ttype == 4) y1err = QVector<double>::fromStdVector (spectra->fParticles[id].GetYErrorVector(ttype) );
            else z1 = QVector<double>::fromStdVector (spectra->fParticles[id].GetZVector(ttype) );
            plotDistr->graph(0)->setData(x1, y1);
            //plotDistr->graph(0)->setDataValueError(x1,y1,y1err);
        }

        x2 = QVector<double>::fromStdVector (spectra->fParticles[id].GetModelX() );
        y2 = QVector<double>::fromStdVector (spectra->fParticles[id].GetModelY() );
    }

    //qDebug() << "munlock";


    if (((index>=0 && index<3) || index == 4) && id<spectra->fParticles.size()) {
      double tr = 2.;
      if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1)
        tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;
      //printf("%lf\n", tr);
      replot(x1,y1,y1err,x2,y2,index, tr);
    }

    if (index==3 && id<spectra->fParticles.size()) {
      double tr = 2.;
      if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1)
        tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;
      replot2D(x1,y1,z1,index, tr);
    }

    mutex.unlock();
}

void EventGeneratorTab::setModel(ThermalModelBase *modelop) {
    *model = *modelop;
}

void EventGeneratorTab::quantityDoubleClick(const QModelIndex & index) {
    int row = index.row();
    if (row>=0) {
    }
}

void EventGeneratorTab::chooseOutputFile() {
  QString path = QFileDialog::getSaveFileName(this, tr("Choose file to write events"), leFilePath->text());
  if (path.length()>0) leFilePath->setText(path);
}

void EventGeneratorTab::loadAcceptance() {
    QString path = QFileDialog::getExistingDirectory(this, tr("Open folder with acceptance data"), leAcceptancePath->text());
    if (path.length()>0) leAcceptancePath->setText(path);
}

void EventGeneratorTab::changePlot() {
    int ttype = comboDistr->currentIndex();
    int id = getCurrentRow();
    if (id<0) id = 0;
    mutex.lock();
    if (id<spectra->fParticles.size()) {
        int ttype = comboDistr->currentIndex();
        if ((ttype>=0 && ttype<=2) || ttype == 4) {
            double tl = 0.;
            if (ttype==1) tl = -3.;
            double tr = 2.;
            if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1)
              tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;
            if (ttype==1) tr = 3.;
            //if (ttype==2) tr = 2.;
            spectra->fParticles[id].FillModelDistribution(tl, tr, 100, ttype);
        }
    }
    mutex.unlock();
    if (!fRunning) replot();
}

void EventGeneratorTab::modelChanged()
{
  const ThermalModelConfig& config = configWidget->currentConfig;
  if (config.Ensemble != ThermalModelConfig::EnsembleCE) {
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
    spinmuB->setEnabled(!config.CanonicalB);
    spinmuS->setEnabled(!config.CanonicalS);
    spinmuQ->setEnabled(!config.CanonicalQ);
    spinmuC->setEnabled(!config.CanonicalC);
    spinB->setEnabled(config.CanonicalB);
    spinS->setEnabled(config.CanonicalS);
    spinQ->setEnabled(config.CanonicalQ);
    spinC->setEnabled(config.CanonicalC);

    spinVolumeRSC->setEnabled(true);
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
  else
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
  
  if (radioSR->isChecked()) {
    spinBeta->setEnabled(true);
    spinBetat->setEnabled(false);
    spinEtaMax->setEnabled(false);
    spinn->setEnabled(false);
  }
  else {
    spinBeta->setEnabled(false);
    spinBetat->setEnabled(true);
    spinEtaMax->setEnabled(true);
    if (radioSSH->isChecked())
      spinn->setEnabled(true);
    else
      spinn->setEnabled(false);
  }

  if (labelBeta->text() == "β:") {
    prevBeta = spinBeta->value();
  }
  else {
    prevRmax = spinBeta->value();
  }

  if (radioSR->isChecked()) {
    labelBeta->setText("β:");
    spinBeta->setToolTip(tr("Radial flow velocity"));
    spinBeta->setMaximum(1.0);
    spinBeta->setValue(prevBeta);
  }
  else {
    labelBeta->setText("R<sub>T</sub>(fm):");
    spinBeta->setToolTip(tr("Maximum transverse radius"));
    spinBeta->setMaximum(25.);
    spinBeta->setValue(prevRmax);
    if (radioSSH->isChecked()) {
      spinBeta->setEnabled(true);
    }
  }

  if (radioCracow->isChecked()) {
    labelBetat->setText("R/τ<sub>H</sub>:");
  }
  else {
    labelBetat->setText("〈β〉<sub>T</sub>:");
  }

  if (checkFile->isChecked()) {
    leFilePath->setEnabled(true);
    buttonChooseFile->setEnabled(true);
  }
  else {
    leFilePath->setEnabled(false);
    buttonChooseFile->setEnabled(false);
  }

  if (configWidget->currentConfig.UsePCE)
    spinTkin->setValue(configWidget->currentConfig.Tkin * 1.e3);
}

void EventGeneratorTab::resetTPS() {
    spectra->Reset(model);
    myModel->setSpectra(spectra);
    myModel->updateAll();
    tableSpectra->resizeColumnsToContents();
    
    configWidget->setModel(model);
}


void EventGeneratorTab::generateEvents(const ThermalModelConfig & config)
{
    ThermalModelBase *modelnew;

    if (config.Ensemble == ThermalModelConfig::EnsembleGCE) 
      modelnew = new ThermalModelIdeal(model->TPS());
    else if (config.Ensemble == ThermalModelConfig::EnsembleCE)
      modelnew = new ThermalModelCanonical(model->TPS());
    else if (config.Ensemble == ThermalModelConfig::EnsembleSCE) 
      modelnew = new ThermalModelCanonicalStrangeness(model->TPS());
    else 
      modelnew = new ThermalModelCanonicalCharm(model->TPS());

    myModel->setSpectra(NULL);

    configWidget->setModel(modelnew);

    if (model != NULL) delete model;
    model = modelnew;

    SetThermalModelConfiguration(model, config);
    SetThermalModelParameters(model, config);

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


    ThermalModelBase *modelEVVDW;
    ThermalParticleSystem TPSt = *model->TPS();
    if (config.InteractionModel == ThermalModelConfig::InteractionQVDW)
      modelEVVDW = new ThermalModelVDW(&TPSt);
    else if (config.InteractionModel == ThermalModelConfig::InteractionEVCrossterms)
      modelEVVDW = new ThermalModelEVCrossterms(&TPSt);
    else if (config.InteractionModel == ThermalModelConfig::InteractionEVDiagonal)
      modelEVVDW = new ThermalModelEVDiagonal(&TPSt);
    else
      modelEVVDW = new ThermalModelIdeal(&TPSt);

    SetThermalModelConfiguration(modelEVVDW, config);
    SetThermalModelInteraction(modelEVVDW, config);
    SetThermalModelParameters(modelEVVDW, config);

    // Takes time but not really needed
    //model->CalculateDensitiesGCE();

    // Convert the mean transverse velocity into the one at the surface
    double betaS = (2. + spinn->value()) / 2. * spinBetat->value();

    ParticleSpectraConfig config_spectra;
    config_spectra.fT = spinTkin->value() * 1.e-3;
    config_spectra.fEtaMax = spinEtaMax->value();
    config_spectra.fNPow = spinn->value();
    if (radioSR->isChecked()) {
      config_spectra.fBeta = spinBeta->value();
      config_spectra.fDistrType = 0;
    }
    else if (radioSSH->isChecked()) {
      config_spectra.fBeta = betaS;
      config_spectra.fDistrType = 1;
    }
    else {
      config_spectra.fBeta = betaS;
      config_spectra.fDistrType = 2;
    }

    config_spectra.fStableOnly = partsConfig.type;
    config_spectra.fPdgCodes = partsConfig.pdgCodes;

    config_spectra.fBins  = binConfig.bins1D;
    config_spectra.fBinsX = binConfig.bins2D_x;
    config_spectra.fBinsY = binConfig.bins2D_y;

    spectra->Reset(model, config_spectra);

    //if (radioSR->isChecked()) 
    //  spectra->Reset(model, spinTkin->value() * 1.e-3, spinBeta->value());
    //else if (radioSSH->isChecked()) 
    //  spectra->Reset(model, spinTkin->value() * 1.e-3, betaS, 1, spinEtaMax->value(), spinn->value());
    //else 
    //  spectra->Reset(model, spinTkin->value() * 1.e-3, betaS, 2, spinEtaMax->value());

    int id = getCurrentRow();
    if (id<0) id = 0;
    if (id<spectra->fParticles.size()) {
      int ttype = comboDistr->currentIndex();
      if ((ttype >= 0 && ttype <= 2) || ttype == 4) {
        double tl = 0.;
        if (ttype == 1) tl = -3.;
        double tr = 2.;
        if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1)
          tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;
        if (ttype == 1) tr = 3.;
        if (ttype == 1 && spectra->fDistributionType == 1) {
          tl = -3. - spectra->fEtaMax;
          tr = 3. + spectra->fEtaMax;
        }
        //if (ttype==2) tr = 2.;
        spectra->fParticles[id].FillModelDistribution(tl, tr, 100, ttype);
      }
    }

    replot();

    myModel->setSpectra(spectra);

    EventGeneratorConfiguration configMC;
    configMC.fModelType = EventGeneratorConfiguration::PointParticle;

    if (config.InteractionModel == ThermalModelConfig::InteractionEVDiagonal)
      configMC.fModelType = EventGeneratorConfiguration::DiagonalEV;
    if (config.InteractionModel == ThermalModelConfig::InteractionEVCrossterms)
      configMC.fModelType = EventGeneratorConfiguration::CrosstermsEV;
    if (config.InteractionModel == ThermalModelConfig::InteractionQVDW)
      configMC.fModelType = EventGeneratorConfiguration::QvdW;

    configMC.fEnsemble = EventGeneratorConfiguration::GCE;
    if (model->Ensemble() == ThermalModelBase::CE)
      configMC.fEnsemble = EventGeneratorConfiguration::CE;
    if (model->Ensemble() == ThermalModelBase::SCE)
      configMC.fEnsemble = EventGeneratorConfiguration::SCE;
    if (model->Ensemble() == ThermalModelBase::CCE)
      configMC.fEnsemble = EventGeneratorConfiguration::CCE;

    configMC.CFOParameters = model->Parameters();
    configMC.B = model->Parameters().B;
    configMC.Q = model->Parameters().Q;
    configMC.S = model->Parameters().S;
    configMC.C = model->Parameters().C;
    configMC.CanonicalB = config.CanonicalB;
    configMC.CanonicalQ = config.CanonicalQ;
    configMC.CanonicalS = config.CanonicalS;
    configMC.CanonicalC = config.CanonicalC;

    configMC.bij = CuteHRGHelper::bijMatrix(modelEVVDW);

    configMC.aij = CuteHRGHelper::aijMatrix(modelEVVDW);

    if (config.UsePCE) {
      modelEVVDW->SolveChemicalPotentials(
        config.B, config.Q, config.S, config.C,
        config.muB, config.muQ, config.muS, config.muC,
        config.CanonicalB, config.CanonicalQ, config.CanonicalS, config.CanonicalC
        );

      ThermalModelPCE modelpce(modelEVVDW);
      modelpce.SetStabilityFlags(ThermalModelPCE::ComputePCEStabilityFlags(modelEVVDW->TPS(), config.PCESahaForNuclei, config.PCEFreezeLongLived, config.PCEWidthCut));
      modelpce.SetChemicalFreezeout(modelEVVDW->Parameters(), modelEVVDW->ChemicalPotentials());
      modelpce.CalculatePCE(config.Tkin);

      configMC.fUsePCE = true;
      configMC.fPCEChems = modelpce.ChemicalPotentials();
      configMC.CFOParameters = modelEVVDW->Parameters();
    }

    configMC.fUseEVRejectionMultiplicity = config.fUseEVRejectionMultiplicity;
    configMC.fUseEVRejectionCoordinates = config.fUseEVRejectionCoordinates;
    configMC.fUseEVUseSPRApproximation = config.fUseEVUseSPRApproximation;

    if (radioSR->isChecked())
      generator = new SphericalBlastWaveEventGenerator(model->TPS(), configMC, spinTkin->value() * 1.e-3, spinBeta->value());
    else if (radioSSH->isChecked())
      generator = new CylindricalBlastWaveEventGenerator(model->TPS(), configMC, spinTkin->value() * 1.e-3, betaS, spinEtaMax->value(), spinn->value(), spinBeta->value());
    else
      generator = new CracowFreezeoutEventGenerator(model->TPS(), configMC, spinTkin->value() * 1.e-3, betaS, spinEtaMax->value());

    delete modelEVVDW;
    
    //std::vector<Acceptance::AcceptanceFunction>& tacc = generator->GetAcceptance();
    //for (int i = 0; i<spectra->fParticles.size(); ++i) {
    //  int tind = model->TPS()->PdgToId(spectra->fParticles[i].GetPDGID());
    //  if (tacc.size()>tind && tacc[tind].init) spectra->fParticles[i].SetAcceptance(true);
    //  else spectra->fParticles[i].SetAcceptance(false);
    //}

    timer.start();

    fTotalSize = spinEvents->value();
    fCurrentSize = 0;
    fStop = 0;


    std::string filepath = "";
    if (checkFile->isChecked())
      filepath = leFilePath->text().toStdString();

    EventGeneratorWorker *wrk = new EventGeneratorWorker(generator, spectra, &mutex, fTotalSize, &fCurrentSize, &fStop, &nE, checkDecays->isChecked(), filepath);

    connect(wrk, SIGNAL(calculated()), this, SLOT(finalize()));
    connect(wrk, SIGNAL(finished()), wrk, SLOT(deleteLater()));

    progBar->setFormat("");

    wrk->start();

    buttonCalculate->setText(tr("Stop"));
    fRunning = true;
    progBar->setVisible(true);

    calcTimer->start(150);
}

void EventGeneratorTab::changeVolumeRSC(double VRSC)
{
  spinVolumeRSC->setValue(VRSC);
  double R = spinVolumeR->value();
  double V = 4. / 3. * xMath::Pi() * R * R * R;
  //spinVolumeR->setToolTip(tr("System radius, current value corresponds to the volume V = ")
  //  + QString("%1").arg(V)
  //  + tr(" fm^3"));
  labelVolumeVal->setText(QString("%1").arg(V));

  double Rc = spinVolumeRSC->value();
  double Vc = 4. / 3. * xMath::Pi() * Rc * Rc * Rc;
  //spinVolumeRSC->setToolTip(tr("System correlation radius, current value corresponds to the volume Vc = ")
  //  + QString("%1").arg(Vc)
  //  + tr(" fm^3"));
  labelVolumeSCVal->setText(QString("%1").arg(Vc));
}

void EventGeneratorTab::changeTkin(double Tch)
{
  if (!configWidget->currentConfig.UsePCE) 
    spinTkin->setValue(Tch);
}

void EventGeneratorTab::contextMenuRequestPlot1D(QPoint pos)
{
  QMenu* menu = new QMenu(this);
  menu->setAttribute(Qt::WA_DeleteOnClose);

  menu->addAction("Save as pdf", this, SLOT(saveAsPdf1D()));
  menu->addAction("Save as png", this, SLOT(saveAsPng1D()));
  menu->addAction("Write computed values to file", this, SLOT(saveAsAscii1D()));

  menu->popup(plotDistr->mapToGlobal(pos));
}

void EventGeneratorTab::saveAs1D(int type)
{
  QString tname = plotDistr->yAxis->label();
  for (int i = 0; i < tname.size(); ++i) {
    if (tname[i] == '/')
      tname[i] = '.';
  }

  {
    int id = getCurrentRow();
    if (id >= 0 && id < spectra->fParticles.size()) {
      if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1) {
        tname = QString::fromStdString(TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Name() + ".") + tname;
      }
    }
  }

  QVector<QString> exts;
  exts.push_back("pdf");
  exts.push_back("png");
  exts.push_back("dat");

  QString listpathprefix = cpath + "/" + tname + "." + exts[type];
  QString path = QFileDialog::getSaveFileName(this, "Save plot as " + exts[type], listpathprefix, "*." + exts[type]);

  if (path.length() > 0)
  {
    if (type == 0)
      plotDistr->savePdf(path, plotDistr->width(), plotDistr->height());
    else if (type == 1)
      plotDistr->savePng(path, plotDistr->width(), plotDistr->height());
    else {
      QFile fout(path);

      if (fout.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&fout);
        out.setFieldWidth(15);
        out.setFieldAlignment(QTextStream::AlignLeft);
        out << plotDistr->xAxis->label();
        out << plotDistr->yAxis->label();
        out << "error";
        out << qSetFieldWidth(0) << endl << qSetFieldWidth(15);
        for (int i = 0; i < plotDistr->graph(0)->data()->size(); ++i) {
          out.setFieldWidth(15);
          out.setFieldAlignment(QTextStream::AlignLeft);
          double x = plotDistr->graph(0)->data()->at(i)->key;
          double y = plotDistr->graph(0)->data()->at(i)->value;
          double yerr = errorBars->data()->at(i).errorPlus;
          out << x << y << yerr;
          out << qSetFieldWidth(0) << endl << qSetFieldWidth(15);
        }
      }
    }

    QFileInfo saved(path);
    cpath = saved.absolutePath();
  }
}

void EventGeneratorTab::contextMenuRequestPlot2D(QPoint pos)
{
  QMenu* menu = new QMenu(this);
  menu->setAttribute(Qt::WA_DeleteOnClose);

  menu->addAction("Save as pdf", this, SLOT(saveAsPdf2D()));
  menu->addAction("Save as png", this, SLOT(saveAsPng2D()));
  menu->addAction("Write computed values to file", this, SLOT(saveAsAscii2D()));

  menu->popup(plot2D->mapToGlobal(pos));
}

void EventGeneratorTab::saveAs2D(int type)
{
  QString plotName = "d2N/dpTdy";
  QString tname = plotName;
  for (int i = 0; i < tname.size(); ++i) {
    if (tname[i] == '/')
      tname[i] = '.';
  }

  {
    int id = getCurrentRow();
    if (id >= 0 && id < spectra->fParticles.size()) {
      if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1) {
        tname = QString::fromStdString(TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Name() + ".") + tname;
      }
    }
  }

  QVector<QString> exts;
  exts.push_back("pdf");
  exts.push_back("png");
  exts.push_back("dat");

  QString listpathprefix = cpath + "/" + tname + "." + exts[type];
  QString path = QFileDialog::getSaveFileName(this, "Save plot as " + exts[type], listpathprefix, "*." + exts[type]);

  if (path.length() > 0)
  {
    if (type == 0)
      plot2D->savePdf(path, plotDistr->width(), plotDistr->height());
    else if (type == 1)
      plot2D->savePng(path, plotDistr->width(), plotDistr->height());
    else {
      QFile fout(path);

      if (fout.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&fout);
        out.setFieldWidth(15);
        out.setFieldAlignment(QTextStream::AlignLeft);
        out << plot2D->xAxis->label();
        out << plot2D->yAxis->label();
        out << plotName;
        out << "error";
        out << qSetFieldWidth(0) << endl << qSetFieldWidth(15);
        for (int i = 0; i < fXv.size(); ++i) {
          out.setFieldWidth(15);
          out.setFieldAlignment(QTextStream::AlignLeft);
          out << fXv[i] << fYv[i] << fZv[i] << fZvErr[i];
          out << qSetFieldWidth(0) << endl << qSetFieldWidth(15);
        }
      }
    }

    QFileInfo saved(path);
    cpath = saved.absolutePath();
  }
}

void EventGeneratorTab::changeParticles()
{
  ParticlesAnalyzeDialog dialog(&partsConfig, this);
  dialog.setWindowFlags(Qt::Window);
  dialog.exec();
}


void EventGeneratorTab::changeBinning()
{
  BinningDialog dialog(&binConfig, this);
  dialog.setWindowFlags(Qt::Window);
  dialog.exec();
}
