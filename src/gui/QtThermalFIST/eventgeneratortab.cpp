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
#include "HRGEventGenerator/SREventGenerator.h"
#include "HRGEventGenerator/SSHEventGenerator.h"

#include "DebugText.h"
#include "spectramodel.h"
#include "particlespectra.h"
#include "qcustomplot.h"

using namespace thermalfist;

void EventGeneratorWorker::run()
{
     if (mutex!=NULL) {
         for(int i=0;i<events && !(*stop);++i) {
             SimpleEvent ev = generator->GetEvent(performDecays);
             mutex->lock();
             wsum += ev.weight;
             w2sum += ev.weight*ev.weight;
             spectra->ProcessEvent(ev);
             (*eventsProcessed)++;

             if (fout.is_open())
               ev.writeToFile(fout, (*eventsProcessed));

             mutex->unlock();
         }
     }
     if (fout.is_open())
      fout.close();
     *nE = wsum * wsum / w2sum;
     emit calculated();
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

    //fitcopy = NULL;
    dbgstrm.setString(&dbgstr);

    TPS = modelop->TPS();
    model = new ThermalModelIdeal(modelop->TPS());
    *model = *modelop;
    spectra = new ParticlesSpectra(model);


    QHBoxLayout *DataEditLay = new QHBoxLayout();

    QVBoxLayout *dataLayv = new QVBoxLayout();
    QVBoxLayout *dataLayv2 = new QVBoxLayout();

    QLabel *labelParticleList = new QLabel(tr("Particle list:"));
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
    selLay->addWidget(labelSel);
    selLay->addWidget(comboDistr);

    plotDistr = new QCustomPlot();
    plotDistr->xAxis->setLabel("p (GeV)");
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

    errorBars = new QCPErrorBars(plotDistr->xAxis, plotDistr->yAxis);
    errorBars->removeFromLegend();
    errorBars->setErrorType(QCPErrorBars::etValueError);
    errorBars->setPen(QPen(Qt::blue));
    errorBars->setSymbolGap(5);
    //errorBars->set
    //errorBars->setDataPlottable(plotDistr->graph(0));

    plot2D = new QCustomPlot();
    plot2D->xAxis->setLabel(tr("y"));
    plot2D->yAxis->setLabel(tr("pT (GeV)"));

    colormap = new QCPColorMap(plot2D->xAxis, plot2D->yAxis);

    //plot2D->addPlottable(colormap);

    colorScale = new QCPColorScale(plot2D);
    plot2D->plotLayout()->addElement(0, 1, colorScale);

    plot = new QStackedWidget();
    plot->addWidget(plotDistr);

    plot->addWidget(plot2D);
    if (comboDistr->currentIndex()==3) plot->setCurrentIndex(1);
    else plot->setCurrentIndex(0);

    dataLayv->addWidget(labelParticleList);
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


    QGroupBox *grParameters = new QGroupBox(tr("Thermal model parameters:"));

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
    QLabel *labelVolume = new QLabel(tr("R (fm):"));
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
    spinB->setValue(2);
    labelS = new QLabel(tr("S:"));
    spinS = new QSpinBox();
    spinS->setMinimum(-1000);
    spinS->setMaximum(1000);
    spinS->setValue(0);
    labelQ = new QLabel(tr("Q:"));
    spinQ = new QSpinBox();
    spinQ->setMinimum(-1000);
    spinQ->setMaximum(1000);
    spinQ->setValue(2);
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
    layParameters->addWidget(labelVolume, 2, 0, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinVolumeR, 2, 1);
    layParameters->addWidget(labelVolumeRSC, 2, 2, 1, 1, Qt::AlignRight);
    layParameters->addWidget(spinVolumeRSC, 2, 3);

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
    radioSSH->setToolTip(tr("Schnedermann-Sollfrank-Heinz prescription (note that initialization takes while!)"));
    radioSR->setChecked(true);
    connect(radioSR, SIGNAL(clicked()), this, SLOT(modelChanged()));
    connect(radioSSH, SIGNAL(clicked()), this, SLOT(modelChanged()));


    layMom1->addWidget(radioSR);
    layMom1->addWidget(radioSSH);


    QHBoxLayout *layMom2 = new QHBoxLayout();
    layMom2->setAlignment(Qt::AlignLeft);

    QLabel *labelTkin = new QLabel(tr("T<sub>kin</sub> (MeV):"));
    spinTkin = new QDoubleSpinBox();
    spinTkin->setMinimum(1.);
    spinTkin->setMaximum(10000.);
    spinTkin->setValue(model->Parameters().T * 1e3);
    spinTkin->setToolTip(tr("Kinetic freeze-out temperature which fixes the momentum spectrum"));
    QLabel *labelBeta = new QLabel(tr("β:"));
    spinBeta = new QDoubleSpinBox();
    spinBeta->setMinimum(0.);
    spinBeta->setMaximum(1.);
    spinBeta->setDecimals(3);
    spinBeta->setValue(0.5);
    spinBeta->setToolTip(tr("Radial flow velocity"));
    QLabel *labelBetat = new QLabel(tr("β<sub>T</sub>:"));
    spinBetat = new QDoubleSpinBox();
    spinBetat->setMinimum(0.);
    spinBetat->setMaximum(1.);
    spinBetat->setDecimals(3);
    spinBetat->setValue(0.5);
    spinBetat->setToolTip(tr("Transverse radial flow parameter"));
    QLabel *labelEtaMax = new QLabel(tr("η<sub>max</sub>:"));
    spinEtaMax = new QDoubleSpinBox();
    spinEtaMax->setMinimum(0.);
    spinEtaMax->setMaximum(100.);
    spinEtaMax->setDecimals(3);
    spinEtaMax->setValue(1.);
    spinEtaMax->setToolTip(tr("Longitudinal space-time rapidity cut-off"));
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
}



EventGeneratorTab::~EventGeneratorTab()
{
    delete spectra;
    delete model;
}

ThermalModelConfig EventGeneratorTab::getConfig()
{
  ThermalModelConfig ret = configWidget->updatedConfig();

  ret.QuantumStatistics = 0;
  ret.QuantumStatisticsType = 0;
  ret.QuantumStatisticsInclude = 0;

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
        if (ttype>=0 && ttype<=2) {
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
    if (index>=0 && index<3 && id<spectra->fParticles.size()) {

        double tmin = 1.e5, tmax = 0.;
        double tl = 0., tr = 2.;
        if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1)
          tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;

        plotDistr->graph(0)->data()->clear();
        plotDistr->graph(0)->setName(paramnames[index]);
        plotDistr->graph(1)->data()->clear();
        plotDistr->graph(1)->setName(paramnames[index]);

        int ttype = comboDistr->currentIndex();
        if (ttype>=0 && ttype<=2) {
            QVector<double> x1,y1,y1err;
            x1    = QVector<double>::fromStdVector (spectra->fParticles[id].GetXVector(ttype) );
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
        plot2D->xAxis->setRange(-3., 3.);
        if (spectra->fDistributionType==1) plot2D->xAxis->setRange(-3. - spectra->fEtaMax, 3. + spectra->fEtaMax);
        plot2D->yAxis->setRange( 0., 2.);

        QCPColorMapData *data = new QCPColorMapData(spectra->fParticles[id].d2ndptdy.szx,
                                                    spectra->fParticles[id].d2ndptdy.szy,
                                                    plot2D->xAxis->range(),
                                                    plot2D->yAxis->range());

        std::vector<double> xv = spectra->fParticles[id].GetXVector(3);
        std::vector<double> yv = spectra->fParticles[id].GetYVector(3);
        std::vector<double> zv = spectra->fParticles[id].GetZVector(3);
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
    plot2D->xAxis->setRange(-3., 3.);
    //if (spectra->fDistributionType==1) plot2D->xAxis->setRange(-3. - spectra->fEtaMax, 3. + spectra->fEtaMax);
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
        if (ttype>=0 && ttype<=3) {
            x1 = QVector<double>::fromStdVector (spectra->fParticles[id].GetXVector(ttype) );
            y1 = QVector<double>::fromStdVector (spectra->fParticles[id].GetYVector(ttype) );
            if (ttype<=2) y1err = QVector<double>::fromStdVector (spectra->fParticles[id].GetYErrorVector(ttype) );
            else z1 = QVector<double>::fromStdVector (spectra->fParticles[id].GetZVector(ttype) );
            plotDistr->graph(0)->setData(x1, y1);
            //plotDistr->graph(0)->setDataValueError(x1,y1,y1err);
        }

        x2 = QVector<double>::fromStdVector (spectra->fParticles[id].GetModelX() );
        y2 = QVector<double>::fromStdVector (spectra->fParticles[id].GetModelY() );
    }

    //qDebug() << "munlock";


    if (index>=0 && index<3 && id<spectra->fParticles.size()) {
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
        if (ttype>=0 && ttype<=2) {
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
    spinmuB->setEnabled(false);
    spinmuS->setEnabled(false);
    spinmuQ->setEnabled(false);
    spinmuC->setEnabled(false);
    spinB->setEnabled(true);
    spinS->setEnabled(true);
    spinQ->setEnabled(true);
    spinC->setEnabled(true);

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
    spinn->setEnabled(true);
  }

  if (checkFile->isChecked()) {
    leFilePath->setEnabled(true);
    buttonChooseFile->setEnabled(true);
  }
  else {
    leFilePath->setEnabled(false);
    buttonChooseFile->setEnabled(false);
  }
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

    if (config.Ensemble == ThermalModelConfig::EnsembleGCE) modelnew = new ThermalModelIdeal(model->TPS());
    else if (config.Ensemble == ThermalModelConfig::EnsembleCE) {
      modelnew = new ThermalModelCanonical(model->TPS());
      ThermalModelParameters params;
      params.T = 0.100;
      params.gammaC = 1.;
      params.gammaS = 1.;
      params.gammaq = 1.;
      params.V = 20.;
      params.B = config.B;
      params.Q = config.Q;
      params.S = config.S;
      params.C = config.C;
      modelnew->SetParameters(params);
    }
    else if (config.Ensemble == ThermalModelConfig::EnsembleSCE) modelnew = new ThermalModelCanonicalStrangeness(model->TPS());
    else modelnew = new ThermalModelCanonicalCharm(model->TPS());

    myModel->setSpectra(NULL);

    configWidget->setModel(modelnew);

    if (model != NULL) delete model;
    model = modelnew;

    SetThermalModelConfiguration(model, config);


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

    SetThermalModelInteraction(modelEVVDW, config);


    model->CalculateDensitiesGCE();



    if (radioSR->isChecked()) spectra->Reset(model, spinTkin->value() * 1.e-3, spinBeta->value());
    else if (radioSSH->isChecked()) spectra->Reset(model, spinTkin->value() * 1.e-3, spinBetat->value(), 1, spinEtaMax->value(), spinn->value());

    int id = getCurrentRow();
    if (id<0) id = 0;
    if (id<spectra->fParticles.size()) {
      int ttype = comboDistr->currentIndex();
      if (ttype >= 0 && ttype <= 2) {
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

    configMC.bij = CuteHRGHelper::bijMatrix(modelEVVDW);

    configMC.aij = CuteHRGHelper::aijMatrix(modelEVVDW);

    if (radioSR->isChecked())
      generator = new SphericalBlastWaveEventGenerator(model->TPS(), configMC, spinTkin->value() * 1.e-3, spinBeta->value());
    else 
      generator = new CylindricalBlastWaveEventGenerator(model->TPS(), configMC, spinTkin->value() * 1.e-3, spinBetat->value(), spinEtaMax->value(), spinn->value());


    //if (radioSR->isChecked()) generator = new SphericalBlastWaveEventGenerator(model, spinTkin->value() * 1.e-3, spinBeta->value(), false, EV, modelEVVDW);
    //else generator = new CylindricalBlastWaveEventGenerator(model, spinTkin->value() * 1.e-3, spinBetat->value(), spinEtaMax->value(), spinn->value(), false, EV, modelEVVDW);

    delete modelEVVDW;
    
    std::vector<Acceptance::AcceptanceFunction>& tacc = generator->GetAcceptance();
    for (int i = 0; i<spectra->fParticles.size(); ++i) {
      int tind = model->TPS()->PdgToId(spectra->fParticles[i].GetPDGID());
      if (tacc.size()>tind && tacc[tind].init) spectra->fParticles[i].SetAcceptance(true);
      else spectra->fParticles[i].SetAcceptance(false);
    }



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
}

void EventGeneratorTab::changeTkin(double Tch)
{
  spinTkin->setValue(Tch);
}