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
#include "HRGVDW/ThermalModelVDWFull.h"
#include "HRGBase/ThermalModelCanonical.h"
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonicalCharm.h"
#include "HRGEventGenerator/SREventGenerator.h"
#include "HRGEventGenerator/SSHEventGenerator.h"

#include "DebugText.h"
#include "spectramodel.h"
#include "particlespectra.h"
#include "QCustomPlot/qcustomplot.h"

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
    plotDistr->graph(0)->setErrorType(QCPGraph::etValue);
    plotDistr->addGraph();
    plotDistr->graph(1)->setName(comboDistr->currentText());
    plotDistr->graph(1)->setPen(QPen(Qt::blue, 2, Qt::DashLine));

    plot2D = new QCustomPlot();
    plot2D->xAxis->setLabel(tr("y"));
    plot2D->yAxis->setLabel(tr("pt (GeV)"));

    colormap = new QCPColorMap(plot2D->xAxis, plot2D->yAxis);

    plot2D->addPlottable(colormap);

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


		QHBoxLayout *layModelEnsemble = new QHBoxLayout();
		layModelEnsemble->setAlignment(Qt::AlignLeft);

    QGroupBox *grModels = new QGroupBox(tr("HRG Model:"));

    QHBoxLayout *layModels = new QHBoxLayout();
    layModels->setAlignment(Qt::AlignLeft);
    radIdeal      = new QRadioButton(tr("Ideal"));
		radEVD        = new QRadioButton(tr("Diagonal EV"));
		radQVDW       = new QRadioButton(tr("QvdW"));
    connect(radIdeal,				SIGNAL(clicked()), this, SLOT(modelChanged()));
    connect(radEVD,					SIGNAL(clicked()), this, SLOT(modelChanged()));
		connect(radQVDW,					SIGNAL(clicked()), this, SLOT(modelChanged()));

    layModels->addWidget(radIdeal);
    layModels->addWidget(radEVD);
		layModels->addWidget(radQVDW);

    grModels->setLayout(layModels);

		radIdeal->setChecked(true);

    QGroupBox *grEnsemble = new QGroupBox(tr("Ensemble:"));

    QHBoxLayout *layEnsemble = new QHBoxLayout();
    layEnsemble->setAlignment(Qt::AlignLeft);
    radGCE   = new QRadioButton(tr("GCE"));
    radCE    = new QRadioButton(tr("CE" ));
    radSCE   = new QRadioButton(tr("SCE"));
    connect(radGCE, SIGNAL(clicked()), this, SLOT(modelChanged()));
    connect(radCE,  SIGNAL(clicked()), this, SLOT(modelChanged()));
    connect(radSCE, SIGNAL(clicked()), this, SLOT(modelChanged()));

    layEnsemble->addWidget(radGCE);
    layEnsemble->addWidget(radCE);
    layEnsemble->addWidget(radSCE);

    grEnsemble->setLayout(layEnsemble);

    radGCE->setChecked(true);

		layModelEnsemble->addWidget(grModels);
		layModelEnsemble->addWidget(grEnsemble);

    QGroupBox *grEV = new QGroupBox(tr("Excluded volume/van der Waals:"));

    QHBoxLayout *layRadius = new QHBoxLayout();
    layRadius->setAlignment(Qt::AlignLeft);
    QLabel *labelRadius = new QLabel(tr("Radius (fm):"));
    spinRadius = new QDoubleSpinBox();
    spinRadius->setMinimum(0.);
    spinRadius->setMaximum(100.);
    spinRadius->setValue(0.3);
    radioUniform  = new QRadioButton(tr("Same for all"));
    radioBaglike  = new QRadioButton(tr("Bag-like"));
    radioMesons   = new QRadioButton(tr("Point-like mesons"));
		radioCustomEV = new QRadioButton(tr("Custom..."));
		strEVPath = "";
		connect(radioCustomEV, SIGNAL(clicked()), this, SLOT(loadEVFromFile()));
    layRadius->addWidget(labelRadius);
    layRadius->addWidget(spinRadius);
    layRadius->addWidget(radioUniform);
    layRadius->addWidget(radioBaglike);
    layRadius->addWidget(radioMesons);
		layRadius->addWidget(radioCustomEV);
    radioUniform->setChecked(true);

    grEV->setLayout(layRadius);

    QHBoxLayout *layEnergy = new QHBoxLayout();
    layEnergy->setAlignment(Qt::AlignLeft);
    QLabel *labelEnergy = new QLabel(tr("e<sup>kin</sup> (GeV):"));
    spinEnergy = new QDoubleSpinBox();
    spinEnergy->setMinimum(0.);
    spinEnergy->setMaximum(1000000.);
    spinEnergy->setValue(25.);
    connect(spinEnergy, SIGNAL(valueChanged(double)), this, SLOT(updateThermalParameters()));
    layEnergy->addWidget(labelEnergy);
    layEnergy->addWidget(spinEnergy);


    QGroupBox *grParameters = new QGroupBox(tr("Thermal model parameters:"));

		QGridLayout *layParameters = new QGridLayout();
		layParameters->setAlignment(Qt::AlignLeft);
		//layQuant->setHorizontalSpacing(15);
		QLabel *labelTemperature = new QLabel(tr("T<sub>ch</sub> (MeV):"));
		spinTemperature = new QDoubleSpinBox();
		spinTemperature->setMinimum(1.);
		spinTemperature->setMaximum(10000.);
		spinTemperature->setValue(model->Parameters().T * 1e3);
		connect(spinTemperature, SIGNAL(valueChanged(double)), this, SLOT(changeTkin(double)));
		QLabel *labelmuB = new QLabel(tr("μ<sub>B</sub> (MeV):"));
		spinmuB = new QDoubleSpinBox();
		spinmuB->setMinimum(-1000.);
		spinmuB->setMaximum(1000.);
		spinmuB->setValue(model->Parameters().muB * 1e3);
		QLabel *labelgammaq = new QLabel(tr("γ<sub>q</sub>:"));
		spingammaq = new QDoubleSpinBox();
		spingammaq->setMinimum(0.);
		spingammaq->setMaximum(10.);
		spingammaq->setDecimals(4);
		spingammaq->setValue(model->Parameters().gammaq);
		QLabel *labelgammaS = new QLabel(tr("γ<sub>S</sub>:"));
		spingammaS = new QDoubleSpinBox();
		spingammaS->setMinimum(0.);
		spingammaS->setMaximum(10.);
		spingammaS->setDecimals(4);
		spingammaS->setValue(model->Parameters().gammaS);

		labelmuS = new QLabel(tr("μ<sub>S</sub> (MeV):"));
		spinmuS = new QDoubleSpinBox();
		spinmuS->setMinimum(-1000.);
		spinmuS->setMaximum(1000.);
		spinmuS->setValue(model->Parameters().muS * 1e3);
		QLabel *labelmuQ = new QLabel(tr("μ<sub>Q</sub> (MeV):"));
		spinmuQ = new QDoubleSpinBox();
		spinmuQ->setMinimum(-1000.);
		spinmuQ->setMaximum(1000.);
		spinmuQ->setValue(model->Parameters().muQ * 1e3);
		labelmuC = new QLabel(tr("μ<sub>C</sub> (MeV):"));
		spinmuC = new QDoubleSpinBox();
		spinmuC->setMinimum(-1000.);
		spinmuC->setMaximum(1000.);
		spinmuC->setValue(model->Parameters().muC * 1e3);
		QLabel *labelVolume = new QLabel(tr("R (fm):"));
		spinVolumeR = new QDoubleSpinBox();
		spinVolumeR->setMinimum(0.);
		spinVolumeR->setMaximum(25.);
		spinVolumeR->setDecimals(4);
		spinVolumeR->setValue(8.);
		connect(spinVolumeR, SIGNAL(valueChanged(double)), this, SLOT(changeVolumeRSC(double)));

		QLabel *labelB = new QLabel(tr("B:"));
		spinB = new QSpinBox();
		spinB->setMinimum(-1000);
		spinB->setMaximum(1000);
		spinB->setValue(2);
		QLabel *labelS = new QLabel(tr("S:"));
		spinS = new QSpinBox();
		spinS->setMinimum(-1000);
		spinS->setMaximum(1000);
		spinS->setValue(0);
		QLabel *labelQ = new QLabel(tr("Q:"));
		spinQ = new QSpinBox();
		spinQ->setMinimum(-1000);
		spinQ->setMaximum(1000);
		spinQ->setValue(2);

		layParameters->addWidget(labelTemperature, 0, 0, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinTemperature, 0, 1);
		layParameters->addWidget(labelmuB, 1, 0, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinmuB, 1, 1);
		layParameters->addWidget(labelgammaq, 0, 2, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spingammaq, 0, 3);
		layParameters->addWidget(labelgammaS, 0, 4, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spingammaS, 0, 5);

		layParameters->addWidget(labelmuS, 1, 4, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinmuS, 1, 5);
		layParameters->addWidget(labelmuQ, 1, 2, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinmuQ, 1, 3);
		layParameters->addWidget(labelmuC, 1, 6, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinmuC, 1, 7);
		layParameters->addWidget(labelVolume, 0, 6, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinVolumeR, 0, 7);

		layParameters->addWidget(labelB, 2, 0, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinB, 2, 1);
		layParameters->addWidget(labelS, 2, 4, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinS, 2, 5);
		layParameters->addWidget(labelQ, 2, 2, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinQ, 2, 3);

		grParameters->setLayout(layParameters);

		QHBoxLayout *layQB = new QHBoxLayout();
		layQB->setAlignment(Qt::AlignLeft);

		checkFixMuQ = new QCheckBox(tr("Constrain μQ"));
		checkFixMuQ->setChecked(false); 

		QLabel *labelQB = new QLabel(tr("Q/B ratio:"));
		spinQBRatio = new QDoubleSpinBox();
		spinQBRatio->setDecimals(3);
		spinQBRatio->setMinimum(0.);
		spinQBRatio->setValue(model->QoverB());

		checkFixMuS = new QCheckBox(tr("Constrain μS"));
		checkFixMuS->setChecked(false);  
		checkFixMuC = new QCheckBox(tr("Constrain μC"));
		checkFixMuC->setChecked(false);  

    //connect(checkFixMuQ, SIGNAL(clicked()), this, SLOT(modelChanged()));
    //connect(checkFixMuS, SIGNAL(clicked()), this, SLOT(modelChanged()));
    //connect(checkFixMuC, SIGNAL(clicked()), this, SLOT(modelChanged()));

		QLabel *labelVolumeRSC = new QLabel(tr("R<sub>SC</sub>:"));
		spinVolumeRSC = new QDoubleSpinBox();
		spinVolumeRSC->setDecimals(4);
		spinVolumeRSC->setMinimum(0.);
		spinVolumeRSC->setMaximum(25.);
		spinVolumeRSC->setValue(spinVolumeR->value());
		spinVolumeRSC->setEnabled(false);


		layQB->addWidget(labelVolumeRSC);
		layQB->addWidget(spinVolumeRSC);


    QGroupBox *grParametersMom = new QGroupBox(tr("Blast-wave momentum spectrum:"));
    QVBoxLayout *layMom = new QVBoxLayout();

    QHBoxLayout *layMom1 = new QHBoxLayout();
    layMom1->setAlignment(Qt::AlignLeft);
    radioSR = new QRadioButton(tr("Spherically symmetric"));
    radioSSH = new QRadioButton(tr("Cylindrically symmetric"));
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
    QLabel *labelBeta = new QLabel(tr("β:"));
    spinBeta = new QDoubleSpinBox();
    spinBeta->setMinimum(0.);
    spinBeta->setMaximum(1.);
    spinBeta->setDecimals(3);
    spinBeta->setValue(0.5);
    QLabel *labelBetat = new QLabel(tr("β<sub>T</sub>:"));
    spinBetat = new QDoubleSpinBox();
    spinBetat->setMinimum(0.);
    spinBetat->setMaximum(1.);
    spinBetat->setDecimals(3);
    spinBetat->setValue(0.5);
    QLabel *labelEtaMax = new QLabel(tr("η<sub>max</sub>:"));
    spinEtaMax = new QDoubleSpinBox();
    spinEtaMax->setMinimum(0.);
    spinEtaMax->setMaximum(100.);
    spinEtaMax->setDecimals(3);
    spinEtaMax->setValue(1.);
		QLabel *labeln = new QLabel(tr("n:"));
		spinn = new QDoubleSpinBox();
		spinn->setMinimum(0.);
		spinn->setMaximum(10.);
		spinn->setDecimals(3);
		spinn->setValue(1.);

		layMom2->addWidget(labelTkin);
		layMom2->addWidget(spinTkin);
		layMom2->addSpacing(20);
    layMom2->addWidget(labelBeta);
    layMom2->addWidget(spinBeta);
    layMom2->addSpacing(20);
    layMom2->addWidget(labelBetat);
    layMom2->addWidget(spinBetat);
    layMom2->addSpacing(5);
    layMom2->addWidget(labelEtaMax);
    layMom2->addWidget(spinEtaMax);
		layMom2->addSpacing(5);
		layMom2->addWidget(labeln);
		layMom2->addWidget(spinn);

    layMom->addLayout(layMom1);
    layMom->addLayout(layMom2);

    grParametersMom->setLayout(layMom);

    QHBoxLayout *layFlags = new QHBoxLayout();
    layFlags->setAlignment(Qt::AlignLeft);
		QLabel *labelWidth = new QLabel(tr("Resonance widths:"));
		comboWidth = new QComboBox();
		comboWidth->addItem(tr("Zero-width"));
		comboWidth->addItem(tr("Breit-Wigner"));
		comboWidth->addItem(tr("eBW"));
		comboWidth->setCurrentIndex(0);
    checkBratio = new QCheckBox(tr("Renormalize branching ratios"));
    checkBratio->setChecked(false);

		checkDecays = new QCheckBox(tr("Perform decays"));
		checkDecays->setChecked(false);

		layFlags->addWidget(labelWidth);
		layFlags->addWidget(comboWidth);
    layFlags->addWidget(checkBratio);
		layFlags->addWidget(checkDecays);
    

		QHBoxLayout *layFile = new QHBoxLayout();
		layFile->setAlignment(Qt::AlignLeft);
		checkFile = new QCheckBox(tr("Write events to file"));
		checkFile->setChecked(false);
		connect(checkFile, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));

		leFilePath = new QLineEdit("");
		leFilePath->setReadOnly(true);
		leFilePath->setText(QDir::currentPath() + "/events.dat");
		leFilePath->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
		buttonChooseFile = new QPushButton(tr("Choose path..."));
		connect(buttonChooseFile, SIGNAL(clicked()), this, SLOT(chooseOutputFile()));
		layFile->addWidget(checkFile);
		layFile->addWidget(leFilePath);
		layFile->addWidget(buttonChooseFile);


    QGroupBox *grParallel = new QGroupBox(tr("Paralellization:"));

    QHBoxLayout *layParallel = new QHBoxLayout();
    layParallel->setAlignment(Qt::AlignLeft);

    checkOMP = new QCheckBox(tr("OpenMP"));
    checkOMP->setChecked(false);
    layParallel->addWidget(checkOMP);

    grParallel->setLayout(layParallel);

    QHBoxLayout *layEvents = new QHBoxLayout();
    layEvents->setAlignment(Qt::AlignLeft);
    QLabel *labelEvents = new QLabel(tr("Events:"));
    spinEvents = new QSpinBox();
    spinEvents->setMinimum(0);
    spinEvents->setMaximum(1000000000);
    spinEvents->setValue(10000);

    layEvents->addWidget(labelEvents);
    layEvents->addWidget(spinEvents);

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
    teDebug->setFontPointSize(10);

		editorLay->addLayout(layModelEnsemble);
    editorLay->addWidget(grEV);
    editorLay->addWidget(grParameters);
    editorLay->addLayout(layQB);
    editorLay->addWidget(grParametersMom);
    editorLay->addLayout(layFlags);
		editorLay->addLayout(layFile);
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

		spinTemperature->setValue(160.0);
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
	ThermalModelConfig ret;

	ret.ModelType = ThermalModelConfig::Ideal;
	if (radGCE->isChecked()) {
		if (radEVD->isChecked())
			ret.ModelType = ThermalModelConfig::DiagonalEV;
		if (radQVDW->isChecked())
			ret.ModelType = ThermalModelConfig::QvdW;

		ret.Ensemble = ThermalModelConfig::EnsembleGCE;
	}

	if (radCE->isChecked()) {
		ret.ModelType = ThermalModelConfig::CE;
		ret.Ensemble = ThermalModelConfig::EnsembleCE;
	}

	if (radSCE->isChecked()) {
		ret.ModelType = ThermalModelConfig::SCE;
		if (radEVD->isChecked())
			ret.ModelType = ThermalModelConfig::EVSCE;
		if (radQVDW->isChecked())
			ret.ModelType = ThermalModelConfig::VDWSCE;

		ret.Ensemble = ThermalModelConfig::EnsembleSCE;
	}

	ret.InteractionModel = ThermalModelConfig::InteractionIdeal;
	if (radEVD->isChecked())
		ret.InteractionModel = ThermalModelConfig::InteractionEVDiagonal;
	if (radQVDW->isChecked())
		ret.InteractionModel = ThermalModelConfig::InteractionQVDW;

	ret.QuantumStatistics = 0;// static_cast<int>(!radioBoltz->isChecked());
	ret.QuantumStatisticsType = 0;// static_cast<int>(CBQuadratures->isChecked());
	ret.QuantumStatisticsInclude = 0;

	if (radioUniform->isChecked())
		ret.Interaction = 0;
	else if (radioBaglike->isChecked())
		ret.Interaction = 1;
	else if (radioMesons->isChecked())
		ret.Interaction = 2;
	else
		ret.Interaction = 3;

	ret.EVRadius = spinRadius->value();
	ret.InteractionInput = strEVPath.toStdString();

	ret.T = spinTemperature->value() * 1.e-3;
	ret.muB = spinmuB->value() * 1.e-3;
	ret.muQ = spinmuQ->value() * 1.e-3;
	ret.muS = spinmuS->value() * 1.e-3;
	ret.muC = spinmuC->value() * 1.e-3;
	ret.gq = spingammaq->value();
	ret.gS = spingammaS->value();
	ret.gC = 1.;
	ret.VolumeR = spinVolumeR->value();
	ret.VolumeRSC = spinVolumeRSC->value();
	ret.B = spinB->value();
	ret.Q = spinQ->value();
	ret.S = spinS->value();
	ret.C = 0;

	ret.QoverB = spinQBRatio->value();
	ret.ConstrainMuQ = checkFixMuQ->isChecked();
	ret.ConstrainMuS = checkFixMuS->isChecked();
	ret.ConstrainMuC = checkFixMuC->isChecked();

	ret.FiniteWidth = comboWidth->currentIndex();
	ret.RenormalizeBR = checkBratio->isChecked();
	ret.ComputeFluctations = false;

	return ret;
}

void EventGeneratorTab::loadEVFromFile()
{
	QString listpathprefix = QString(INPUT_FOLDER) + "/interaction";
	if (strEVPath.size() != 0)
		listpathprefix = QString(strEVPath);
	QString path = QFileDialog::getOpenFileName(this, tr("Open file with EV/vdW parameters"), listpathprefix);
	if (path.length() > 0)
	{
		strEVPath = path;
	}
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

        double tmin = 0., tmax = 0.;
				double tl = 0., tr = 2.;
				if (TPS->PdgToId(spectra->fParticles[id].GetPDGID()) != -1)
					tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;

        plotDistr->graph(0)->clearData();
        plotDistr->graph(0)->setName(paramnames[index]);
        plotDistr->graph(1)->clearData();
        plotDistr->graph(1)->setName(paramnames[index]);

        int ttype = comboDistr->currentIndex();
        if (ttype>=0 && ttype<=2) {
            QVector<double> x1,y1,y1err;
            x1    = QVector<double>::fromStdVector (spectra->fParticles[id].GetXVector(ttype) );
            y1    = QVector<double>::fromStdVector (spectra->fParticles[id].GetYVector(ttype) );
            y1err = QVector<double>::fromStdVector (spectra->fParticles[id].GetYErrorVector(ttype) );
            for(int i=0;i<x1.size();++i) {
                tmin = std::min(tmin, y1[i]);
                tmax = std::max(tmax, y1[i]);
            }
            plotDistr->graph(0)->setDataValueError(x1,y1,y1err);
        }


        plotDistr->xAxis->setLabel(paramnamesx[index]);
        plotDistr->yAxis->setLabel(paramnames[index]);

        QVector<double> x2,y2;
        x2 = QVector<double>::fromStdVector (spectra->fParticles[id].GetModelX() );
        y2 = QVector<double>::fromStdVector (spectra->fParticles[id].GetModelY() );
        for(int i=0;i<x2.size();++i) {
            tmin = std::min(tmin, y2[i]);
            tmax = std::max(tmax, y2[i]);
        }
        plotDistr->graph(1)->setData(x2,y2);

        plotDistr->xAxis->setRange(tl, tr);
        if (x2.size()>1) plotDistr->xAxis->setRange(x2[0]-0.5*(x2[1]-x2[0]), x2[x2.size()-1]+0.5*(x2[1]-x2[0]));

        if (comboDistr->currentIndex()==2) {
            plotDistr->yAxis->setScaleType(QCPAxis::stLogarithmic);
            plotDistr->yAxis->setScaleLogBase(100);
            plotDistr->yAxis->setRange(1.e-1*tmin, 1.e1*tmax);
            plotDistr->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
            plotDistr->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
            plotDistr->yAxis->setSubTickCount(10);
        }
        else {
            plotDistr->yAxis->setScaleType(QCPAxis::stLinear);
            plotDistr->yAxis->setRange(0., 1.1*tmax);
            plotDistr->yAxis->setNumberFormat("g");
            plotDistr->yAxis->setNumberPrecision(5); // makes sure "1*10^4" is displayed only as "10^4"
            plotDistr->yAxis->setSubTickCount(5);
        }

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
        double tmin = 0., tmax = 0.;
        double tl = 0., tr = rightlimit;
				//double tl = 0., tr = TPS->ParticleByPDG(spectra->fParticles[id].GetPDGID()).Mass() + 2.;

        plotDistr->graph(0)->clearData();
        plotDistr->graph(0)->setName(paramnames[index]);
        plotDistr->graph(1)->clearData();
        plotDistr->graph(1)->setName(paramnames[index]);

        for(int i=0;i<x1.size();++i) {
            tmin = std::min(tmin, y1[i]);
            tmax = std::max(tmax, y1[i]);
        }
        plotDistr->graph(0)->setDataValueError(x1,y1,y1err);

        plotDistr->xAxis->setLabel(paramnamesx[index]);
        plotDistr->yAxis->setLabel(paramnames[index]);
        for(int i=0;i<x2.size();++i) {
            tmin = std::min(tmin, y2[i]);
            tmax = std::max(tmax, y2[i]);
        }
        plotDistr->graph(1)->setData(x2,y2);

        if (comboDistr->currentIndex()==2) {
            plotDistr->yAxis->setScaleType(QCPAxis::stLogarithmic);
            plotDistr->yAxis->setScaleLogBase(100);
            plotDistr->yAxis->setRange(1.e-1*tmin, 1.e1*tmax);
            plotDistr->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
            plotDistr->yAxis->setNumberPrecision(0);
            plotDistr->yAxis->setSubTickCount(10);
        }
        else {
            plotDistr->yAxis->setScaleType(QCPAxis::stLinear);
            plotDistr->yAxis->setRange(0., 1.1*tmax);
            plotDistr->yAxis->setNumberFormat("g");
            plotDistr->yAxis->setNumberPrecision(5); // makes sure "1*10^4" is displayed only as "10^4"
            plotDistr->yAxis->setSubTickCount(5);
        }

        plotDistr->xAxis->setRange(tl, tr);
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
    dbgstrm << "CE acceptance rate " << EventGeneratorBase::fCEAccepted / (double)(EventGeneratorBase::fCETotal) << " events" << endl;
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
            plotDistr->graph(0)->setDataValueError(x1,y1,y1err);
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
	bool tfl = !radIdeal->isChecked();

	if (tfl) {
		spinRadius->setEnabled(true);
		radioUniform->setEnabled(true);
		radioMesons->setEnabled(true);
		radioBaglike->setEnabled(true);
		radioCustomEV->setEnabled(true);

		radCE->setEnabled(false);

		if (radCE->isChecked())
			radGCE->setChecked(true);

		radSCE->setEnabled(true);
	}
	else {
		spinRadius->setEnabled(false);
		radioUniform->setEnabled(false);
		radioMesons->setEnabled(false);
		radioBaglike->setEnabled(false);
		radioCustomEV->setEnabled(false);

		radCE->setEnabled(true);
		radSCE->setEnabled(true);
	}

	if (!radCE->isChecked()) {
		spinmuB->setEnabled(true);
		spinmuS->setEnabled(true);
		spinmuQ->setEnabled(true);
		spinmuC->setEnabled(true);
		spinB->setEnabled(false);
		spinS->setEnabled(false);
		spinQ->setEnabled(false);
		//spinQBRatio->setEnabled(checkFixMuQ->isChecked());

		//checkFixMuQ->setEnabled(true);
		//checkFixMuS->setEnabled(true);
		//checkFixMuC->setEnabled(true);


	}
	else {
		spinmuB->setEnabled(false);
		spinmuS->setEnabled(false);
		spinmuQ->setEnabled(false);
		spinmuC->setEnabled(false);
		spinB->setEnabled(true);
		spinS->setEnabled(true);
		spinQ->setEnabled(true);
		//spinQBRatio->setEnabled(false);

		//checkFixMuQ->setEnabled(false);
		//checkFixMuS->setEnabled(false);
		//checkFixMuC->setEnabled(false);
	}


	if (radSCE->isChecked()) {
		//checkFixMuS->setEnabled(false);
		spinmuS->setEnabled(false);
		spinVolumeRSC->setEnabled(true);
	}
	else spinVolumeRSC->setEnabled(false);
	//else spinmuS->setEnabled(true);


	spinmuS->setVisible(model->TPS()->hasStrange());
	spinmuC->setVisible(model->TPS()->hasCharmed());
	labelmuS->setVisible(model->TPS()->hasStrange());
	labelmuC->setVisible(model->TPS()->hasCharmed());


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
    //tableSpectra->reset();
}

void EventGeneratorTab::updateThermalParameters() {
    double ssqrt = sqrt(2.*xMath::mnucleon()*(spinEnergy->value() + 2.*xMath::mnucleon()));
    spinTemperature->setValue(Tss(ssqrt)*1.e3);
    spinmuB->setValue(muBss(ssqrt)*1.e3);
    spingammaS->setValue(gammaSss(ssqrt));
}

void EventGeneratorTab::generateEvents(const ThermalModelConfig & config)
{
		ThermalModelBase *modelnew;

		if (config.Ensemble == ThermalModelConfig::EnsembleGCE) modelnew = new ThermalModelIdeal(model->TPS());
		else if (config.Ensemble == ThermalModelConfig::EnsembleCE) {
			modelnew = new ThermalModelCanonical(model->TPS());
			ThermalModelParameters params;
			params.T = 0.100;
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

		if (model != NULL) delete model;
		model = modelnew;

		model->SetNormBratio(config.RenormalizeBR);

		if (config.FiniteWidth == 0)
			model->SetUseWidth(ThermalParticle::ZeroWidth);
		else if (config.FiniteWidth == 1)
			model->SetUseWidth(ThermalParticle::BWTwoGamma);
		else if (config.FiniteWidth == 2)
			model->SetUseWidth(ThermalParticle::eBW);
		else
			model->SetUseWidth(ThermalParticle::ZeroWidth);

		model->SetStatistics(config.QuantumStatistics);
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
		model->SetStrangenessCanonicalVolumeRadius(config.VolumeRSC);


		ThermalModelBase *modelEVVDW;
		ThermalParticleSystem TPSt = *model->TPS();
		if (config.InteractionModel == ThermalModelConfig::InteractionQVDW)
			modelEVVDW = new ThermalModelVDWFull(&TPSt);
		else if (config.InteractionModel == ThermalModelConfig::InteractionEVDiagonal)
			modelEVVDW = new ThermalModelEVDiagonal(&TPSt);
		else
			modelEVVDW = new ThermalModelIdeal(&TPSt);

		std::vector<double> radii(model->TPS()->Particles().size(), 0.);

		// Uniform EV
		if (config.Interaction == 0) {
			modelEVVDW->SetRadius(config.EVRadius);
			std::fill(radii.begin(), radii.end(), config.EVRadius);
		}

		// Bag model EV
		if (config.Interaction == 1) {
			for (int i = 0; i < modelEVVDW->TPS()->Particles().size(); ++i) {
				ThermalParticle &part = modelEVVDW->TPS()->Particle(i);
				radii[i] = config.EVRadius * pow(part.Mass() / xMath::mnucleon(), 1. / 3.);
			}
			modelEVVDW->FillVirial(radii);
		}

		// Two-component EV
		if (config.Interaction == 2) {
			for (int i = 0; i < modelEVVDW->TPS()->Particles().size(); ++i) {
				ThermalParticle &part = modelEVVDW->TPS()->Particle(i);
				if (part.BaryonCharge() != 0)
					radii[i] = config.EVRadius * pow(abs(part.BaryonCharge()), 1. / 3.);
				else
					radii[i] = 0.;
			}
			modelEVVDW->FillVirial(radii);
		}

		modelEVVDW->FillVirial(radii); // Just in case

		// Read from file
		if (config.Interaction == 3) {
			modelEVVDW->ReadInteractionParameters(config.InteractionInput);
		}


		model->SetQoverB(config.QoverB);
		model->ConstrainMuQ(config.ConstrainMuQ);
		model->ConstrainMuS(config.ConstrainMuS);
		model->ConstrainMuC(config.ConstrainMuC);
		model->FixParameters();
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



		EventGeneratorConfiguration::ModelType EV = EventGeneratorConfiguration::PointParticle;
		if (config.InteractionModel == ThermalModelConfig::InteractionEVDiagonal)
			EV = EventGeneratorConfiguration::DiagonalEV;
		if (config.InteractionModel == ThermalModelConfig::InteractionQVDW)
			EV = EventGeneratorConfiguration::QvdW;


		if (radioSR->isChecked()) generator = new SiemensRasmussenEventGenerator(model, spinTkin->value() * 1.e-3, spinBeta->value(), false, EV, modelEVVDW);
		else generator = new SSHEventGenerator(model, spinTkin->value() * 1.e-3, spinBetat->value(), spinEtaMax->value(), spinn->value(), false, EV, modelEVVDW);

		delete modelEVVDW;

		generator->SetCollisionKineticEnergy(spinEnergy->value());
		
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