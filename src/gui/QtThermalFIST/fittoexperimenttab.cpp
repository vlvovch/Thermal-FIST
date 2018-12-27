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
#include "chi2dialog.h"
#include "chi2ProfileDialog.h"

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
		cpath = QString(INPUT_FOLDER) + "/data";
    dbgstrm.setString(&dbgstr);

    TPS = modelop->TPS();
    model = new ThermalModelIdeal(modelop->TPS());
    *model = *modelop;

		fitcopy = new ThermalModelFit(model);

    QHBoxLayout *DataEditLay = new QHBoxLayout();

    QVBoxLayout *dataLayv = new QVBoxLayout();

    quantities.resize(0);
    quantities.push_back(FittedQuantity(ExperimentRatio(-211, 211, 1.09898, 0.0798674)));
    quantities.push_back(FittedQuantity(ExperimentRatio(-321, 321, 0.324873, 0.027183)));
    quantities.push_back(FittedQuantity(ExperimentRatio(321, 211, 0.201706, 0.0160556)));
    quantities.push_back(FittedQuantity(ExperimentMultiplicity(-211, 300., 30.)));

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
    //connect(tableQuantities->selectionModel(), SIGNAL(selectionChanged(QItemSelection,QItemSelection)), this, SLOT(changedRow()));

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

    QHBoxLayout *layMisc = new QHBoxLayout();
    layMisc->setAlignment(Qt::AlignLeft);

    buttonResults = new QPushButton(tr("Equation of state..."));
    connect(buttonResults, SIGNAL(clicked()), this, SLOT(showResults()));
    buttonChi2Map = new QPushButton(tr("Show chi2 map..."));
    connect(buttonChi2Map, SIGNAL(clicked()), this, SLOT(showChi2Map()));
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

    //dataLayv->addWidget(labelQuantities);
    dataLayv->addLayout(layoutTop);
    dataLayv->addWidget(tableQuantities);
    dataLayv->addLayout(layEditQuantities);
    dataLayv->addWidget(labelParams);
    dataLayv->addWidget(tableParameters);
    dataLayv->addLayout(layMisc);

   // QGroupBox *grEditor = new QGroupBox(tr("Editor"));

    QVBoxLayout *editorLay = new QVBoxLayout();
    editorLay->setContentsMargins(15, 0, 0, 0);
    //editorLay->setMargin(10);
    editorLay->setAlignment(Qt::AlignTop);


		QHBoxLayout *layModelEnsemble = new QHBoxLayout();
		layModelEnsemble->setAlignment(Qt::AlignLeft);

		QGroupBox *grModels = new QGroupBox(tr("HRG Model:"));

		QHBoxLayout *layModels = new QHBoxLayout();
		layModels->setAlignment(Qt::AlignLeft);


		radIdeal = new QRadioButton(tr("Ideal"));
		radEVD = new QRadioButton(tr("Diagonal EV"));
		radEVCRS = new QRadioButton(tr("Crossterms EV"));
		radQVDW = new QRadioButton(tr("QvdW"));

    radIdeal->setToolTip(tr("Point-particle ideal gas"));
    radEVD->setToolTip(tr("Diagonal excluded volume model"));
    radEVCRS->setToolTip(tr("Crossterms (non-diagonal) excluded volume model"));
    radQVDW->setToolTip(tr("Quantum van der Waals HRG model"));

		connect(radIdeal, SIGNAL(clicked()), this, SLOT(modelChanged()));
		connect(radEVD, SIGNAL(clicked()), this, SLOT(modelChanged()));
		connect(radEVCRS, SIGNAL(clicked()), this, SLOT(modelChanged()));
		connect(radQVDW, SIGNAL(clicked()), this, SLOT(modelChanged()));

		layModels->addWidget(radIdeal);
		layModels->addWidget(radEVD);
		layModels->addWidget(radEVCRS);
		layModels->addWidget(radQVDW);

		grModels->setLayout(layModels);


		QGroupBox *grEnsemble = new QGroupBox(tr("Ensemble:"));

		QHBoxLayout *layEnsemble = new QHBoxLayout();
		layEnsemble->setAlignment(Qt::AlignLeft);

		radGCE = new QRadioButton(tr("GCE"));
		radCE = new QRadioButton(tr("CE"));
		radSCE = new QRadioButton(tr("SCE"));

    radGCE->setToolTip(tr("Grand canonical ensemble"));
    radCE->setToolTip(tr("Canonical ensemble"));
    radSCE->setToolTip(tr("Strangeness-canonical ensemble"));

		connect(radGCE, SIGNAL(clicked()), this, SLOT(modelChanged()));
		connect(radCE, SIGNAL(clicked()), this, SLOT(modelChanged()));
		connect(radSCE, SIGNAL(clicked()), this, SLOT(modelChanged()));

		layEnsemble->addWidget(radGCE);
		layEnsemble->addWidget(radCE);
		layEnsemble->addWidget(radSCE);

		grEnsemble->setLayout(layEnsemble);

		layModelEnsemble->addWidget(grModels);
		layModelEnsemble->addWidget(grEnsemble);

		radIdeal->setChecked(true);
		radGCE->setChecked(true);

    QGroupBox *grStats = new QGroupBox(tr("Statistics:"));

    QHBoxLayout *layStats = new QHBoxLayout();
    layStats->setAlignment(Qt::AlignLeft);
    radioBoltz  = new QRadioButton(tr("Boltzmann"));
    radioQuant  = new QRadioButton(tr("Quantum"));

    radioBoltz->setToolTip(tr("Maxwell-Boltzmann"));
    radioQuant->setToolTip(tr("Fermi-Dirac/Bose-Einstein"));

		CBBoseOnly  = new QCheckBox(tr("Mesons only"));
		CBPionsOnly = new QCheckBox(tr("Pions only"));
		CBQuadratures = new QCheckBox(tr("Use quadratures"));
		CBQuadratures->setChecked(true);

    CBBoseOnly->setToolTip(tr("Include quantum statistics for mesons only"));
    CBPionsOnly->setToolTip(tr("Include quantum statistics for pions only"));
    CBQuadratures->setToolTip(tr("Use Gauss-Laguerre quadratures (slower but more reliable) or cluster expansion (faster but unreliable for large mu)"));

    connect(radioBoltz, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));
		//connect(radioQuant, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));
		//connect(CBBoseOnly, SIGNAL(stateChanged(bool)), this, SLOT(modelChanged()));
		//connect(CBPionsOnly, SIGNAL(stateChanged(bool)), this, SLOT(modelChanged()));

    layStats->addWidget(radioBoltz);
    layStats->addWidget(radioQuant);
		layStats->addWidget(CBBoseOnly);
		layStats->addWidget(CBPionsOnly);
		layStats->addWidget(CBQuadratures);

    grStats->setLayout(layStats);

    radioQuant->setChecked(true);

    QGroupBox *grEV = new QGroupBox(tr("Excluded volume/van der Waals:"));

    QHBoxLayout *layRadius = new QHBoxLayout();
    layRadius->setAlignment(Qt::AlignLeft);
    QLabel *labelRadius = new QLabel(tr("Radius (fm):"));
    spinRadius = new QDoubleSpinBox();
    spinRadius->setMinimum(0.);
    spinRadius->setMaximum(100.);
    spinRadius->setValue(0.3);
    layRadius->addWidget(labelRadius);
    layRadius->addWidget(spinRadius);
    radioUniform = new QRadioButton(tr("Same for all"));
    radioUniform->setToolTip(tr("Same eigenvolume for all particles"));
    radioBaglike = new QRadioButton(tr("Bag-like"));
    radioBaglike->setToolTip(tr("Eigenvolumes scale linearly with mass. The input radius fixes the radius parameter of protons"));
    radioMesons = new QRadioButton(tr("Point-like mesons"));
    radioMesons->setToolTip(tr("Eigenvolumes scale linearly with absolute baryon number. The input radius fixes the radius parameter of baryons"));
    radioCustomEV = new QRadioButton(tr("Custom..."));
    radioCustomEV->setToolTip(tr("Load EV/QvdW parameters for different (pairs of) particles from file"));
    strEVPath = "";
		connect(radioCustomEV, SIGNAL(clicked()), this, SLOT(loadEVFromFile()));
    layRadius->addWidget(radioUniform);
    layRadius->addWidget(radioBaglike);
    layRadius->addWidget(radioMesons);
		layRadius->addWidget(radioCustomEV);
    radioUniform->setChecked(true);

    grEV->setLayout(layRadius);


    QGroupBox *grParameters = new QGroupBox(tr("Parameters:"));

    QGridLayout *layParameters = new QGridLayout();
    layParameters->setAlignment(Qt::AlignLeft);
    //layQuant->setHorizontalSpacing(15);
    QLabel *labelTemperature = new QLabel(tr("T (MeV):"));
    spinTemperature = new QDoubleSpinBox();
    spinTemperature->setMinimum(20.);
    spinTemperature->setMaximum(2000.);
    spinTemperature->setValue(model->Parameters().T * 1e3);
    QLabel *labelTmin = new QLabel(tr("T<sub>min</sub> (MeV):"));
    spinTmin = new QDoubleSpinBox();
    spinTmin->setMinimum(20.);
    spinTmin->setMaximum(2000.);
    spinTmin->setValue(20.);
    QLabel *labelTmax = new QLabel(tr("T<sub>max</sub> (MeV):"));
    spinTmax = new QDoubleSpinBox();
    spinTmax->setMinimum(20.);
    spinTmax->setMaximum(2000.);
    spinTmax->setValue(500.);
    CBTemperature = new QCheckBox(tr("Fit"));
    CBTemperature->setChecked(true);
    QLabel *labelmuB = new QLabel(tr("μB (MeV):"));
    spinmuB = new QDoubleSpinBox();
    spinmuB->setMinimum(-3000.);
    spinmuB->setMaximum(3000.);
    spinmuB->setValue(model->Parameters().muB * 1e3);
    QLabel *labelmuBmin = new QLabel(tr("μB<sub>min</sub> (MeV):"));
    spinmuBmin = new QDoubleSpinBox();
    spinmuBmin->setMinimum(-3000.);
    spinmuBmin->setMaximum(3000.);
    spinmuBmin->setValue(0.);
    QLabel *labelmuBmax = new QLabel(tr("μB<sub>max</sub> (MeV):"));
    spinmuBmax = new QDoubleSpinBox();
    spinmuBmax->setMinimum(-3000.);
    spinmuBmax->setMaximum(3000.);
    spinmuBmax->setValue(950.);
    CBmuB = new QCheckBox(tr("Fit"));
    CBmuB->setChecked(true);
    QLabel *labelgammaq = new QLabel(tr("γq:"));
    spingammaq = new QDoubleSpinBox();
    spingammaq->setMinimum(0.);
    spingammaq->setMaximum(100.);
		spingammaq->setDecimals(5);
    spingammaq->setValue(model->Parameters().gammaq);
    QLabel *labelgqmin = new QLabel(tr("γq<sub>min</sub>:"));
    spingqmin = new QDoubleSpinBox();
    spingqmin->setMinimum(0.);
    spingqmin->setMaximum(100.);
		spingqmin->setDecimals(5);
    spingqmin->setValue(0.);
    QLabel *labelgqmax = new QLabel(tr("γq<sub>max</sub>:"));
    spingqmax = new QDoubleSpinBox();
    spingqmax->setMinimum(0.);
    spingqmax->setMaximum(100.);
		spingqmax->setDecimals(5);
    spingqmax->setValue(3.);
    CBgammaq = new QCheckBox(tr("Fit"));
    CBgammaq->setChecked(false);
    QLabel *labelgammaS = new QLabel(tr("γS:"));
    spingammaS = new QDoubleSpinBox();
    spingammaS->setMinimum(0.);
		spingammaS->setDecimals(5);
    spingammaS->setMaximum(100.);
    spingammaS->setValue(model->Parameters().gammaS);
    QLabel *labelgsmin = new QLabel(tr("γS<sub>min</sub>:"));
    spingsmin = new QDoubleSpinBox();
    spingsmin->setMinimum(0.);
    spingsmin->setMaximum(100.);
		spingsmin->setDecimals(5);
    spingsmin->setValue(0.);
    QLabel *labelgsmax = new QLabel(tr("γS<sub>max</sub>:"));
    spingsmax = new QDoubleSpinBox();
    spingsmax->setMinimum(0.);
    spingsmax->setMaximum(100.);
		spingsmax->setDecimals(5);
    spingsmax->setValue(3.);
    CBgammaS = new QCheckBox(tr("Fit"));
    CBgammaS->setChecked(false);

    QLabel *labelVolume = new QLabel(tr("R (fm):"));
    spinVolumeR = new QDoubleSpinBox();
    spinVolumeR->setMinimum(0.);
    spinVolumeR->setMaximum(1.e4);
    spinVolumeR->setValue(8.);
    spinVolumeR->setDecimals(5);
    QLabel *labelVRmin = new QLabel(tr("R<sub>min</sub> (fm):"));
    spinVRmin = new QDoubleSpinBox();
    spinVRmin->setMinimum(0.);
    spinVRmin->setMaximum(1.e4);
    spinVRmin->setValue(0.);
    spinVRmin->setDecimals(5);
    QLabel *labelVRmax = new QLabel(tr("R<sub>max</sub> (fm):"));
    spinVRmax = new QDoubleSpinBox();
    spinVRmax->setMinimum(0.);
    spinVRmax->setMaximum(1.e4);
    spinVRmax->setValue(25.);
    spinVRmax->setDecimals(5);
    CBVolumeR = new QCheckBox(tr("Fit"));
    CBVolumeR->setChecked(true);

    layParameters->addWidget(labelTemperature, 0, 0);
    layParameters->addWidget(spinTemperature, 0, 1);
    layParameters->addWidget(labelTmin, 0, 2);
    layParameters->addWidget(spinTmin, 0, 3);
    layParameters->addWidget(labelTmax, 0, 4);
    layParameters->addWidget(spinTmax, 0, 5);
    layParameters->addWidget(CBTemperature, 0, 6);
    layParameters->addWidget(labelmuB, 1, 0);
    layParameters->addWidget(spinmuB, 1, 1);
    layParameters->addWidget(labelmuBmin, 1, 2);
    layParameters->addWidget(spinmuBmin, 1, 3);
    layParameters->addWidget(labelmuBmax, 1, 4);
    layParameters->addWidget(spinmuBmax, 1, 5);
    layParameters->addWidget(CBmuB, 1, 6);
    layParameters->addWidget(labelgammaq, 2, 0);
    layParameters->addWidget(spingammaq, 2, 1);
    layParameters->addWidget(labelgqmin, 2, 2);
    layParameters->addWidget(spingqmin, 2, 3);
    layParameters->addWidget(labelgqmax, 2, 4);
    layParameters->addWidget(spingqmax, 2, 5);
    layParameters->addWidget(CBgammaq, 2, 6);
    layParameters->addWidget(labelgammaS, 3, 0);
    layParameters->addWidget(spingammaS, 3, 1);
    layParameters->addWidget(labelgsmin, 3, 2);
    layParameters->addWidget(spingsmin, 3, 3);
    layParameters->addWidget(labelgsmax, 3, 4);
    layParameters->addWidget(spingsmax, 3, 5);
    layParameters->addWidget(CBgammaS, 3, 6);


    layParameters->addWidget(labelVolume, 4, 0);
    layParameters->addWidget(spinVolumeR, 4, 1);
    layParameters->addWidget(labelVRmin, 4, 2);
    layParameters->addWidget(spinVRmin, 4, 3);
    layParameters->addWidget(labelVRmax, 4, 4);
    layParameters->addWidget(spinVRmax, 4, 5);
    layParameters->addWidget(CBVolumeR, 4, 6);

    grParameters->setLayout(layParameters);


    QLabel *labelQB = new QLabel(tr("Q/B ratio:"));
    spinQBRatio = new QDoubleSpinBox();
    spinQBRatio->setDecimals(3);
    spinQBRatio->setMinimum(0.);
    spinQBRatio->setValue(model->QoverB());

		checkFixMuQ = new QCheckBox(tr("Constrain μQ"));
		checkFixMuQ->setChecked(true);
    checkFixMuQ->setToolTip(tr("Constrain μQ to reproduce the needed Q/B ratio, otherwise fit μQ as free parameter"));
		connect(checkFixMuQ, SIGNAL(clicked()), this, SLOT(modelChanged()));

		checkFixMuS = new QCheckBox(tr("Constrain μS"));
		checkFixMuS->setChecked(true);
    checkFixMuS->setToolTip(tr("Constrain μS to reproduce the needed Q/B ratio, otherwise fit μS as free parameter"));


    QHBoxLayout *layCE = new QHBoxLayout();
    layCE->setAlignment(Qt::AlignLeft);
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

    spinB->setToolTip(tr("Total baryon number in CE calculation"));
    spinQ->setToolTip(tr("Total electric charge in CE calculation"));
    spinS->setToolTip(tr("Total strangeness in CE calculation"));

    layCE->addWidget(labelQB);
    layCE->addWidget(spinQBRatio);
		layCE->addWidget(checkFixMuQ);
		layCE->addWidget(checkFixMuS);
    layCE->addSpacing(20);
    layCE->addWidget(labelB);
    layCE->addWidget(spinB);
    layCE->addWidget(labelS);
    layCE->addWidget(spinS);
    layCE->addWidget(labelQ);
    layCE->addWidget(spinQ);

    QHBoxLayout *layFlags = new QHBoxLayout();
    layFlags->setAlignment(Qt::AlignLeft);

		QLabel *labelWidth = new QLabel(tr("Resonance widths:"));
		comboWidth = new QComboBox();
		comboWidth->addItem(tr("Zero-width"));
		comboWidth->addItem(tr("Const Breit-Wigner"));
		comboWidth->addItem(tr("eBW"));
    comboWidth->addItem(tr("eBW (const BRs)"));
		comboWidth->setCurrentIndex(static_cast<int>(model->TPS()->ResonanceWidthIntegrationType()));
    comboWidth->setToolTip(tr("Prescription for treatment of resonance widths"));
    checkBratio = new QCheckBox(tr("Renormalize branching ratios"));
    checkBratio->setChecked(false);
    checkBratio->setToolTip(tr("Renormalize branching ratios of all particle to sum to 100\%"));

	checkFitRc  = new QCheckBox(tr("Fit (Str.-)Can. radius"));
	checkFitRc->setChecked(true);
  checkFitRc->setToolTip(tr("Fit the radius of correlation volume, otherwise assume correlation volume = system volume"));

	layFlags->addWidget(labelWidth);
	layFlags->addWidget(comboWidth);
  layFlags->addWidget(checkBratio);
	layFlags->addWidget(checkFitRc);

    QGroupBox *grParallel = new QGroupBox(tr("Paralellization:"));

    QHBoxLayout *layParallel = new QHBoxLayout();
    layParallel->setAlignment(Qt::AlignLeft);

    checkOMP = new QCheckBox(tr("OpenMP"));
    checkOMP->setChecked(false);
    layParallel->addWidget(checkOMP);

    grParallel->setLayout(layParallel);

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
    teDebug->setFontPointSize(10);

		editorLay->addLayout(layModelEnsemble);
    editorLay->addWidget(grStats);
    editorLay->addWidget(grEV);
    editorLay->addWidget(grParameters);
    editorLay->addLayout(layCE);
    editorLay->addLayout(layFlags);
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

		QString datapathprefix = QString(INPUT_FOLDER) + "/data";
		quantities = ThermalModelFit::loadExpDataFromFile((QString(INPUT_FOLDER) + "/data/ALICE-PbPb2.76TeV-0-5-1512.08046.dat").toStdString());
    myModel->setQuantities(&quantities);
    tableQuantities->resizeColumnsToContents();

    lastconfig = getConfig();
}



FitToExperimentTab::~FitToExperimentTab()
{
    if (fitcopy!=NULL) { delete fitcopy; }
    delete model;
}

ThermalModelConfig FitToExperimentTab::getConfig()
{
	ThermalModelConfig ret;

	ret.ModelType = ThermalModelConfig::Ideal;
	if (radGCE->isChecked()) {
		if (radEVD->isChecked())
			ret.ModelType = ThermalModelConfig::DiagonalEV;
		if (radEVCRS->isChecked())
			ret.ModelType = ThermalModelConfig::CrosstermsEV;
		if (radQVDW->isChecked())
			ret.ModelType = ThermalModelConfig::QvdW;
	}

	if (radCE->isChecked()) {
		ret.ModelType = ThermalModelConfig::CE;
	}

	if (radSCE->isChecked()) {
		ret.ModelType = ThermalModelConfig::SCE;
		if (radEVD->isChecked())
			ret.ModelType = ThermalModelConfig::EVSCE;
		if (radQVDW->isChecked())
			ret.ModelType = ThermalModelConfig::VDWSCE;
	}

	ret.QuantumStatistics = static_cast<int>(!radioBoltz->isChecked());
	ret.QuantumStatisticsType = static_cast<int>(CBQuadratures->isChecked());
	ret.QuantumStatisticsInclude = 0;
	if (CBBoseOnly->isChecked())
		ret.QuantumStatisticsInclude = 1;
	if (CBPionsOnly->isChecked())
		ret.QuantumStatisticsInclude = 2;

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
	ret.muQ = 0.;// spinmuQ->value() * 1.e-3;
	ret.muS = 0.;// spinmuS->value() * 1.e-3;
	ret.muC = 0.;// spinmuC->value() * 1.e-3;
	ret.gq = spingammaq->value();
	ret.gS = spingammaS->value();
	ret.gC = 1.;
	ret.VolumeR = spinVolumeR->value();
	ret.VolumeRSC = spinVolumeR->value();//spinVolumeRSC->value();
	ret.B = spinB->value();
	ret.Q = spinQ->value();
	ret.S = spinS->value();
	ret.C = 0;

	ret.QoverB = spinQBRatio->value();
	ret.ConstrainMuQ = checkFixMuQ->isChecked();
	ret.ConstrainMuS = checkFixMuS->isChecked();
	ret.ConstrainMuC = true;// checkFixMuC->isChecked();

	ret.FiniteWidth = comboWidth->currentIndex();
	ret.RenormalizeBR = checkBratio->isChecked();
	ret.ComputeFluctations = false; //checkFluctuations->isChecked();

	return ret;
}

ThermalModelFitParameters FitToExperimentTab::getFitParameters()
{
	ThermalModelFitParameters ret;

	ret.T = FitParameter("T", CBTemperature->isChecked(), spinTemperature->value() * 1.e-3, 0.05,
		spinTmin->value() * 1.e-3, spinTmax->value() * 1.e-3);
	ret.muB = FitParameter("muB", CBmuB->isChecked(), spinmuB->value() * 1.e-3, 0.05,
		spinmuBmin->value() * 1.e-3, spinmuBmax->value() * 1.e-3);
	ret.gammaq = FitParameter("gammaq", CBgammaq->isChecked(), spingammaq->value(), 0.5, 
		spingqmin->value(), spingqmax->value());
	ret.gammaS = FitParameter("gammaS", CBgammaS->isChecked(), spingammaS->value(), 0.5,
		spingsmin->value(), spingsmax->value());
	ret.R = FitParameter("R", CBVolumeR->isChecked(), spinVolumeR->value(), 1.0, 
		spinVRmin->value(), spinVRmax->value());

	ret.muQ.toFit = !checkFixMuQ->isChecked();
	ret.muS.toFit = !checkFixMuS->isChecked();
	ret.muC.toFit = false;
	ret.Rc.toFit  = checkFitRc->isChecked() && (radSCE->isChecked() || radCE->isChecked());
	if (!ret.Rc.toFit)
		ret.Rc.value = ret.R.value;

	return ret;
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
	model->SetStrangenessCanonicalVolumeRadius(config.VolumeRSC);

	model->SetBaryonCharge(config.B);
	model->SetElectricCharge(config.Q);
	model->SetStrangeness(config.S);
	model->SetCharm(config.C);

	timer.start();

	ThermalModelFit *fit = new ThermalModelFit(model);

	fit->SetParameters(params);

	if (config.FiniteWidth == 0)
		model->SetUseWidth(ThermalParticle::ZeroWidth);
	else if (config.FiniteWidth == 1)
		model->SetUseWidth(ThermalParticle::BWTwoGamma);
	else if (config.FiniteWidth == 2)
		model->SetUseWidth(ThermalParticle::eBW);
  else if (config.FiniteWidth == 3)
    model->SetUseWidth(ThermalParticle::eBWconstBR);
	else
		model->SetUseWidth(ThermalParticle::ZeroWidth);

	model->SetResonanceWidthShape(ThermalParticle::RelativisticBreitWiger);
	//model->SetResonanceWidthShape(ThermalParticle::NonRelativisticBreitWiger);

	fit->SetQBConstraint(config.QoverB);

	model->ConstrainMuQ(config.ConstrainMuQ);
	model->ConstrainMuS(config.ConstrainMuS);
	model->ConstrainMuC(config.ConstrainMuC);

	model->SetNormBratio(config.RenormalizeBR);

	if (config.QuantumStatisticsType)
		model->SetCalculationType(IdealGasFunctions::Quadratures);
	else
		model->SetCalculationType(IdealGasFunctions::ClusterExpansion);

	model->SetStatistics(config.QuantumStatistics);

	if (config.QuantumStatisticsInclude == 1 || config.QuantumStatisticsInclude == 2) {
		for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
			ThermalParticle &part = model->TPS()->Particle(i);
			if (config.QuantumStatisticsInclude == 2) {
				if (part.PdgId() != 211 && part.PdgId() != 111 && part.PdgId() != -211) {
					part.UseStatistics(false);
				}
			}
			else if (config.QuantumStatisticsInclude == 1) {
				if (part.BaryonCharge() != 0) {
					part.UseStatistics(false);
				}
			}
		}
	}

	std::vector<double> radii(model->TPS()->Particles().size(), 0.);

	// Uniform EV
	if (config.Interaction == 0) {
		model->SetRadius(config.EVRadius);
		std::fill(radii.begin(), radii.end(), config.EVRadius);
	}

	// Bag model EV
	if (config.Interaction == 1) {
		for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
			ThermalParticle &part = model->TPS()->Particle(i);
			radii[i] = config.EVRadius * pow(part.Mass() / xMath::mnucleon(), 1. / 3.);
		}
		model->FillVirial(radii);
	}

	// Two-component EV
	if (config.Interaction == 2) {
		for (int i = 0; i < model->TPS()->Particles().size(); ++i) {
			ThermalParticle &part = model->TPS()->Particle(i);
			if (part.BaryonCharge() != 0)
				radii[i] = config.EVRadius * pow(abs(part.BaryonCharge()), 1. / 3.);
			else
				radii[i] = 0.;
		}
		model->FillVirial(radii);
	}

	model->FillVirial(radii); // Just in case

  // Read from file
	if (config.Interaction == 3) {
		model->ReadInteractionParameters(config.InteractionInput);
	}

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

void FitToExperimentTab::showChi2Map() {
  if (fitcopy!=NULL) {
      chi2Dialog dialog(this, fitcopy);
      dialog.setWindowFlags(Qt::Window);
      dialog.setMinimumSize(QSize(800, 400));
      dialog.exec();
  }
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
    qDebug() << selectedList.count() << "\n";
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

void FitToExperimentTab::loadEVFromFile()
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
    if (!radIdeal->isChecked()) {
        spinRadius->setEnabled(true);
        radioUniform->setEnabled(true);
        radioMesons->setEnabled(true);
        radioBaglike->setEnabled(true);
				radioCustomEV->setEnabled(true);

				radCE->setEnabled(false);
				if (radCE->isChecked())
					radGCE->setChecked(true);

				if (radEVCRS->isChecked()) {
					radSCE->setEnabled(false);
					if (radSCE->isChecked())
						radGCE->setChecked(true);
				}
				else {
					radSCE->setEnabled(true);
				}
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
        CBmuB->setEnabled(true);
        spinB->setEnabled(false);
        spinS->setEnabled(false);
        spinQ->setEnabled(false);
        spinQBRatio->setEnabled(checkFixMuQ->isChecked());
        radioBoltz->setEnabled(true);
        radioQuant->setEnabled(true);

				checkFixMuQ->setEnabled(true);
				checkFixMuS->setEnabled(true);

				CBQuadratures->setEnabled(true);
    }
    else {
        spinmuB->setEnabled(false);
        CBmuB->setEnabled(false);
        spinB->setEnabled(true);
        spinS->setEnabled(true);
        spinQ->setEnabled(true);
        spinQBRatio->setEnabled(false);

				checkFixMuQ->setEnabled(false);
				checkFixMuS->setEnabled(false);
				//checkFixMuC->setEnabled(false);

				CBQuadratures->setChecked(false);
				CBQuadratures->setEnabled(false);
    }

		if (radSCE->isChecked()) {
			checkFitRc->setEnabled(true);
			checkFixMuS->setEnabled(false);
		}
		else if (radCE->isChecked())
			checkFitRc->setEnabled(true);
		else
			checkFitRc->setEnabled(false);

	checkFixMuS->setVisible(model->TPS()->hasStrange());

	if (radioBoltz->isChecked()) {
		CBBoseOnly->setEnabled(false);
		CBPionsOnly->setEnabled(false);
		CBQuadratures->setEnabled(false);
	}
	else {
		CBBoseOnly->setEnabled(true);
		CBPionsOnly->setEnabled(true);
		CBQuadratures->setEnabled(true);
	}
}

void FitToExperimentTab::resetTPS() {
    model->ChangeTPS(model->TPS());
    myModel->updateAll();

		labelValid->setVisible(false);
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
  dbgstrm << "T\t= " << model->Parameters().T * 1.e3 << " MeV" << endl;
	if (!(model->Ensemble() == ThermalModelBase::CE)) {
	dbgstrm << "muB\t= " << model->Parameters().muB * 1.e3 << " MeV" << endl;
      if (!(model->Ensemble() == ThermalModelBase::SCE))
				dbgstrm << "muS\t= " << model->Parameters().muS * 1.e3 << " MeV" << endl;
      dbgstrm << "muQ\t= " << model->Parameters().muQ * 1.e3 << " MeV" << endl;
  }
  else {
      dbgstrm << "B\t= " << fitcopy->BT() << endl;
      dbgstrm << "S\t= " << fitcopy->ST() << endl;
      dbgstrm << "Q\t= " << fitcopy->QT() << endl;
			dbgstrm << "C\t= " << fitcopy->CT() << endl;
  }
	if (fitcopy->Parameters().gammaq.toFit == true)
		dbgstrm << "gammaq\t= " << model->Parameters().gammaq << endl;
	if (fitcopy->Parameters().gammaS.toFit == true)
		dbgstrm << "gammaS\t= " << model->Parameters().gammaS << endl;
  dbgstrm << "V\t= " << model->Parameters().V << " fm^3" << endl;
	if (model->Ensemble() == ThermalModelBase::SCE || model->Ensemble() == ThermalModelBase::CE)
		dbgstrm << "Vc\t= " << model->Parameters().SVc << " fm^3" << endl;
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

    dbgstrm << "T\t= " << model->Parameters().T * 1.e3 << " MeV" << endl;
    dbgstrm << "muB\t= " << model->Parameters().muB * 1.e3 << " MeV" << endl;
		if (!(model->Ensemble() == ThermalModelBase::CE)) {
			if (!(model->Ensemble() == ThermalModelBase::SCE))
				dbgstrm << "muS\t= " << model->Parameters().muS * 1.e3 << " MeV" << endl;
      dbgstrm << "muQ\t= " << model->Parameters().muQ * 1.e3 << " MeV" << endl;
      if (model->TPS()->hasCharmed() && !(model->Ensemble() == ThermalModelBase::SCE || model->Ensemble() == ThermalModelBase::CCE))
        dbgstrm << "muC\t= " << model->Parameters().muC * 1.e3 << " MeV" << endl;
    }
    else {
        dbgstrm << "B\t= " << model->CalculateBaryonDensity()      * model->Volume() << endl;
        dbgstrm << "S\t= " << model->CalculateStrangenessDensity() * model->Volume() << endl;
        dbgstrm << "Q\t= " << model->CalculateChargeDensity()      * model->Volume() << endl;
				dbgstrm << "C\t= " << model->CalculateCharmDensity()       * model->Volume() << endl;
    }
    dbgstrm << "gammaq\t= " << model->Parameters().gammaq << endl;
    dbgstrm << "gammaS\t= " << model->Parameters().gammaS << endl;
    dbgstrm << "gammaC\t= " << model->Parameters().gammaC << endl;
    dbgstrm << "V\t= " << model->Parameters().V << " fm^3" << endl;
    dbgstrm << endl;
    dbgstrm << "Particle density\t= " << model->CalculateHadronDensity() << " fm^-3" << endl;
    dbgstrm << "Net baryon density\t= " << model->CalculateBaryonDensity() << " fm^-3" << endl;
    dbgstrm << "Net baryon number\t= " << model->CalculateBaryonDensity() * model->Parameters().V << endl;
    if (model->TPS()->hasCharged())
      dbgstrm << "Net electric charge\t= " << model->CalculateChargeDensity() * model->Parameters().V << endl;
    if (model->TPS()->hasStrange())
      dbgstrm << "Net strangeness\t= " << model->CalculateStrangenessDensity() * model->Parameters().V << endl;
    if (model->TPS()->hasCharmed())
      dbgstrm << "Net charm\t= " << model->CalculateCharmDensity() * model->Parameters().V << endl;
    dbgstrm << "E/N\t\t= " << model->CalculateEnergyDensity() / model->CalculateHadronDensity() << endl;
    dbgstrm << "S/Nb\t\t= " << model->CalculateEntropyDensity() / model->CalculateBaryonDensity() << endl;
    dbgstrm << "EV/V\t\t= " << model->CalculateEigenvolumeFraction() << endl;
    dbgstrm << "Q/B\t\t= " << model->CalculateChargeDensity() / model->CalculateBaryonDensity() << endl;
    if (model->TPS()->hasStrange())
      dbgstrm << "S/|S|\t\t= " << model->CalculateStrangenessDensity() / model->CalculateAbsoluteStrangenessDensity() << endl;
    if (model->TPS()->hasCharmed())
      dbgstrm << "C/|C|\t\t= " << model->CalculateCharmDensity() / model->CalculateAbsoluteCharmDensity() << endl;
    dbgstrm << endl;
    dbgstrm << "chi2/ndf\t\t= " << result.chi2ndf * fitcopy->Ndf() << "/" << fitcopy->Ndf() << " = " << result.chi2ndf << endl;
    dbgstrm << endl;
    dbgstrm << "Calculation time = " << timer.elapsed() << " ms" << endl;
    dbgstrm << "----------------------------------------------------------" << endl;
    teDebug->clear();
    teDebug->append(dbgstr);
    dbgstr.clear();

    myModel->updateAll();

		tableQuantities->resizeColumnsToContents();

    tableParameters->clearContents();

    tableParameters->setColumnCount(3);
    int RowCount = 17;

		if (model->Ensemble() == ThermalModelBase::CE)
			RowCount -= 3;

		if (model->Ensemble() == ThermalModelBase::SCE)
			RowCount -= 1;

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

	if (!(model->Ensemble() == ThermalModelBase::CE)) {
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
    
	cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("γs"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.gammaS.value)));
	if (result.gammaS.toFit==true) 
		tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.gammaS.error)));
	else
		tableParameters->setItem(cRow, 2, new QTableWidgetItem("--"));

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

	if (model->Ensemble() == ThermalModelBase::SCE || model->Ensemble() == ThermalModelBase::CE) {
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

    cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("chi2/dof"));
	tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.chi2) + "/" + QString::number(result.ndf)));
    tableParameters->setItem(cRow, 2, new QTableWidgetItem(""));

	cRow++;
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("chi2/dof"));
	tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.chi2ndf)));
    tableParameters->setItem(cRow, 2, new QTableWidgetItem(""));

		if (!(model->Ensemble() == ThermalModelBase::CE)
			&& !(model->Ensemble() == ThermalModelBase::SCE)) {
        cRow++;
        tableParameters->setItem(cRow, 0, new QTableWidgetItem("μS (MeV)"));
        tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.muS.value*1.e3)));
        tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.muS.error*1.e3)));
    }
		if (!(model->Ensemble() == ThermalModelBase::CE)) {
        cRow++;
        tableParameters->setItem(cRow, 0, new QTableWidgetItem("μQ (MeV)"));
        tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(result.muQ.value*1.e3)));
        tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(result.muQ.error*1.e3)));
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
    tableParameters->setItem(cRow, 0, new QTableWidgetItem("P (MeV/fm^3)"));
    tableParameters->setItem(cRow, 1, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().pressure.value*1.e3)));
    tableParameters->setItem(cRow, 2, new QTableWidgetItem(QString::number(fitcopy->ExtendedParameters().pressure.error*1.e3)));

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

