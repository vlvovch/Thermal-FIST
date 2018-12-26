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

#include "DebugText.h"
#include "tablemodel.h"
#include "particledialog.h"
#include "resultdialog.h"
#include "fittoexperimenttab.h"

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

    labelHint = new QLabel(tr("Hint: double-click on particle for more info"));
    QFont tmpf = QApplication::font();
    tmpf.setPointSize(9);
    labelHint->setFont(tmpf);

		labelValid = new QPushButton(tr("Calculation valid!"));
		labelValid->setFlat(true);
		labelValid->setVisible(false);
		connect(labelValid, SIGNAL(clicked()), this, SLOT(showValidityCheckLog()));


    layMisc->addWidget(checkOnlyStable);
    layMisc->addWidget(buttonResults);
    layMisc->addWidget(labelHint);
		layMisc->addStretch(1);
		layMisc->addWidget(labelValid, 0, Qt::AlignRight);

    dataLayv->addWidget(tableParticles);
    dataLayv->addLayout(layMisc);


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
		radEVD   = new QRadioButton(tr("Diagonal EV"));
		radEVCRS = new QRadioButton(tr("Crossterms EV"));
		radQVDW  = new QRadioButton(tr("QvdW"));

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
		radCE  = new QRadioButton(tr("CE"));
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
    radioBoltz = new QRadioButton(tr("Boltzmann"));
    radioQuant = new QRadioButton(tr("Quantum"));

    radioBoltz->setToolTip(tr("Maxwell-Boltzmann"));
    radioQuant->setToolTip(tr("Fermi-Dirac/Bose-Einstein"));

		CBBoseOnly = new QCheckBox(tr("Mesons only"));
		CBPionsOnly = new QCheckBox(tr("Pions only"));
		CBQuadratures = new QCheckBox(tr("Use quadratures"));

    CBBoseOnly->setToolTip(tr("Include quantum statistics for mesons only"));
    CBPionsOnly->setToolTip(tr("Include quantum statistics for pions only"));
    CBQuadratures->setToolTip(tr("Use Gauss-Laguerre quadratures (slower but more reliable) or cluster expansion (faster but unreliable for large mu)"));

		connect(radioBoltz, SIGNAL(toggled(bool)), this, SLOT(modelChanged()));

    layStats->addWidget(radioBoltz);
    layStats->addWidget(radioQuant);
		layStats->addWidget(CBBoseOnly);
		layStats->addWidget(CBPionsOnly);
		layStats->addWidget(CBQuadratures);

    grStats->setLayout(layStats);

    radioQuant->setChecked(true);
		CBQuadratures->setChecked(true);

    QGroupBox *grEV = new QGroupBox(tr("Excluded volume/van der Waals:"));

    QHBoxLayout *layRadius = new QHBoxLayout();
    layRadius->setAlignment(Qt::AlignLeft);
    QLabel *labelRadius = new QLabel(tr("Radius (fm):"));
    spinRadius = new QDoubleSpinBox();
    spinRadius->setMinimum(0.);
    spinRadius->setMaximum(100.);
    spinRadius->setValue(0.3);
    radioUniform  = new QRadioButton(tr("Same for all"));
    radioUniform->setToolTip(tr("Same eigenvolume for all particles"));
    radioBaglike  = new QRadioButton(tr("Bag-like"));
    radioBaglike->setToolTip(tr("Eigenvolumes scale linearly with mass. The input radius fixes the radius parameter of protons"));
    radioMesons   = new QRadioButton(tr("Point-like mesons"));
    radioMesons->setToolTip(tr("Eigenvolumes scale linearly with absolute baryon number. The input radius fixes the radius parameter of baryons"));
		radioCustomEV = new QRadioButton(tr("Custom..."));
    radioCustomEV->setToolTip(tr("Load EV/QvdW parameters for different (pairs of) particles from file"));
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


    QGroupBox *grParameters = new QGroupBox(tr("Parameters:"));

    QGridLayout *layParameters = new QGridLayout();
    layParameters->setAlignment(Qt::AlignLeft);
    QLabel *labelTemperature = new QLabel(tr("T (MeV):"));
    spinTemperature = new QDoubleSpinBox();
    spinTemperature->setMinimum(1.);
    spinTemperature->setMaximum(10000.);
    spinTemperature->setValue(model->Parameters().T * 1e3);
    spinTemperature->setToolTip(tr("Temperature"));
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
		QLabel *labelgammaS = new QLabel(tr("γ<sub>S</sub>:"));
    spingammaS = new QDoubleSpinBox();
    spingammaS->setMinimum(0.);
    spingammaS->setMaximum(10.);
    spingammaS->setDecimals(4);
    spingammaS->setValue(model->Parameters().gammaS);
    spingammaS->setToolTip(tr("Chemical non-equilibrium factor for strange quarks"));

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
		layParameters->addWidget(labelC, 2, 6, 1, 1, Qt::AlignRight);
		layParameters->addWidget(spinC, 2, 7);

    grParameters->setLayout(layParameters);

    QHBoxLayout *layQB = new QHBoxLayout();
    layQB->setAlignment(Qt::AlignLeft);

		checkFixMuQ = new QCheckBox(tr("Constrain μQ"));
		checkFixMuQ->setChecked(true);
    checkFixMuQ->setToolTip(tr("Constrain μQ to reproduce the needed Q/B ratio"));
		connect(checkFixMuQ, SIGNAL(clicked()), this, SLOT(modelChanged()));
		//connect(checkFixMuS, SIGNAL(clicked()), this, SLOT(modelChanged()));
		//connect(checkFixMuC, SIGNAL(clicked()), this, SLOT(modelChanged()));

    QLabel *labelQB = new QLabel(tr("Q/B ratio:"));
    spinQBRatio = new QDoubleSpinBox();
    spinQBRatio->setDecimals(3);
    spinQBRatio->setMinimum(0.);
    spinQBRatio->setValue(model->QoverB());

		checkFixMuS = new QCheckBox(tr("Constrain μS"));
		checkFixMuS->setChecked(true);
    checkFixMuS->setToolTip(tr("Constrain μS to obtain zero net strangeness"));
		checkFixMuC = new QCheckBox(tr("Constrain μC"));
		checkFixMuC->setChecked(true);
    checkFixMuC->setToolTip(tr("Constrain μS to obtain zero net charm"));


		QLabel *labelVolumeRSC = new QLabel(tr("R<sub>C</sub>:"));
		spinVolumeRSC = new QDoubleSpinBox();
		spinVolumeRSC->setDecimals(4);
		spinVolumeRSC->setMinimum(0.);
		spinVolumeRSC->setMaximum(25.);
		spinVolumeRSC->setValue(spinVolumeR->value());
		spinVolumeRSC->setEnabled(false);
    spinVolumeRSC->setToolTip(tr("Correlation radius: the (canonical) correlation volume is a sphere of this radius"));

		
    layQB->addWidget(labelQB);
    layQB->addWidget(spinQBRatio);
		layQB->addWidget(checkFixMuQ);
		layQB->addWidget(checkFixMuS);
		layQB->addWidget(checkFixMuC);
		layQB->addWidget(labelVolumeRSC);
		layQB->addWidget(spinVolumeRSC);


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
		checkFluctuations = new QCheckBox(tr("Fluctuations"));
		checkFluctuations->setChecked(false);
    checkFluctuations->setToolTip(tr("Compute fluctuation observables (viewable in \"Equation of state\" dialog)"));

		layFlags->addWidget(labelWidth);
		layFlags->addWidget(comboWidth);
    layFlags->addWidget(checkBratio);
		layFlags->addWidget(checkFluctuations);

    QGroupBox *grParallel = new QGroupBox(tr("Paralellization:"));

    QHBoxLayout *layParallel = new QHBoxLayout();
    layParallel->setAlignment(Qt::AlignLeft);

    checkOMP = new QCheckBox(tr("OpenMP"));
    checkOMP->setChecked(false);
    layParallel->addWidget(checkOMP);

    grParallel->setLayout(layParallel);

    QHBoxLayout *layButtons = new QHBoxLayout();
    layButtons->setAlignment(Qt::AlignLeft);

    buttonCalculate = new QPushButton(tr("Calculate"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));

    buttonCalculateFitted = new QPushButton(tr("Calculate from fit tab"));
    buttonCalculateFitted->setToolTip(tr("Use parameters found from a thermal fit in the other tab (here no rounding errors in input due to spin boxes)"));
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
    teDebug->setFontPointSize(10);

		editorLay->addLayout(layModelEnsemble);
    editorLay->addWidget(grStats);
    editorLay->addWidget(grEV);
    editorLay->addWidget(grParameters);
    editorLay->addLayout(layQB);
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
}

void ModelTab::loadEVFromFile()
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

ThermalModelConfig ModelTab::getConfig()
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

	ret.QuantumStatistics        = static_cast<int>(!radioBoltz->isChecked());
	ret.QuantumStatisticsType    = static_cast<int>(CBQuadratures->isChecked());
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

	ret.T         = spinTemperature->value() * 1.e-3;
	ret.muB       = spinmuB->value() * 1.e-3;
	ret.muQ       = spinmuQ->value() * 1.e-3;
	ret.muS       = spinmuS->value() * 1.e-3;
	ret.muC       = spinmuC->value() * 1.e-3;
	ret.gq        = spingammaq->value();
	ret.gS        = spingammaS->value();
	ret.gC        = 1.;
	ret.VolumeR   = spinVolumeR->value();
	ret.VolumeRSC = spinVolumeRSC->value();
	ret.B         = spinB->value();
	ret.Q         = spinQ->value();
	ret.S         = spinS->value();
	ret.C         = spinC->value();

	ret.QoverB       = spinQBRatio->value();
	ret.ConstrainMuQ = checkFixMuQ->isChecked();
	ret.ConstrainMuS = checkFixMuS->isChecked();
	ret.ConstrainMuC = checkFixMuC->isChecked();

	ret.FiniteWidth        = comboWidth->currentIndex();
	ret.RenormalizeBR      = checkBratio->isChecked();
	ret.ComputeFluctations = checkFluctuations->isChecked();

	return ret;
}

ThermalModelConfig ModelTab::getConfigFromFit(thermalfist::ThermalModelFit * fit, const ThermalModelConfig & configfit)
{
  ThermalModelConfig ret = configfit;

  /*ret.ModelType = static_cast<int>(fit->model()->InteractionModel());
  ThermalModelBase::ThermalModelEnsemble ens = fit->model()->Ensemble();
  if (ens != ThermalModelBase::GCE) {
    if (ens == ThermalModelBase::CE)
      ret.ModelType = ThermalModelConfig::CE;

    if (ens == ThermalModelBase::SCE) {
      if (fit->model()->InteractionModel() == ThermalModelBase::Ideal)
        ret.ModelType = ThermalModelConfig::SCE;
      if (fit->model()->InteractionModel() == ThermalModelBase::DiagonalEV)
        ret.ModelType = ThermalModelConfig::EVSCE;
      if (fit->model()->InteractionModel() == ThermalModelBase::QvdW)
        ret.ModelType = ThermalModelConfig::VDWSCE;
    }

    if (ens == ThermalModelBase::CCE)
      ret.ModelType = ThermalModelConfig::CCE;
  }

  ret.QuantumStatistics = static_cast<int>( fit->model()->QuantumStatistics() );
  ret.QuantumStatisticsType = static_cast<int>(fit->model()->TPS()->QStatsCalculationType());

  bool pionsOnly = true, mesonsOnly = true;
  for (int i = 0; i < fit->model()->TPS()->Particles().size(); ++i) {
    const ThermalParticle &part = fit->model()->TPS()->Particles()[i];
    if (part.Statistics() != 0) {
      if (part.PdgId() != 211 && part.PdgId() != 111 && part.PdgId() != -211) {
        pionsOnly = false;
      }
      if (part.BaryonCharge() != 0) {
        mesonsOnly = false;
      }
    }
  }



  if (radioUniform->isChecked())
    ret.Interaction = 0;
  else if (radioBaglike->isChecked())
    ret.Interaction = 1;
  else if (radioMesons->isChecked())
    ret.Interaction = 2;
  else
    ret.Interaction = 3;

  ret.EVRadius = spinRadius->value();
  ret.InteractionInput = strEVPath.toStdString();*/

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

  ret.ComputeFluctations = checkFluctuations->isChecked();

  return ret;
}

void ModelTab::updateControlsWithConfig(const ThermalModelConfig & config)
{
  if (config.ModelType == ThermalModelConfig::DiagonalEV
    || config.ModelType == ThermalModelConfig::EVSCE)
    radEVD->setChecked(true);
  else if (config.ModelType == ThermalModelConfig::CrosstermsEV)
    radEVCRS->setChecked(true);
  else if (config.ModelType == ThermalModelConfig::QvdW
    || config.ModelType == ThermalModelConfig::VDWSCE)
    radQVDW->setChecked(true);
  else
    radIdeal->setChecked(true);

  if (config.ModelType == ThermalModelConfig::CE)
    radCE->setChecked(true);
  else if (config.ModelType == ThermalModelConfig::SCE
    || config.ModelType == ThermalModelConfig::EVSCE
    || config.ModelType == ThermalModelConfig::VDWSCE)
    radSCE->setChecked(true);
  else
    radGCE->setChecked(true);

  if (config.QuantumStatistics)
    radioQuant->setChecked(true);
  else
    radioBoltz->setChecked(true);

  if (config.QuantumStatisticsInclude == 1) {
    CBBoseOnly->setChecked(true);
    CBPionsOnly->setChecked(false);
  }
  else if (config.QuantumStatisticsInclude == 2) {
    CBBoseOnly->setChecked(false);
    CBPionsOnly->setChecked(true);
  }

  CBQuadratures->setChecked(config.QuantumStatisticsType);

  spinRadius->setValue(config.EVRadius);

  if (config.Interaction == 1)
    radioBaglike->setChecked(true);
  else if (config.Interaction == 2)
    radioMesons->setChecked(true);
  else if (config.Interaction == 3)
    radioCustomEV->setChecked(true);
  else
    radioUniform->setChecked(true);

  spinTemperature->setValue(config.T * 1.e3);
  spinmuB->setValue(config.muB * 1.e3);
  spinmuQ->setValue(config.muQ * 1.e3);
  spinmuS->setValue(config.muS * 1.e3);
  spinmuC->setValue(config.muC * 1.e3);
  spingammaq->setValue(config.gq);
  spingammaS->setValue(config.gS);
  spinVolumeR->setValue(config.VolumeR);
  spinVolumeRSC->setValue(config.VolumeRSC);
  spinB->setValue(config.B);
  spinQ->setValue(config.Q);
  spinS->setValue(config.S);
  spinC->setValue(config.C);

  spinQBRatio->setValue(config.QoverB);
  checkFixMuQ->setChecked(config.ConstrainMuQ);
  checkFixMuS->setChecked(config.ConstrainMuS);
  checkFixMuC->setChecked(config.ConstrainMuC);
  comboWidth->setCurrentIndex(config.FiniteWidth);
  checkBratio->setChecked(config.RenormalizeBR);
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

	if (model != NULL) 
		modelnew->SetNormBratio(config.RenormalizeBR);

	myModel->setModel(modelnew);
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
	model->SetStrangenessCanonicalVolumeRadius(config.VolumeRSC);

	dbgstrm << "T\t= " << model->Parameters().T * 1.e3 << " MeV" << endl;
	dbgstrm << "muB\t= " << model->Parameters().muB * 1.e3 << " MeV" << endl;

	model->SetBaryonCharge(config.B);
	model->SetElectricCharge(config.Q);
	model->SetStrangeness(config.S);
	model->SetCharm(config.C);

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

	model->SetQoverB(config.QoverB);
	model->ConstrainMuQ(config.ConstrainMuQ);
	model->ConstrainMuS(config.ConstrainMuS);
	model->ConstrainMuC(config.ConstrainMuC);

	model->SetNormBratio(config.RenormalizeBR);

	model->SetStatistics(config.QuantumStatistics);
	//model->SetClusterExpansionOrder(10);

	if (config.QuantumStatisticsType)
		model->SetCalculationType(IdealGasFunctions::Quadratures);
	else
		model->SetCalculationType(IdealGasFunctions::ClusterExpansion);


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

	QElapsedTimer timervdw;
	timervdw.start();

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

	printf("Parameters time = %d ms\n", timervdw.elapsed());

	// If fluctuations are calculated within the CE one needs a twice larger range of quantum numbers
	if (config.ModelType == ThermalModelConfig::CE) {
		static_cast<ThermalModelCanonical*>(model)->CalculateQuantumNumbersRange(config.ComputeFluctations);
	}

	printf("Initialization time = %d ms\n", timerc.elapsed());

	timerc.restart();

	//model->FixParameters();
	model->ConstrainChemicalPotentials();

	model->CalculateDensities();

	printf("Densities time = %d ms\n", timerc.elapsed());

	timerc.restart();

	if (config.ComputeFluctations) {
		model->CalculateFluctuations();

		computeHigherOrderFluctuations();
	}

	printf("Fluctuations time = %d ms\n", timerc.elapsed());

	timerc.restart();

	if (config.ModelType != ThermalModelConfig::CE) {
		if (config.ModelType != ThermalModelConfig::SCE && config.ModelType != ThermalModelConfig::EVSCE && config.ModelType != ThermalModelConfig::VDWSCE)
			dbgstrm << "muS\t= " << model->Parameters().muS * 1.e3 << " MeV" << endl;

		dbgstrm << "muQ\t= " << model->Parameters().muQ * 1.e3 << " MeV" << endl;
		if (model->TPS()->hasCharmed() && config.ModelType != ThermalModelConfig::CCE)
			dbgstrm << "muC\t= " << model->Parameters().muC * 1.e3 << " MeV" << endl;
	}
	else {
		dbgstrm << "B\t= " << model->CalculateBaryonDensity()      * model->Volume() << endl;
		dbgstrm << "S\t= " << model->CalculateStrangenessDensity() * model->Volume() << endl;
		dbgstrm << "Q\t= " << model->CalculateChargeDensity()      * model->Volume() << endl;
		if (model->TPS()->hasCharmed())
			dbgstrm << "Q\t= " << model->CalculateCharmDensity()      * model->Volume() << endl;
	}
	dbgstrm << "gammaq\t= " << model->Parameters().gammaq << endl;
	dbgstrm << "gammaS\t= " << model->Parameters().gammaS << endl;
	if (model->TPS()->hasCharmed())
		dbgstrm << "gammaC\t= " << model->Parameters().gammaC << endl;
	dbgstrm << "V\t= " << model->Volume() << " fm^3" << endl;
	dbgstrm << endl;
	dbgstrm << "Particle density\t= " << model->CalculateHadronDensity() << " fm^-3" << endl;
	dbgstrm << "Net baryon density\t= " << model->CalculateBaryonDensity() << " fm^-3" << endl;
	dbgstrm << "Net baryon number\t= " << model->CalculateBaryonDensity() * model->Volume() << endl;
	dbgstrm << "Net electric charge\t= " << model->CalculateChargeDensity() * model->Volume() << endl;
	dbgstrm << "Net strangeness\t= " << model->CalculateStrangenessDensity() * model->Volume() << endl;
	if (model->TPS()->hasCharmed())
		dbgstrm << "Net charm\t\t= " << model->CalculateCharmDensity() * model->Volume() << endl;
	dbgstrm << "E/N\t\t= " << model->CalculateEnergyDensity() / model->CalculateHadronDensity() << endl;
	dbgstrm << "E/Nb\t\t= " << model->CalculateEnergyDensity() / model->CalculateBaryonDensity() << endl;
	dbgstrm << "S/Nb\t\t= " << model->CalculateEntropyDensity() / model->CalculateBaryonDensity() << endl;
	dbgstrm << "S/N\t\t= " << model->CalculateEntropyDensity() / model->CalculateHadronDensity() << endl;
	dbgstrm << "Q/B\t\t= " << model->CalculateChargeDensity() / model->CalculateBaryonDensity() << endl;
	if (model->TPS()->hasStrange())
		dbgstrm << "S/|S|\t\t= " << model->CalculateStrangenessDensity() / model->CalculateAbsoluteStrangenessDensity() << endl;
	if (model->TPS()->hasCharmed())
		dbgstrm << "C/|C|\t\t= " << model->CalculateCharmDensity() / model->CalculateAbsoluteCharmDensity() << endl;
	dbgstrm << endl;
	dbgstrm << endl;
	dbgstrm << "Calculation time = " << timer.elapsed() << " ms" << endl;
	dbgstrm << "----------------------------------------------------------" << endl;
	teDebug->append(dbgstr);
	dbgstr.clear();

	printf("Finalizing time = %d ms\n", timerc.elapsed());

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
	performCalculation(getConfig());

	if (!radCE->isChecked()) {
		if (!(radSCE->isChecked()))
			spinmuS->setValue(model->Parameters().muS*1.e3);
		spinmuQ->setValue(model->Parameters().muQ*1.e3);
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
                fprintf(f, "%20s%15d%8d%8d%8d%8d%8d%15E%15E\n",
                        model->TPS()->Particles()[i].Name().c_str(),
                        model->TPS()->Particles()[i].PdgId(),
                        model->TPS()->Particles()[i].IsStable(),
                        model->TPS()->Particles()[i].BaryonCharge(),
                        model->TPS()->Particles()[i].ElectricCharge(),
                        model->TPS()->Particles()[i].Strangeness(),
                        model->TPS()->Particles()[i].Charm(),
                        model->GetParticlePrimordialDensity(i) * model->Volume(),
                        model->GetParticleTotalDensity(i) * model->Volume());
            }
            fclose(f);
        }
    }
}

double mnucl = 0.938;

double ssqrt(double elab) {
    return sqrt(2*mnucl*(elab+mnucl));
}

double muBss(double ss) {
    return 1.308 / (1. + 0.273 * ss);
}


double Tss(double ss) {
    double tmpmuB = muBss(ss);
    return 0.166 - 0.139 * tmpmuB * tmpmuB - 0.053 * tmpmuB * tmpmuB * tmpmuB * tmpmuB;
}

double gammaSss(double ss) {
    double tmpmuB = muBss(ss);
    double tmpT = Tss(ss);
    return 1. - 0.396 * exp(-1.23 * tmpT / tmpmuB);
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
}

void ModelTab::modelChanged()
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

				if (radEVCRS->isChecked()) {
					radSCE->setEnabled(false);
					if (radSCE->isChecked())
						radGCE->setChecked(true);
				}
				else
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
				spinC->setEnabled(false);
        spinQBRatio->setEnabled(checkFixMuQ->isChecked());
        radioBoltz->setEnabled(true);
        radioQuant->setEnabled(true);

				checkFixMuQ->setEnabled(true);
				checkFixMuS->setEnabled(true);
				checkFixMuC->setEnabled(true);

				CBQuadratures->setEnabled(true);
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
        spinQBRatio->setEnabled(false);

				checkFixMuQ->setEnabled(false);
				checkFixMuS->setEnabled(false);
				checkFixMuC->setEnabled(false);

				CBQuadratures->setChecked(false);
				CBQuadratures->setEnabled(false);

				spinVolumeRSC->setEnabled(true);
    }


		if (radSCE->isChecked()) {
			checkFixMuS->setEnabled(false);
			spinmuS->setEnabled(false);
			spinVolumeRSC->setEnabled(true);
		}
		else if (!radCE->isChecked()) spinVolumeRSC->setEnabled(false);


		spinmuS->setVisible(model->TPS()->hasStrange());
		spinmuC->setVisible(model->TPS()->hasCharmed());
		labelmuS->setVisible(model->TPS()->hasStrange());
		labelmuC->setVisible(model->TPS()->hasCharmed());
		checkFixMuS->setVisible(model->TPS()->hasStrange());
		checkFixMuC->setVisible(model->TPS()->hasCharmed());

		labelQ->setVisible(model->TPS()->hasCharged());
		spinQ->setVisible(model->TPS()->hasCharged());
		labelS->setVisible(model->TPS()->hasStrange());
		spinS->setVisible(model->TPS()->hasStrange());
		labelC->setVisible(model->TPS()->hasCharmed());
		spinC->setVisible(model->TPS()->hasCharmed());

		if (radioBoltz->isChecked()) {
			CBBoseOnly->setEnabled(false);
			CBPionsOnly->setEnabled(false);
		}
		else {
			CBBoseOnly->setEnabled(true);
			CBPionsOnly->setEnabled(true);
		}
}
