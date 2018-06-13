#include "energydependencetab.h"

#include <algorithm>

#include <QLayout>
#include <QGroupBox>
#include <QLabel>

#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGEV/ThermalModelEVCrossterms.h"
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonicalCharm.h"
#include "HRGBase/xMath.h"

EnergyDependenceTab::EnergyDependenceTab(QWidget *parent, ThermalModelBase *modelop) :
    QWidget(parent)
{
    TPS = modelop->TPS();
    model = new ThermalModelIdeal(modelop->TPS());
    *model = *modelop;

    int index = 0;
    QString tname;

    tname = "T (MeV)";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "muB (MeV)";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "muS (MeV)";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "muQ (MeV)";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "gammaS";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "n (fm^-3)";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "rhoB (fm^-3)";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "rhoQ (e/fm^3)";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "e (MeV/fm^3)";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    //tname = "s (fm^-3)";
    tname = "s/T^3";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    //tname = "P (MeV/fm^3)";
    tname = "P/T^4";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "(e-3P)/T^4";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "eta/s";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "w";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    tname = "K+/pi+";
    paramnames.push_back(tname);
    parammap[tname] = index;
    index++;

    QHBoxLayout *mainLayout = new QHBoxLayout();

    QVBoxLayout *layLeft = new QVBoxLayout();

    QHBoxLayout *selLay = new QHBoxLayout();
    selLay->setAlignment(Qt::AlignLeft);
    QLabel *labelSel = new QLabel(tr("Quantity:"));
    comboQuantity = new QComboBox();
    for(int i=0;i<index;++i) comboQuantity->addItem(paramnames[i]);
    comboQuantity->setCurrentIndex(0);
    connect(comboQuantity, SIGNAL(currentIndexChanged(int)), this, SLOT(replot()));
    selLay->addWidget(labelSel);
    selLay->addWidget(comboQuantity);

    plotDependence = new QCustomPlot();
    plotDependence->xAxis->setLabel("Ecm (GeV)");
    plotDependence->yAxis->setLabel(comboQuantity->currentText());

    plotDependence->addGraph();
    plotDependence->graph(0)->setName(comboQuantity->currentText());

    layLeft->addLayout(selLay);
    layLeft->addWidget(plotDependence);

    QVBoxLayout *layRight = new QVBoxLayout();
    layRight->setAlignment(Qt::AlignTop);


    QGroupBox *grModels = new QGroupBox(tr("Model:"));

    //QHBoxLayout *layQuant = new QHBoxLayout();
    QHBoxLayout *layModels = new QHBoxLayout();
    layModels->setAlignment(Qt::AlignLeft);
    radioIdeal = new QRadioButton(tr("Ideal gas"));
    radioEVMF = new QRadioButton(tr("MF VdW gas"));
    radioIdealCanonStrangeness = new QRadioButton(tr("IG canon. strangeness"));
    radioIdealCanonCharm = new QRadioButton(tr("IG canon. charm"));

    layModels->addWidget(radioIdeal);
    layModels->addWidget(radioEVMF);
    layModels->addWidget(radioIdealCanonStrangeness);
    layModels->addWidget(radioIdealCanonCharm);

    grModels->setLayout(layModels);

    radioIdeal->setChecked(true);

    QGroupBox *grStats = new QGroupBox(tr("Statistics:"));

    QHBoxLayout *layStats = new QHBoxLayout();
    layStats->setAlignment(Qt::AlignLeft);
    radioBoltz = new QRadioButton(tr("Boltzmann"));
    radioQuant = new QRadioButton(tr("Quantum"));

    layStats->addWidget(radioBoltz);
    layStats->addWidget(radioQuant);

    grStats->setLayout(layStats);

    radioQuant->setChecked(true);

    QHBoxLayout *layRadius = new QHBoxLayout();
    layRadius->setAlignment(Qt::AlignLeft);
    QLabel *labelRadius = new QLabel(tr("Hadron radius (fm):"));
    spinRadius = new QDoubleSpinBox();
    spinRadius->setMinimum(0.);
    spinRadius->setMaximum(100.);
    spinRadius->setValue(0.5);
    layRadius->addWidget(labelRadius);
    layRadius->addWidget(spinRadius);

    QHBoxLayout *layQB = new QHBoxLayout();
    layQB->setAlignment(Qt::AlignLeft);
    QLabel *labelQB = new QLabel(tr("Q/B ratio:"));
    spinQBRatio = new QDoubleSpinBox();
    spinQBRatio->setDecimals(3);
    spinQBRatio->setMinimum(0.);
    spinQBRatio->setValue(model->QoverB());

    layQB->addWidget(labelQB);
    layQB->addWidget(spinQBRatio);

    QHBoxLayout *layFlags = new QHBoxLayout();
    layFlags->setAlignment(Qt::AlignLeft);
    checkFiniteWidth = new QCheckBox(tr("Finite resonance width"));
    checkFiniteWidth->setChecked(model->UseWidth());
    checkBratio = new QCheckBox(tr("Renormalize branching ratios"));
    checkBratio->setChecked(false);

    layFlags->addWidget(checkFiniteWidth);
    layFlags->addWidget(checkBratio);

    QGroupBox *grParallel = new QGroupBox(tr("Paralellization:"));

    //QHBoxLayout *layQuant = new QHBoxLayout();
    QHBoxLayout *layParallel = new QHBoxLayout();
    layParallel->setAlignment(Qt::AlignLeft);

    checkOMP = new QCheckBox(tr("OpenMP"));
    checkOMP->setChecked(false);
    layParallel->addWidget(checkOMP);

    grParallel->setLayout(layParallel);

    QGroupBox *grRange = new QGroupBox(tr("Energy range:"));

    QGridLayout *layBounds = new QGridLayout();
    layBounds->setAlignment(Qt::AlignLeft);
    QLabel *labelEnMin = new QLabel(tr("Ecm<sub>min</sub> (GeV):"));
    spinEnMin = new QDoubleSpinBox();
    spinEnMin->setDecimals(3);
    spinEnMin->setMinimum(2*0.938);
    spinEnMin->setMaximum(1e5);
    spinEnMin->setValue(2*0.938);

    QLabel *labelEnMax = new QLabel(tr("Ecm<sub>max</sub> (GeV):"));
    spinEnMax = new QDoubleSpinBox();
    spinEnMax->setDecimals(3);
    spinEnMax->setMinimum(2*0.938);
    spinEnMax->setMaximum(1e5);
    spinEnMax->setValue(20);

    QLabel *labelEnIters = new QLabel(tr("Iters:"));
    spinEnIters = new QSpinBox();
    spinEnIters->setMinimum(1);
    spinEnIters->setMaximum(100000);
    spinEnIters->setValue(20);

    layBounds->addWidget(labelEnMin, 0, 0);
    layBounds->addWidget(spinEnMin, 0, 1);
    layBounds->addWidget(labelEnMax, 1, 0);
    layBounds->addWidget(spinEnMax, 1, 1);
    layBounds->addWidget(labelEnIters, 2, 0);
    layBounds->addWidget(spinEnIters, 2, 1);

    grRange->setLayout(layBounds);

    buttonCalculate = new QPushButton(tr("Calculate"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));

    layRight->addWidget(grModels);
    layRight->addWidget(grStats);
    layRight->addLayout(layRadius);
    layRight->addLayout(layQB);
    layRight->addLayout(layFlags);
    layRight->addWidget(grParallel);
    layRight->addWidget(grRange);
    layRight->addWidget(buttonCalculate, 0, Qt::AlignLeft);

    mainLayout->addLayout(layLeft, 1);
    mainLayout->addLayout(layRight);

    setLayout(mainLayout);

    fRunning = false;

    calcTimer = new QTimer(this);
    connect(calcTimer, SIGNAL(timeout()), this, SLOT(replot()));
}

EnergyDependenceTab::~EnergyDependenceTab() {
    if (model!=NULL) delete model;
}

void EnergyDependenceTab::calculate() {
    if (!fRunning) {
        ThermalModelBase *modelnew;

        if (radioIdeal->isChecked()) modelnew = new ThermalModelIdeal(model->TPS());
        else if (radioEVMF->isChecked()) modelnew = new ThermalModelEVDiagonal(model->TPS());
        else if (radioIdealCanonStrangeness->isChecked()) modelnew = new ThermalModelCanonicalStrangeness(model->TPS());
        else modelnew = new ThermalModelCanonicalCharm(model->TPS());

        if (model!=NULL) delete model;
        model = modelnew;

        params.resize(0);
        varvalues.resize(0);

        double de = (spinEnMax->value()-spinEnMin->value()) / spinEnIters->value();

        model->SetStatistics(!radioBoltz->isChecked());
        model->SetUseWidth(checkFiniteWidth->isChecked());
        model->SetQoverB(spinQBRatio->value());
        model->SetOMP(checkOMP->isChecked());
        if (checkBratio->isChecked()) model->TPS()->NormalizeBranchingRatios();
        else model->TPS()->RestoreBranchingRatios();
        model->SetNormBratio(checkBratio->isChecked());

        fCurrentSize = 0;

        params.resize(spinEnIters->value());
        varvalues.resize(spinEnIters->value());

        fStop = 0;

        EnergyDependenceWorker *wrk = new EnergyDependenceWorker(model, spinEnMin->value(), de, spinEnIters->value(), spinRadius->value(),
                                                                 &params, &varvalues, &fCurrentSize, &fStop, this);
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

void EnergyDependenceTab::finalize() {
    calcTimer->stop();
    buttonCalculate->setText(tr("Calculate"));
    fRunning = false;
    replot();
}

void EnergyDependenceTab::replot() {
    int index = comboQuantity->currentIndex();
    if (index>=0 && index<paramnames.size()) {
        plotDependence->graph(0)->clearData();
        plotDependence->graph(0)->setName(paramnames[index]);
        plotDependence->yAxis->setLabel(paramnames[index]);
        double tmin = 0., tmax = 0.;
        //int tsz = varvalues.size();
        for(int i=0;i<fCurrentSize;++i) {
            tmin = std::min(tmin, params[i][index]);
            tmax = std::max(tmax, params[i][index]);
            plotDependence->graph(0)->addData(varvalues[i], params[i][index]);
        }
        plotDependence->xAxis->setRange(spinEnMin->value(), spinEnMax->value());
        plotDependence->yAxis->setRange(tmin, 1.1*tmax);
        plotDependence->replot();
    }
}

std::vector<double> EnergyDependenceTab::getParams(ThermalModelBase *modelop) {
    std::vector<double> ret(0);
    ret.push_back(modelop->Parameters().T);
    ret.push_back(modelop->Parameters().muB);
    ret.push_back(modelop->Parameters().muS);
    ret.push_back(modelop->Parameters().muQ);
    ret.push_back(modelop->Parameters().gammaS);
    ret.push_back(modelop->CalculateHadronDensity());
    ret.push_back(modelop->CalculateBaryonDensity());
    ret.push_back(modelop->CalculateChargeDensity());
    ret.push_back(modelop->CalculateEnergyDensity());
    ret.push_back(modelop->CalculateEntropyDensity()/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/xMath::GeVtoifm()/xMath::GeVtoifm()/xMath::GeVtoifm());//ret.push_back(modelop->CalculateEntropyDensity());
    //ret.push_back(modelop->CalculatePressure());
    ret.push_back(modelop->CalculatePressure()/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/xMath::GeVtoifm()/xMath::GeVtoifm()/xMath::GeVtoifm());
    ret.push_back((modelop->CalculateEnergyDensity()-3.*modelop->CalculatePressure())/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/xMath::GeVtoifm()/xMath::GeVtoifm()/xMath::GeVtoifm());
    ret.push_back(modelop->CalculateShearViscosity() / modelop->CalculateEntropyDensity());
    ret.push_back(modelop->CalculateHadronScaledVariance());
		if (modelop->TPS()->PdgToId(321) != -1 && modelop->TPS()->PdgToId(211) != -1)
			ret.push_back(modelop->GetDensity(321, 1) / modelop->GetDensity(211, 1));
		else
			ret.push_back(0.); 
		return ret;
}

std::vector<double> EnergyDependenceWorker::getParams(ThermalModelBase *modelop) {
    std::vector<double> ret(0);
    ret.push_back(modelop->Parameters().T);
    ret.push_back(modelop->Parameters().muB);
    ret.push_back(modelop->Parameters().muS);
    ret.push_back(modelop->Parameters().muQ);
    ret.push_back(modelop->Parameters().gammaS);
    ret.push_back(modelop->CalculateHadronDensity());
    ret.push_back(modelop->CalculateBaryonDensity());
    ret.push_back(modelop->CalculateChargeDensity());
    ret.push_back(modelop->CalculateEnergyDensity());
    ret.push_back(modelop->CalculateEntropyDensity()/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/xMath::GeVtoifm()/xMath::GeVtoifm()/xMath::GeVtoifm());//ret.push_back(modelop->CalculateEntropyDensity());
    //ret.push_back(modelop->CalculatePressure());
    ret.push_back(modelop->CalculatePressure()/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/xMath::GeVtoifm()/xMath::GeVtoifm()/xMath::GeVtoifm());
    ret.push_back((modelop->CalculateEnergyDensity()-3.*modelop->CalculatePressure())/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/modelop->Parameters().T/xMath::GeVtoifm()/xMath::GeVtoifm()/xMath::GeVtoifm());
    ret.push_back(modelop->CalculateShearViscosity() / modelop->CalculateEntropyDensity());
    ret.push_back(modelop->CalculateHadronScaledVariance());
		if (modelop->TPS()->PdgToId(321) != -1 && modelop->TPS()->PdgToId(211) != -1)
			ret.push_back(modelop->GetDensity(321, 1) / modelop->GetDensity(211, 1));
		else
			ret.push_back(0.); 
		return ret;
}
