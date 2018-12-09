#include "contourplottab.h"

#include <QLayout>
#include <QGroupBox>
#include <QLabel>
#include <QDebug>

#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonicalCharm.h"
#include "HRGBase/xMath.h"

using namespace thermalfist;

ContourPlotTab::ContourPlotTab(QWidget *parent, ThermalModelBase *modelop) :
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

    plot = new QCustomPlot();
    plot->xAxis->setLabel(tr("muB (MeV)"));
    plot->yAxis->setLabel(tr("T (MeV)"));

    colormap = new QCPColorMap(plot->xAxis, plot->yAxis);

    plot->addPlottable(colormap);

    colorScale = new QCPColorScale(plot);
    plot->plotLayout()->addElement(0, 1, colorScale);

    //plot->addGraph();
    //plot->graph(0)->setName(comboQuantity->currentText());

    layLeft->addLayout(selLay);
    layLeft->addWidget(plot);

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

    QGroupBox *grRange = new QGroupBox(tr("Range:"));

    QGridLayout *layBounds = new QGridLayout();
    layBounds->setAlignment(Qt::AlignLeft);
    QLabel *labelTMin = new QLabel(tr("T<sub>min</sub> (MeV):"));
    spinTMin = new QDoubleSpinBox();
    spinTMin->setDecimals(3);
    spinTMin->setMinimum(0.);
    spinTMin->setMaximum(300.);
    spinTMin->setValue(50.);

    QLabel *labelTMax = new QLabel(tr("T<sub>max</sub> (MeV):"));
    spinTMax = new QDoubleSpinBox();
    spinTMax->setDecimals(3);
    spinTMax->setMinimum(0.);
    spinTMax->setMaximum(300.);
    spinTMax->setValue(180.);

    QLabel *labelTIters = new QLabel(tr("Iters:"));
    spinTIters = new QSpinBox();
    spinTIters->setMinimum(1);
    spinTIters->setMaximum(100000);
    spinTIters->setValue(20);

    QLabel *labelmuMin = new QLabel(tr("μB<sub>min</sub> (MeV):"));
    spinmuMin = new QDoubleSpinBox();
    spinmuMin->setDecimals(3);
    spinmuMin->setMinimum(0.);
    spinmuMin->setMaximum(1000.);
    spinmuMin->setValue(50.);

    QLabel *labelmuMax = new QLabel(tr("μB<sub>max</sub> (MeV):"));
    spinmuMax = new QDoubleSpinBox();
    spinmuMax->setDecimals(3);
    spinmuMax->setMinimum(0.);
    spinmuMax->setMaximum(1000.);
    spinmuMax->setValue(800.);

    QLabel *labelmuIters = new QLabel(tr("Iters:"));
    spinmuIters = new QSpinBox();
    spinmuIters->setMinimum(1);
    spinmuIters->setMaximum(100000);
    spinmuIters->setValue(20);

    layBounds->addWidget(labelTMin, 0, 0);
    layBounds->addWidget(spinTMin, 0, 1);
    layBounds->addWidget(labelTMax, 0, 2);
    layBounds->addWidget(spinTMax, 0, 3);
    layBounds->addWidget(labelTIters, 0, 4);
    layBounds->addWidget(spinTIters, 0, 5);
    layBounds->addWidget(labelmuMin, 1, 0);
    layBounds->addWidget(spinmuMin, 1, 1);
    layBounds->addWidget(labelmuMax, 1, 2);
    layBounds->addWidget(spinmuMax, 1, 3);
    layBounds->addWidget(labelmuIters, 1, 4);
    layBounds->addWidget(spinmuIters, 1, 5);

    grRange->setLayout(layBounds);

    buttonCalculate = new QPushButton(tr("Calculate"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));

    progBar = new QProgressBar();
    progBar->setVisible(false);

    layRight->addWidget(grModels);
    layRight->addWidget(grStats);
    layRight->addLayout(layRadius);
    layRight->addLayout(layQB);
    layRight->addLayout(layFlags);
    layRight->addWidget(grParallel);
    layRight->addWidget(grRange);
    layRight->addWidget(buttonCalculate, 0, Qt::AlignLeft);
    layRight->addWidget(progBar);

    mainLayout->addLayout(layLeft, 1);
    mainLayout->addLayout(layRight);

    setLayout(mainLayout);

    fRunning = false;

    calcTimer = new QTimer(this);
    connect(calcTimer, SIGNAL(timeout()), this, SLOT(updateProgress()));
}

ContourPlotTab::~ContourPlotTab() {
    if (model!=NULL) delete model;
}

void ContourPlotTab::calculate() {
    if (!fRunning) {
        ThermalModelBase *modelnew;

        if (radioIdeal->isChecked()) modelnew = new ThermalModelIdeal(model->TPS());
        else if (radioEVMF->isChecked()) modelnew = new ThermalModelEVDiagonal(model->TPS());
        else if (radioIdealCanonStrangeness->isChecked()) modelnew = new ThermalModelCanonicalStrangeness(model->TPS());
        else modelnew = new ThermalModelCanonicalCharm(model->TPS());

        if (model!=NULL) delete model;
        model = modelnew;

        params.resize(spinmuIters->value()*spinTIters->value());
        Tvalues.resize(0);
        muBvalues.resize(0);

        model->SetStatistics(!radioBoltz->isChecked());
        model->SetUseWidth(checkFiniteWidth->isChecked());
        model->SetQoverB(spinQBRatio->value());
        model->SetOMP(checkOMP->isChecked());
        if (checkBratio->isChecked()) model->TPS()->NormalizeBranchingRatios();
        else model->TPS()->RestoreBranchingRatios();
        model->SetNormBratio(checkBratio->isChecked());


        double dT = (spinTMax->value() - spinTMin->value()) / spinTIters->value();
        double dMu = (spinmuMax->value() - spinmuMin->value()) / spinmuIters->value();
        double Tst = spinTMin->value();
        double muBst = spinmuMin->value();

        for(int i=0;i<spinmuIters->value();++i)
            for(int j=0;j<spinTIters->value();++j) {
                double tmpmuB = muBst + (0.5 + i)*dMu;
                double tmpT = Tst + (0.5 + j)*dT;
                //double gammaS = 1.;

                Tvalues.push_back(tmpT);
                muBvalues.push_back(tmpmuB);
            }

        fCurrentSize = 0;
        fTotalSize = spinmuIters->value() * spinTIters->value();
        fStop = 0;

        ContourPlotWorker *wrk = new ContourPlotWorker(model, &Tvalues, &muBvalues, &params,
                                                       spinRadius->value(),
                                                       &fCurrentSize, &fStop, this);

        connect(wrk, SIGNAL(calculated()), this, SLOT(finalize()));
        connect(wrk, SIGNAL(finished()), wrk, SLOT(deleteLater()));

        wrk->start();

        buttonCalculate->setText(tr("Stop"));
        fRunning = true;
        progBar->setVisible(true);

        calcTimer->start(10);

        //replot();
    }
    else {
        fStop = 1;
    }
}

void ContourPlotTab::finalize() {
    calcTimer->stop();
    buttonCalculate->setText(tr("Calculate"));
    fRunning = false;
    progBar->setVisible(false);
    replot();
}

void ContourPlotTab::updateProgress() {
    progBar->setRange(0, fTotalSize);
    progBar->setValue(fCurrentSize);
}

void ContourPlotTab::replot() {
    int index = comboQuantity->currentIndex();
    if (index>=0 && index<paramnames.size()) {

        plot->xAxis->setRange(spinmuMin->value(),spinmuMax->value());
        plot->yAxis->setRange(spinTMin->value(),spinTMax->value());
        QCPColorMapData *data = new QCPColorMapData(spinmuIters->value(),
                                                    spinTIters->value(),
                                                    plot->xAxis->range(),
                                                    plot->yAxis->range());
        for(int i=0;i<fCurrentSize;++i) {
           data->setData(muBvalues[i], Tvalues[i], params[i][index]);
           //qDebug() << Tvalues[i] << " " << muBvalues[i] << " " << params[i][index] ;
        }

        //colormap->setGradient(QCPColorGradient::gpPolar);
        //colormap->setGradient(colormap->gradient().inverted());
        QCPColorGradient grad;
        QMap<double, QColor> gradmap;
        gradmap.insert(0., Qt::yellow);
        gradmap.insert(1., Qt::blue);
        //gradmap.insert(0., Qt::black);
        //gradmap.insert(1., Qt::red);
        grad.setColorStops(gradmap);
        colormap->setGradient(grad);

        colormap->setData(data, true);
        colormap->rescaleDataRange(true);

        colorScale->setGradient(colormap->gradient());
        colorScale->setDataRange(colormap->dataRange());
        //colormap->setColorScale(colorScale);
        //colorScale->setDataRange(colormap->dataRange());

        delete data;

        plot->rescaleAxes(true);
        plot->replot();
    }
}

std::vector<double> ContourPlotTab::getParams(ThermalModelBase *modelop) {
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
		if (modelop->TPS()->PdgToId(321) != -1 && modelop->TPS()->PdgToId(211) != -1)
			ret.push_back(modelop->GetDensity(321, 1) / modelop->GetDensity(211, 1));
		else
			ret.push_back(0.);
		//ret.push_back(modelop->densitiestotal[modelop->TPS()->PDGtoID[321]] / modelop->densitiestotal[modelop->TPS()->PDGtoID[211]]);
    return ret;
}

std::vector<double> ContourPlotWorker::getParams(ThermalModelBase *modelop) {
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
		if (modelop->TPS()->PdgToId(321) != -1 && modelop->TPS()->PdgToId(211) != -1)
			ret.push_back(modelop->GetDensity(321, 1) / modelop->GetDensity(211, 1));
		else
			ret.push_back(0.);
		//ret.push_back(modelop->densitiestotal[modelop->TPS()->PDGtoID[321]] / modelop->densitiestotal[modelop->TPS()->PDGtoID[211]]);
    return ret;
}
