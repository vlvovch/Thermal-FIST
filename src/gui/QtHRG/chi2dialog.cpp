#include "chi2dialog.h"

#include <QLayout>
#include <QLabel>
#include <QDebug>

#include "qcustomplot.h"

using namespace thermalfist;

chi2Dialog::chi2Dialog(QWidget *parent, ThermalModelFit *mod) :
    QDialog(parent), modelFit(mod)
{
    QHBoxLayout *layout = new QHBoxLayout();

    plot = new QCustomPlot;
    plot->xAxis->setLabel("muB (MeV)");
    plot->yAxis->setLabel("T (MeV)");
    plot->xAxis->setRange(50, 700);
    plot->yAxis->setRange(50, 180);

    colormap = new QCPColorMap(plot->xAxis, plot->yAxis);

    plot->addPlottable(colormap);

    colorScale = new QCPColorScale(plot);
    plot->plotLayout()->addElement(0, 1, colorScale);
    //colorScale->setLabel("Some Label Text");

    QVBoxLayout *layInterface = new QVBoxLayout();
    layInterface->setAlignment(Qt::AlignTop);

    QGridLayout *layParams = new QGridLayout();

    QLabel *labelTmin = new QLabel(tr("T<sub>min</sub>:"));
    spinTmin = new QDoubleSpinBox();
    spinTmin->setMinimum(0.);
    spinTmin->setMaximum(800.);
    spinTmin->setValue(50.);
    QLabel *labelTmax = new QLabel(tr("T<sub>max</sub>:"));
    spinTmax = new QDoubleSpinBox();
    spinTmax->setMinimum(0.);
    spinTmax->setMaximum(800.);
    spinTmax->setValue(180.);
    QLabel *labelTiter = new QLabel(tr("T<sub>iters</sub>:"));
    spinTiter = new QSpinBox();
    spinTiter->setMinimum(1);
    spinTiter->setMaximum(100000);
    spinTiter->setValue(10);

    QLabel *labelmuBmin = new QLabel(tr("μB<sub>min</sub>:"));
    spinmuBmin = new QDoubleSpinBox();
    spinmuBmin->setMinimum(0.);
    spinmuBmin->setMaximum(800.);
    spinmuBmin->setValue(50.);
    QLabel *labelmuBmax = new QLabel(tr("μB<sub>max</sub>:"));
    spinmuBmax = new QDoubleSpinBox();
    spinmuBmax->setMinimum(0.);
    spinmuBmax->setMaximum(1000.);
    spinmuBmax->setValue(800.);
    QLabel *labelmuBiter = new QLabel(tr("μB<sub>iters</sub>:"));
    spinmuBiter = new QSpinBox();
    spinmuBiter->setMinimum(1);
    spinmuBiter->setMaximum(100000);
    spinmuBiter->setValue(10);

    layParams->addWidget(labelTmin, 0, 0);
    layParams->addWidget(spinTmin, 0, 1);
    layParams->addWidget(labelTmax, 0, 2);
    layParams->addWidget(spinTmax, 0, 3);
    layParams->addWidget(labelTiter, 0, 4);
    layParams->addWidget(spinTiter, 0, 5);
    layParams->addWidget(labelmuBmin, 1, 0);
    layParams->addWidget(spinmuBmin, 1, 1);
    layParams->addWidget(labelmuBmax, 1, 2);
    layParams->addWidget(spinmuBmax, 1, 3);
    layParams->addWidget(labelmuBiter, 1, 4);
    layParams->addWidget(spinmuBiter, 1, 5);

    QHBoxLayout *laychi2Range = new QHBoxLayout();
    laychi2Range->setAlignment(Qt::AlignLeft);
    QLabel *labelChi2min = new QLabel(tr("χ<sup>2</sup><sub>min</sub>:"));
    spinchi2min = new QDoubleSpinBox();
    spinchi2min->setMinimum(0.);
    spinchi2min->setMaximum(100.);
    spinchi2min->setValue(0.);
    QLabel *labelChi2max = new QLabel(tr("χ<sup>2</sup><sub>max</sub>:"));
    spinchi2max = new QDoubleSpinBox();
    spinchi2max->setMinimum(0.);
    spinchi2max->setMaximum(100.);
    spinchi2max->setValue(5.);
    laychi2Range->addWidget(labelChi2min);
    laychi2Range->addWidget(spinchi2min);
    laychi2Range->addWidget(labelChi2max);
    laychi2Range->addWidget(spinchi2max);


    QHBoxLayout *layButtons = new QHBoxLayout();
    layButtons->setAlignment(Qt::AlignLeft);
    buttonCalculate = new QPushButton(tr("Calculate contour"));
    connect(buttonCalculate, SIGNAL(clicked()), this, SLOT(calculate()));
    buttonReplot = new QPushButton(tr("Replot"));
    connect(buttonReplot, SIGNAL(clicked()), this, SLOT(replot()));
    layButtons->addWidget(buttonCalculate);
    layButtons->addWidget(buttonReplot);

    progBar = new QProgressBar();
    progBar->setVisible(false);

    layInterface->addLayout(layParams);
    layInterface->addLayout(laychi2Range);
    layInterface->addLayout(layButtons);
    layInterface->addWidget(progBar);
    //layInterface->addWidget(buttonCalculate, 0, Qt::AlignLeft);


    layout->addWidget(plot, 1);
    layout->addLayout(layInterface);

    plot->replot();

    setLayout(layout);

    fRunning = false;

    calcTimer = new QTimer(this);
    connect(calcTimer, SIGNAL(timeout()), this, SLOT(updateProgress()));
    connect(calcTimer, SIGNAL(timeout()), this, SLOT(replotpro()));
}

void chi2Dialog::calculate() {
    if (!fRunning) {
        plot->xAxis->setRange(spinmuBmin->value(),spinmuBmax->value());
        plot->yAxis->setRange(spinTmin->value(),spinTmax->value());

        double dT = (spinTmax->value() - spinTmin->value()) / spinTiter->value();
        double dMu = (spinmuBmax->value() - spinmuBmin->value()) / spinmuBiter->value();
        double Tst = spinTmin->value();
        double muBst = spinmuBmin->value();
        params.resize(spinmuBiter->value()*spinTiter->value());
        Tvalues.resize(0);
        muBvalues.resize(0);
        for(int i=0;i<spinmuBiter->value();++i)
            for(int j=0;j<spinTiter->value();++j) {
                double tmpmuB = muBst + (0.5 + i)*dMu;
                double tmpT = Tst + (0.5 + j)*dT;
                Tvalues.push_back(tmpT);
                muBvalues.push_back(tmpmuB);
                //data->setData(tmpmuB, tmpT, modelFit->chi2Ndf(tmpT*1.e-3, tmpmuB*1.e-3));
                //qDebug() << tmpT << " " << tmpmuB << " " << modelFit->chi2Ndf(tmpT*1.e-3, tmpmuB*1.e-3);
            }

        fCurrentSize = 0;
        fPreviousSize = 0;
        fTotalSize = spinmuBiter->value() * spinTiter->value();
        fStop = 0;

        QCPColorMapData *data = new QCPColorMapData(spinmuBiter->value(),
                                                        spinTiter->value(),
                                                        plot->xAxis->range(),
                                                        plot->yAxis->range());
        colormap->setData(data, true);
        delete data;

        chi2Worker *wrk = new chi2Worker(modelFit, &Tvalues, &muBvalues, &params,
                                                       &fCurrentSize, &fStop, this);

        connect(wrk, SIGNAL(calculated()), this, SLOT(finalize()));
        connect(wrk, SIGNAL(finished()), wrk, SLOT(deleteLater()));

        wrk->start();

        buttonCalculate->setText(tr("Stop"));
        buttonReplot->setEnabled(false);
        fRunning = true;
        progBar->setVisible(true);

        calcTimer->start(10);
    }
    else {
        fStop = 1;
    }
}

void chi2Dialog::replot() {
    colorScale->setDataRange(QCPRange(spinchi2min->value(), spinchi2max->value()));
    colormap->setColorScale(colorScale);
    plot->replot();
}

void chi2Dialog::updateProgress() {
    progBar->setRange(0, fTotalSize);
    progBar->setValue(fCurrentSize);
}

void chi2Dialog::replotpro() {
    /*QCPColorMapData *data = new QCPColorMapData(spinmuBiter->value(),
                                                spinTiter->value(),
                                                plot->xAxis->range(),
                                                plot->yAxis->range());*/

    int tsz = fCurrentSize;
    for(int i=fPreviousSize;i<tsz;++i) {
        colormap->data()->setData(muBvalues[i], Tvalues[i], params[i]);
    }

    fPreviousSize = tsz;

    colormap->setGradient(QCPColorGradient::gpPolar);
    colormap->setGradient(colormap->gradient().inverted());
    //colormap->setData(data, true);
    //colormap->rescaleDataRange(true);

    colorScale->setGradient(colormap->gradient());
    colorScale->setDataRange(QCPRange(spinchi2min->value(), spinchi2max->value()));
    colormap->setColorScale(colorScale);
    //colorScale->setDataRange(colormap->dataRange());

    //delete data;

    plot->rescaleAxes(true);
    plot->replot();
}

void chi2Dialog::finalize() {
    calcTimer->stop();
    buttonCalculate->setText(tr("Calculate"));
    buttonReplot->setEnabled(true);
    fRunning = false;
    progBar->setVisible(false);
    replotpro();
}

