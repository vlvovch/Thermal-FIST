#ifndef CONTOURPLOTTAB_H
#define CONTOURPLOTTAB_H

#include <QWidget>

#include <QWidget>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QComboBox>
#include <QRadioButton>
#include <QProgressBar>
#include <QThread>
#include <QTimer>
#include <map>

#include "HRGBase/ThermalModelBase.h"
#include "HRGFit/ThermalModelFit.h"
#include "qcustomplot.h"

// This class is currently not in use


class QuantitiesModel;
//class ThermalParticleSystem;
//class ThermalModelBase;

/*struct CalcParameter {
    std::string name;
    double value;
};*/

class ContourPlotWorker : public QThread
{
    Q_OBJECT

    thermalfist::ThermalModelBase *model;
    double rad;
    int *currentSize;
    int *stop;

    std::vector< std::vector<double> > *params;
    std::vector<double> *Tvalues;
    std::vector<double> *muBvalues;

    void run() Q_DECL_OVERRIDE {
        for(int i=0;i<Tvalues->size() && !(*stop);++i) {
            double tmpmuB = muBvalues->operator [](i);
            double tmpT = Tvalues->operator [](i);
            double gammaS = 1.;

						model->SetParameters(thermalfist::ThermalModelParameters(tmpT*1.e-3, tmpmuB*1.e-3,
							tmpmuB*1.e-3 / 5., -tmpmuB*1.e-3 / 50.,
							gammaS, 4000.));

            //model->SetParameters(tmpT*1.e-3, tmpmuB*1.e-3,
            //                     tmpmuB*1.e-3/5., -tmpmuB*1.e-3/50.,
            //                     gammaS, 4000., rad);

            model->FixParameters();
            model->CalculateDensities();
            params->operator[](i) = getParams(model);
            (*currentSize)++;
        }
        emit calculated();
    }

    std::vector<double> getParams(thermalfist::ThermalModelBase *modelop);

public:
    ContourPlotWorker(thermalfist::ThermalModelBase *mod = NULL,
           std::vector<double> *Tvalueso = NULL,
           std::vector<double> *muBvalueso = NULL,
           std::vector< std::vector<double> > *paramso = NULL,
           double rad = 0.5,
           int *currentSizeo = NULL,
           int *stopo = NULL,
           QObject * parent = 0) :
        QThread(parent), rad(rad) {
            model = mod;
            params = paramso;
            Tvalues = Tvalueso;
            muBvalues = muBvalueso;
            currentSize = currentSizeo;
            stop = stopo;
    }
signals:
    void calculated();
};

class ContourPlotTab : public QWidget
{
    Q_OBJECT

    QComboBox *comboQuantity;
    QCustomPlot *plot;
    QCPColorMap *colormap;
    QCPColorScale *colorScale;

    std::vector< std::vector<double> > params;
    std::vector<double> Tvalues;
    std::vector<double> muBvalues;
    int fCurrentSize;
    int fTotalSize;
    bool fRunning;
    int fStop;

    std::map<QString, int> parammap;
    std::vector<QString> paramnames;

    QRadioButton *radioIdeal, *radioEVMF, *radioIdealCanonStrangeness, *radioIdealCanonCharm;

    QRadioButton *radioBoltz, *radioQuant;

    QCheckBox *checkOnlyStable;

    QDoubleSpinBox *spinQBRatio;
    QDoubleSpinBox *spinRadius;

    QCheckBox *checkFiniteWidth;
    QCheckBox *checkBratio;

    QCheckBox *checkOMP;

    QDoubleSpinBox *spinTMin, *spinTMax;
    QSpinBox *spinTIters;
    QDoubleSpinBox *spinmuMin, *spinmuMax;
    QSpinBox *spinmuIters;

    QPushButton *buttonCalculate;

    QProgressBar *progBar;

    thermalfist::ThermalParticleSystem *TPS;
    thermalfist::ThermalModelBase *model;

    QTimer *calcTimer;
public:
    ContourPlotTab(QWidget *parent = 0, thermalfist::ThermalModelBase *model=NULL);
    ~ContourPlotTab();
signals:

public slots:
    void calculate();
    void replot();
    void updateProgress();
    void finalize();
private:
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

    std::vector<double> getParams(thermalfist::ThermalModelBase *modelop);
};

#endif // CONTOURPLOTTAB_H
