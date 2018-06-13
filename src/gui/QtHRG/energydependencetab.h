#ifndef ENERGYDEPENDENCETAB_H
#define ENERGYDEPENDENCETAB_H

#include <QWidget>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QComboBox>
#include <QRadioButton>
#include <QThread>
#include <QTimer>
#include <map>

#include "HRGFit/ThermalModelFit.h"
#include "qcustomplot.h"

// This class is currently not in use


class QuantitiesModel;
class ThermalParticleSystem;
class ThermalModelBase;

/*struct CalcParameter {
    std::string name;
    double value;
};*/

class EnergyDependenceWorker : public QThread
{
    Q_OBJECT

    ThermalModelBase *model;
    double emin, de;
    double rad;
    int iterse;
    int *currentSize;
    int *stop;

    std::vector< std::vector<double> > *params;
    std::vector<double> *varvalues;

    void run() Q_DECL_OVERRIDE {
        for(int i=0;i<iterse && !(*stop);++i) {
            double ss = emin + de*(0.5+i);
            double TT = Tss(ss);
            double muB = muBss(ss);
            double gammaS = gammaSss(ss);

            //model->SetStatistics(!radioBoltz->isChecked());
            //model->SetUseWidth(checkFiniteWidth->isChecked());
            //model->SetQBgoal(spinQBRatio->value());
            //model->SetOMP(checkOMP->isChecked());
            //if (checkBratio->isChecked()) model->TPS->NormalizeBranchingRatios();
            //else model->TPS->RestoreBranchingRatios();
            //model->SetNormBratio(checkBratio->isChecked());

						model->SetParameters(ThermalModelParameters(TT, muB,
							muB / 5., -muB / 50.,
							gammaS, 4000.));

            //model->SetParameters(TT, muB,
            //                     muB/5., -muB/50.,
            //                     gammaS, 4000., rad);

            model->FixParameters();
            model->CalculateDensities();
            //varvalues->push_back(ss);
            //params->push_back(getParams(model));
            varvalues->operator [](i) = ss;
            params->operator [](i) = getParams(model);
            (*currentSize)++;
        }
        emit calculated();
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

    std::vector<double> getParams(ThermalModelBase *modelop);

public:
    EnergyDependenceWorker(ThermalModelBase *mod = NULL, double emin=2., double de = 1.,
           int iterse=10, double rad=0.5,
           std::vector< std::vector<double> > *paramso = NULL,
           std::vector<double> *varvalueso = NULL,
           int *currentSizeo = NULL,
           int *stopo = NULL,
           QObject * parent = 0) :
        QThread(parent), emin(emin), de(de), iterse(iterse), rad(rad) {
            model = mod;
            params = paramso;
            varvalues = varvalueso;
            currentSize = currentSizeo;
            stop = stopo;
    }
signals:
    void calculated();
};

class EnergyDependenceTab : public QWidget
{
    Q_OBJECT

    QComboBox *comboQuantity;
    QCustomPlot *plotDependence;
    std::vector< std::vector<double> > params;
    std::vector<double> varvalues;
    int fCurrentSize;
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

    QDoubleSpinBox *spinEnMin, *spinEnMax;
    QSpinBox *spinEnIters;

    QPushButton *buttonCalculate;

    ThermalParticleSystem *TPS;
    ThermalModelBase *model;

    QTimer *calcTimer;

public:
    EnergyDependenceTab(QWidget *parent = 0, ThermalModelBase *model=NULL);
    ~EnergyDependenceTab();
signals:

public slots:
    void calculate();
    void replot();
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

    std::vector<double> getParams(ThermalModelBase *modelop);
};


#endif // ENERGYDEPENDENCETAB_H
