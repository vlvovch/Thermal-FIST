#ifndef CHI2DIALOG_H
#define CHI2DIALOG_H

#include <QDialog>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QPushButton>
#include <QProgressBar>
#include <QThread>
#include <QTimer>
#include "HRGFit/ThermalModelFit.h"

//class ThermalModelFit;
class QCPColorMap;
class QCustomPlot;
class QCPColorScale;

class chi2Worker : public QThread
{
    Q_OBJECT

    ThermalModelFit *modelFit;
    int *currentSize;
    int *stop;

    std::vector< double > *params;
    std::vector<double> *Tvalues;
    std::vector<double> *muBvalues;

    void run() Q_DECL_OVERRIDE {
        for(int i=0;i<Tvalues->size() && !(*stop);++i) {
            double tmpmuB = muBvalues->operator [](i);
            double tmpT = Tvalues->operator [](i);
            params->operator[](i) = modelFit->chi2Ndf(tmpT*1.e-3, tmpmuB*1.e-3);
            (*currentSize)++;
        }
        emit calculated();
    }

public:
    chi2Worker(ThermalModelFit *mod = NULL,
           std::vector<double> *Tvalueso = NULL,
           std::vector<double> *muBvalueso = NULL,
           std::vector<double> *paramso = NULL,
           int *currentSizeo = NULL,
           int *stopo = NULL,
           QObject * parent = 0) :
        QThread(parent) {
            modelFit = mod;
            params = paramso;
            Tvalues = Tvalueso;
            muBvalues = muBvalueso;
            currentSize = currentSizeo;
            stop = stopo;
    }
signals:
    void calculated();
};

class chi2Dialog : public QDialog
{
    Q_OBJECT

    //ThermalParticle *fParticle;
    //ThermalParticleSystem *fTPS;
    std::vector<double> params;
    std::vector<double> Tvalues;
    std::vector<double> muBvalues;
    int fPreviousSize;
    int fCurrentSize;
    int fTotalSize;
    bool fRunning;
    int fStop;

    ThermalModelFit *modelFit;

    QCustomPlot *plot;
    QCPColorMap *colormap;
    QCPColorScale *colorScale;

    QDoubleSpinBox *spinTmin, *spinTmax;
    QSpinBox *spinTiter;

    QDoubleSpinBox *spinmuBmin, *spinmuBmax;
    QSpinBox *spinmuBiter;

    QDoubleSpinBox *spinchi2min, *spinchi2max;

    QPushButton *buttonCalculate;
    QPushButton *buttonReplot;

    QProgressBar *progBar;
    QTimer *calcTimer;

    QString GetParameters();
    QString GetResults();
public:
    explicit  chi2Dialog(QWidget *parent = 0, ThermalModelFit *mod = NULL);

signals:
public slots:
    void calculate();
    void replot();
    void replotpro();
    void finalize();
    void updateProgress();
};

#endif // CHI2DIALOG_H

