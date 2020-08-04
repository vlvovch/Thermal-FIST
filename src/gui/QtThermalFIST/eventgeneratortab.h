/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef EVENTGENERATORTAB_H
#define EVENTGENERATORTAB_H

#include <QWidget>
#include <QLineEdit>
#include <QPushButton>
#include <QTableView>
#include <QTableWidget>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QComboBox>
#include <QLabel>
#include <QRadioButton>
#include <QTextEdit>
#include <QProgressBar>
#include <QHBoxLayout>
#include <QThread>
#include <QMutex>
#include <QElapsedTimer>

#include "HRGBase/ThermalModelBase.h"
#include "HRGEventGenerator/EventGeneratorBase.h"
#include "BaseStructures.h"
#include "configwidgets.h"

#include <cmath>
#include <fstream>

class SpectraModel;
class ParticlesSpectra;
class QCustomPlot;
class QCPColorMap;
class QCPColorScale;
class QStackedWidget;
class QCPErrorBars;

class EventGeneratorWorker : public QThread
{
    Q_OBJECT

    thermalfist::EventGeneratorBase *generator;
    ParticlesSpectra *spectra;
    QMutex *mutex;
    int events;
    double wsum, w2sum;
    int *eventsProcessed;
    int *stop;
    double *nE;
    bool performDecays;


    std::ofstream fout;

    void run() Q_DECL_OVERRIDE;

public:
    EventGeneratorWorker(
           thermalfist::EventGeneratorBase *gen=NULL,
           ParticlesSpectra *spec=NULL,
           QMutex *mut=NULL,
           int totalEvents = 0,
           int *evproc = NULL,
           int *stopo = NULL,
           double *nEp = NULL,
           bool pDecays = false,
           std::string fileout = "",
           QObject * parent = 0) :
        QThread(parent), generator(gen), spectra(spec), mutex(mut),
            events(totalEvents), eventsProcessed(evproc), stop(stopo), nE(nEp), performDecays(pDecays)
    {
        wsum = w2sum = 0.;
        fout.clear();
        if (fileout != "")
          fout.open(fileout.c_str());
    }
signals:
    void calculated();
};

class EventGeneratorTab : public QWidget
{
    Q_OBJECT

    QTableView *tableSpectra;

    QComboBox *comboDistr;
    QStackedWidget *plot;
    QCustomPlot *plotDistr;
    QCustomPlot *plot2D;
    QCPColorMap *colormap;
    QCPColorScale *colorScale;
    QCPErrorBars *errorBars;
    std::map<QString, int> parammap;
    std::vector<QString> paramnames, paramnamesx;


    QLabel *labelmuS, *labelmuC, *labelgammaS, *labelgammaC;
    QLabel *labelBetat;
    QDoubleSpinBox *spinTemperature, *spinmuB, *spingammaq, *spingammaS, *spingammaC, *spinmuS, *spinmuQ, *spinmuC, *spinVolumeR;
    QDoubleSpinBox *spinVolumeRSC;
    QLabel *labelB, *labelQ, *labelS, *labelC;
    QSpinBox *spinB, *spinS, *spinQ, *spinC;

    QRadioButton *radioSR, *radioSSH, *radioCracow;

    QDoubleSpinBox *spinTkin;
    QDoubleSpinBox *spinBeta;
    QDoubleSpinBox *spinBetat, *spinEtaMax, *spinn;

    QSpinBox *spinEvents;

    QCheckBox *checkDecays;

    QPushButton *buttonCalculate;

    thermalfist::ThermalParticleSystem *TPS;
    thermalfist::ThermalModelBase *model;
    thermalfist::EventGeneratorBase *generator;
    ParticlesSpectra *spectra;

    QTextEdit *teDebug;

    SpectraModel *myModel;

    QElapsedTimer timer;

    QHBoxLayout *layProgress;
    QProgressBar *progBar;
    QLabel *labelProgress;

    QMutex mutex;

    int getCurrentRow();

    int fCurrentSize;
    int fTotalSize;
    double nE;
    bool fRunning;
    int fStop;
    QTimer *calcTimer;

    QCheckBox *checkAcceptance;
    QLineEdit *leAcceptancePath;

    QCheckBox *checkFile;
    QLineEdit *leFilePath;
    QPushButton *buttonChooseFile;

    ModelConfigWidget *configWidget;
public:
    EventGeneratorTab(QWidget *parent = 0, thermalfist::ThermalModelBase *model=NULL);
    ~EventGeneratorTab();
    ThermalModelConfig getConfig();
signals:

public slots:
    void changedRow();
    void calculate();
    void replot();
    void replot(const QVector<double> &x1, const QVector<double> &y1, const QVector<double> &y1err,
                const QVector<double> &x2, const QVector<double> &y2, int index, double rightlimit = 2.);
    void replot2D(const QVector<double> &xv, const QVector<double> &yv, const QVector<double> &zv, int index, double rightlimit = 2.);
    void quantityDoubleClick(const QModelIndex &);
    void setModel(thermalfist::ThermalModelBase *model);
    void updateProgress();
    void finalize();
    void changePlot();
    void modelChanged();
    void resetTPS();
    void loadAcceptance();
    void chooseOutputFile();
    void changeVolumeRSC(double);
    void changeTkin(double);

    void generateEvents(const ThermalModelConfig & config);
};

#endif // FITTOEXPERIMENTTAB_H
