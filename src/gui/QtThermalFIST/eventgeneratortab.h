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
#include "HRGEventGenerator/HepMCEventWriter.h"

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

    thermalfist::HepMCEventWriter hepmcout;

    void run() Q_DECL_OVERRIDE;

public:
  EventGeneratorWorker(
    thermalfist::EventGeneratorBase* gen = NULL,
    ParticlesSpectra* spec = NULL,
    QMutex* mut = NULL,
    int totalEvents = 0,
    int* evproc = NULL,
    int* stopo = NULL,
    double* nEp = NULL,
    bool pDecays = false,
    std::string fileout = "",
    QObject* parent = 0);

signals:
    void calculated();
};


class BinningDialog : public QDialog
{
  Q_OBJECT
public:
  struct BinningConfig {
    int bins1D;
    int bins2D_x, bins2D_y;
    BinningConfig(int b1d = 500, int b2dx = 40, int b2dy = 40) :
      bins1D(b1d), bins2D_x(b2dx), bins2D_y(b2dy) { }
  };

  BinningConfig* m_config;
  QSpinBox* spinBins;
  QSpinBox* spinBinsX;
  QSpinBox* spinBinsY;

public:
  explicit  BinningDialog(BinningConfig* config = 0, QWidget* parent = 0);

signals:

public slots:
  void OK();
  void Discard() { QDialog::reject(); };
};


class ParticlesAnalyzeDialog : public QDialog
{
  Q_OBJECT
public:
  struct ParticlesAnalyzeConfig {
    int type;
    std::set<long long> pdgCodes;
    ParticlesAnalyzeConfig() :
      type(1) {
      pdgCodes.insert(221);
      pdgCodes.insert(113);
      pdgCodes.insert(313);
    }
  };

  ParticlesAnalyzeConfig* m_config;
  QRadioButton* RBAll;
  QRadioButton* RBStable;
  QRadioButton* RBStablePlus;
  QLineEdit* lePdgs;

public:
  explicit  ParticlesAnalyzeDialog(ParticlesAnalyzeConfig* config = 0, QWidget* parent = 0);

signals:

public slots:
  void OK();
  void Discard() { QDialog::reject(); };
};

class EventGeneratorTab : public QWidget
{
    Q_OBJECT

    QTableView *tableSpectra;
    ParticlesAnalyzeDialog::ParticlesAnalyzeConfig partsConfig;
    //QCheckBox* checkStableOnly;
    //int particlesAnalyzeInclude;
    //std::set<int> pdgCodes;

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
    QLabel *labelBeta, *labelBetat;
    QDoubleSpinBox *spinTemperature, *spinmuB, *spingammaq, *spingammaS, *spingammaC, *spinmuS, *spinmuQ, *spinmuC, *spinVolumeR;
    QDoubleSpinBox *spinVolumeRSC;
    QLabel *labelVolumeVal, *labelVolumeSCVal;
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

    QString cpath;

    std::vector<double> fXv, fYv, fZv;
    std::vector<double> fZvErr;

    BinningDialog::BinningConfig binConfig;

    double prevBeta, prevRmax;
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


    void contextMenuRequestPlot1D(QPoint pos);
    void saveAsPdf1D() { return saveAs1D(0); }
    void saveAsPng1D() { return saveAs1D(1); }
    void saveAsAscii1D() { return saveAs1D(2); }
    // saveAs type: 0 - pdf, 1 - png, 2 - ascii data
    void saveAs1D(int type);

    void contextMenuRequestPlot2D(QPoint pos);
    void saveAsPdf2D() { return saveAs2D(0); }
    void saveAsPng2D() { return saveAs2D(1); }
    void saveAsAscii2D() { return saveAs2D(2); }
    // saveAs type: 0 - pdf, 1 - png, 2 - ascii data
    void saveAs2D(int type);

    void changeParticles();
    void changeBinning();
};

#endif // FITTOEXPERIMENTTAB_H
