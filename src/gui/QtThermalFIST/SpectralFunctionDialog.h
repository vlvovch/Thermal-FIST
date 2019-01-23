/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef SPECTRALFUNCTIONDIALOG_H
#define SPECTRALFUNCTIONDIALOG_H

#include <QDialog>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QSpinBox>
#include <QPushButton>
#include <QProgressBar>
#include <QThread>
#include <QComboBox>
#include <QTimer>
#include <map>

#include "HRGBase/ThermalParticle.h"
#include "BaseStructures.h"

//class ThermalModelFit;
class QCPColorMap;
class QCustomPlot;
class QCPColorScale;

class SpectralFunctionDialog : public QDialog
{
    Q_OBJECT

    //std::vector<double> params;
    //std::vector<double> Avalues;
    //int fPreviousSize;
    //int fCurrentSize;
    //int fTotalSize;

    std::vector< std::vector<double> > vecParams;
    std::vector< std::vector<double> > vecAvalues;
    std::vector< double > vecAleft, vecAright;
    //std::vector< int > vecPreviousSize;
    std::vector< int > vecCurrentSize;
    std::vector< int > vecTotalSize;

    bool fRunning;
    int fStop;

    thermalfist::ThermalParticle *particle;
    
    double T;
    double Mu;
    int WidthScheme; // 0 - BWTwoGamma, 1 - eBW


    QLabel *labelAmin, *labelAmax, *labelAiter;

    QComboBox *comboView;
    QComboBox *comboChannel;
    QCustomPlot *plot;
    
    //QCPColorMap *colormap;
    //QCPColorScale *colorScale;

    QComboBox *comboScheme;

    QDoubleSpinBox *spinT;
    QDoubleSpinBox *spinMu;

    double xleft, xright;
    double norm, normTH;
    std::vector<double> BWm;
    std::vector<double> BWval;
    std::vector<double> BWTHval;
    //std::vector<double> eBWm;
    //std::vector<double> eBWval;

public:
    explicit  SpectralFunctionDialog(QWidget *parent, thermalfist::ThermalParticle *part, double Temp, double Chem = 0., int scheme = 0);

signals:
public slots:
    void calculate();
    void replot();
    void replotSpectralFunctions();
    void replotBranchingRatios();
    void replotPartialWidths();
};

#endif // SPECTRALFUNCTIONDIALOG_H
