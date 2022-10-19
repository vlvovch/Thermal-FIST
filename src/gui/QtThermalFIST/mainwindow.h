/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLineEdit>
#include <QPushButton>
#include <QTableView>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QRadioButton>
#include <QTextEdit>
#include <QTabWidget>

#include "modeltab.h"
#include "eventgeneratortab.h"
#include "equationofstatetab.h"
#include "listeditortab.h"

#include "HRGBase/ThermalModelBase.h"
#include "fittoexperimenttab.h"


class MainWindow : public QMainWindow
{
    Q_OBJECT

    QTabWidget *tabWidget;
    int currentTab;

    ModelTab *tab1;
  
    FitToExperimentTab *tab2;

    EventGeneratorTab *tab5;
    EquationOfStateTab *tabEoS;
    ListEditorTab *tabEditor;

    QLineEdit *leList;
    QPushButton *buttonLoad;
    QPushButton *buttonLoadDecays;

    thermalfist::ThermalParticleSystem *TPS;
    thermalfist::ThermalModelBase *model;

    QString cpath  = "";
    QString clists = "";

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
protected:
#ifndef QT_NO_CONTEXTMENU
//  void contextMenuEvent(QContextMenuEvent *event) override;
#endif // QT_NO_CONTEXTMENU
private:
    void createMenus();
private slots:
    void loadList();
    void loadDecays();
    void tabChanged(int newIndex);
    void about();
    void documentation();
    void quickstartguide();
    void increaseFontSize();
    void decreaseFontSize();
};

#endif // MAINWINDOW_H
