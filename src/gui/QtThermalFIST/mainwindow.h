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
#include <QComboBox>
#include <QTableView>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QRadioButton>
#include <QTextEdit>
#include <QTabWidget>
#include <QStringList>

#include <string>
#include <vector>

#include "trajectoriestab.h"
#include "modeltab.h"
#include "eventgeneratortab.h"
#include "equationofstatetab.h"
#include "listeditortab.h"

#include "HRGBase/ThermalModelBase.h"
#include "fittoexperimenttab.h"

#include "CosmicEoStab.h"

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

    CosmicEoSTab *tabCosmicEoS;

    TrajectoriesTab *tabTrajectories;

    QLineEdit *leList;
    QComboBox *comboPDGEdition;
    QComboBox *comboListVariant;
    QPushButton *buttonLoad;
    QPushButton *buttonLoadDecays;
    QCheckBox   *chkPhaseShifts;
    QPushButton *buttonPhaseShifts;

    thermalfist::ThermalParticleSystem *TPS;
    thermalfist::ThermalModelBase *model;

    QString cpath  = "";
    QString clists = "";

    // S-matrix / phase-shift state
    QStringList m_phaseShiftConfs;                  ///< phase-shift config file(s), applied in order
    std::vector<std::string> m_lastListPaths;      ///< base list files of the current selection
    std::vector<std::string> m_lastDecayPaths;     ///< base decay files of the current selection

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
protected:
#ifndef QT_NO_CONTEXTMENU
//  void contextMenuEvent(QContextMenuEvent *event) override;
#endif // QT_NO_CONTEXTMENU
private:
    void createMenus();
    void applyParticleList(const std::vector<std::string>& listPaths,
                           const std::vector<std::string>& decayPaths,
                           const QString& displayName);
    static QString shortListDisplayName(const QString& fullPath);
    /// (Re)build TPS from the stored base list/decay files, (re)apply the
    /// phase-shift config if enabled, refresh the display and reset all tabs.
    void rebuildCurrentList();
    /// Add the phase-shift config channels to the current TPS (if enabled).
    void applyPhaseShiftsIfEnabled();
    /// Update the particle-list line edit (base name + optional PS tag + count).
    void refreshListDisplay();
    /// Reset every tab's view of the TPS.
    void resetAllTabs();
private slots:
    void loadList();
    void loadDecays();
    void switchParticleList();
    void onPhaseShiftToggled();
    void loadPhaseShiftConf();
    void updateListVariants();
    void tabChanged(int newIndex);
    void about();
    void documentation();
    void quickstartguide();
    void increaseFontSize();
    void decreaseFontSize();
#ifdef Q_OS_WASM
    void toggleFullscreen();
#endif
};

#endif // MAINWINDOW_H
