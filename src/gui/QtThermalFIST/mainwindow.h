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
#include <QSet>

#include <string>
#include <vector>
#include <set>

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

    thermalfist::ThermalParticleSystem *TPS;
    thermalfist::ThermalModelBase *model;

    QString cpath  = "";
    QString clists = "";

    // S-matrix / phase-shift state (configured via the Model configuration dialog)
    bool m_phaseShiftsEnabled = false;             ///< master on/off for phase shifts
    QStringList m_phaseShiftConfs;                  ///< phase-shift config file(s), applied in order
    QSet<QString> m_phaseShiftDisabledChannels;     ///< channels switched off individually
    std::vector<std::string> m_lastListPaths;      ///< base list files of the current selection
    std::vector<std::string> m_lastDecayPaths;     ///< base decay files of the current selection
    std::set<long long> m_basePdgCodes;            ///< PDG codes of the base list captured before
                                                   ///< the last phase-shift apply; lets stripPhaseShifts
                                                   ///< tell created clusters from overridden resonances

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
    /// Used for list (re)loads - discards in-memory edits by design.
    void rebuildCurrentList();
    /// Re-apply the phase-shift selection while PRESERVING the current in-memory
    /// list (including list-editor edits): strips the previous phase-shift
    /// additions in place instead of reloading from files, then re-applies.
    void reapplyPhaseShifts();
    /// Remove the last phase-shift apply's additions from TPS, recovering the base
    /// list: drop created clusters/resonances, clear overriding densities. Editor
    /// edits (no phase-shift density) are untouched. Uses m_basePdgCodes.
    ///
    /// \param fullRollback when true (only valid right after m_basePdgCodes was
    ///   snapshotted, i.e. the catch of a failed apply), ALSO removes any species
    ///   not in the snapshot that lack a phase-shift density - orphan module
    ///   particles left by an apply that threw before attaching densities. Must NOT
    ///   be used for a normal reapply, where post-snapshot editor additions are not
    ///   in m_basePdgCodes yet would be wrongly removed.
    void stripPhaseShifts(bool fullRollback = false);
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
    /// Open the phase-shift configuration dialog (from the Model configuration);
    /// on accept, store the new selection and rebuild the list.
    void openPhaseShiftsDialog();
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
