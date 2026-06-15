/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "mainwindow.h"

#include <QLayout>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QFileDialog>
#include <QLabel>
#include <QApplication>
#include <QMessageBox>
#include <QElapsedTimer>
#include <QDebug>
#include <QDesktopServices>
#include <QUrl>

#include "ThermalFISTConfig.h"
#include "WasmFileIO.h"
#ifdef Q_OS_WASM
#include <emscripten/html5.h>
#endif
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGCanonical/ThermalModelCanonicalStrangeness.h"
#include "HRGCanonical/ThermalModelCanonicalCharm.h"
#include "HRGPhaseShifts/PhaseShiftChannel.h"

#include "aboutdialog.h"

using namespace thermalfist;

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
{
#ifdef Q_OS_WASM
  // WASM: Extract default particle list from embedded resources to sandbox FS
  QString sandboxDir = WasmFileIO::getSandboxTempDir() + "/default";
  QDir().mkpath(sandboxDir);

  // Extract list.dat from resources
  QFile listRes(":/data/list.dat");
  QString listPath = sandboxDir + "/list.dat";
  if (listRes.open(QIODevice::ReadOnly)) {
    QFile listOut(listPath);
    if (listOut.open(QIODevice::WriteOnly)) {
      listOut.write(listRes.readAll());
      listOut.close();
    }
    listRes.close();
  }

  // Extract decays.dat from resources
  QFile decaysRes(":/data/decays.dat");
  QString decaysPath = sandboxDir + "/decays.dat";
  if (decaysRes.open(QIODevice::ReadOnly)) {
    QFile decaysOut(decaysPath);
    if (decaysOut.open(QIODevice::WriteOnly)) {
      decaysOut.write(decaysRes.readAll());
      decaysOut.close();
    }
    decaysRes.close();
  }

  // Extract default experimental data from resources
  QFile expDataRes(":/data/ALICE-PbPb2.76TeV-0-10-all.dat");
  QString expDataPath = sandboxDir + "/ALICE-PbPb2.76TeV-0-10-all.dat";
  if (expDataRes.open(QIODevice::ReadOnly)) {
    QFile expDataOut(expDataPath);
    if (expDataOut.open(QIODevice::WriteOnly)) {
      expDataOut.write(expDataRes.readAll());
      expDataOut.close();
    }
    expDataRes.close();
  }

  // Extract lattice QCD data from resources
  QString lqcdDir = sandboxDir + "/lqcd";
  QDir().mkpath(lqcdDir);

  QStringList lqcdFiles = {
    "WB-EoS.dat",
    "WB-chi2-1112.4416.dat",
    "WB-chi11-1910.14592.dat",
    "WB-chiB-1805.04445.dat",
    "HotQCD-EoS.dat",
    "HotQCD-chi2-1203.0784.dat",
    "HotQCD-chi2-chi8-2001.08530.dat"
  };

  for (const QString& lqcdFile : lqcdFiles) {
    QFile lqcdRes(":/data/lqcd/" + lqcdFile);
    QString lqcdPath = lqcdDir + "/" + lqcdFile;
    if (lqcdRes.open(QIODevice::ReadOnly)) {
      QFile lqcdOut(lqcdPath);
      if (lqcdOut.open(QIODevice::WriteOnly)) {
        lqcdOut.write(lqcdRes.readAll());
        lqcdOut.close();
      }
      lqcdRes.close();
    }
  }

  cpath = listPath;
  QString listpath = listPath;
  TPS = new ThermalParticleSystem(listpath.toStdString(), decaysPath.toStdString());
#else
  //cpath = QString(ThermalFIST_INPUT_FOLDER) + "/list/PDG2025/list-withnuclei.dat";
  cpath = QString(ThermalFIST_DEFAULT_LIST_FILE);

  QString listpath = cpath;
  TPS = new ThermalParticleSystem(listpath.toStdString());
#endif

  //TPS->SetSortMode(ThermalParticleSystem::SortByBaryonAndMassAndPDG);
  model = new ThermalModelIdeal(TPS);

  QWidget *centW = new QWidget;
  setCentralWidget(centW);

  //QVBoxLayout *dirTabLay1 = new QVBoxLayout();
  QHBoxLayout *dataLay = new QHBoxLayout();
  QLabel *labelData = new QLabel(tr("Particle list:"));
  dataLay->setAlignment(Qt::AlignLeft);

  // Combo box for PDG edition
  comboPDGEdition = new QComboBox();
  comboPDGEdition->addItem(tr("PDG2025"), QString("PDG2025"));
  comboPDGEdition->addItem(tr("PDG2020"), QString("PDG2020"));
  comboPDGEdition->addItem(tr("PDG2014"), QString("PDG2014"));
  comboPDGEdition->setCurrentIndex(0);
  comboPDGEdition->setPlaceholderText(tr("Custom"));
  comboPDGEdition->setMinimumWidth(100);
  connect(comboPDGEdition, QOverload<int>::of(&QComboBox::activated), this, [this](int) {
    updateListVariants();
    switchParticleList();
  });

  // Combo box for list variant
  comboListVariant = new QComboBox();
  comboListVariant->setPlaceholderText(tr("Custom"));
  comboListVariant->setMinimumWidth(180);
  updateListVariants();  // Populate variant combo for default PDG2025 (no switch during init)
  connect(comboListVariant, SIGNAL(activated(int)), this, SLOT(switchParticleList()));

  leList = new QLineEdit("");
  leList->setReadOnly(true);
  if (TPS->Particles().size() > 0)
    leList->setText(shortListDisplayName(listpath) + " + decays.dat [" + QString::number(TPS->Particles().size()) + " particles]");

  buttonLoad = new QPushButton(tr("Load list..."));
  connect(buttonLoad, SIGNAL(clicked()), this, SLOT(loadList()));

  buttonLoadDecays = new QPushButton(tr("Load decays..."));
  connect(buttonLoadDecays, SIGNAL(clicked()), this, SLOT(loadDecays()));

  // S-matrix / phase-shift controls. Default config ships with the package.
  m_phaseShiftConf = QString(ThermalFIST_INPUT_FOLDER) + "/list/phaseshifts/pipi.conf";
  chkPhaseShifts = new QCheckBox(tr("Phase shifts"));
  chkPhaseShifts->setToolTip(tr("Add S-matrix / phase-shift channels from a config file\n"
                                "as effective degrees of freedom (Beth-Uhlenbeck)."));
  chkPhaseShifts->setChecked(false);
  connect(chkPhaseShifts, SIGNAL(toggled(bool)), this, SLOT(onPhaseShiftToggled()));

  buttonPhaseShifts = new QPushButton(tr("Phase shifts..."));
  buttonPhaseShifts->setToolTip(tr("Choose the phase-shift config file."));
  connect(buttonPhaseShifts, SIGNAL(clicked()), this, SLOT(loadPhaseShiftConf()));

  dataLay->addWidget(labelData);
  dataLay->addWidget(comboPDGEdition);
  dataLay->addWidget(comboListVariant);
  dataLay->addWidget(leList, 1);
  dataLay->addWidget(buttonLoad);
  dataLay->addWidget(buttonLoadDecays);
  dataLay->addWidget(chkPhaseShifts);
  dataLay->addWidget(buttonPhaseShifts);

  // Remember the initially loaded base list so phase shifts can be (re)applied.
  m_lastListPaths = std::vector<std::string>(1, listpath.toStdString());
  m_lastDecayPaths = std::vector<std::string>();

  tab1 = new ModelTab(NULL, model);
  //tab1->resetTPS();

  tab2 = new FitToExperimentTab(NULL, model);
  tab1->setFitTab(tab2);

  tabEoS = new EquationOfStateTab(NULL, model);

  tabCosmicEoS = new CosmicEoSTab(NULL, model);

  tab5 = new EventGeneratorTab(NULL, model);
  tabEditor = new ListEditorTab(NULL, model);
  tabEditor->setListPath(cpath);

  tabTrajectories = new TrajectoriesTab(NULL, {tabEoS, tabCosmicEoS}, {"HRG model", "Cosmic trajectory"});

  tabWidget = new QTabWidget();
  tabWidget->addTab(tab1, QString(tr("Thermal model")));

  tabWidget->addTab(tab2, QString(tr("Thermal fits")));

  //tabWidget->addTab(tabEoS, QString(tr("Equation of state")));
  tabWidget->addTab(tabTrajectories, QString(tr("Equation of state")));

  tabWidget->addTab(tab5, QString(tr("Event generator")));

  //tabWidget->addTab(tabCosmicEoS, QString(tr("Cosmic trajectory")));

  tabWidget->addTab(tabEditor, QString(tr("Particle list editor")));

  currentTab = tabWidget->currentIndex();
  connect(tabWidget, SIGNAL(currentChanged(int)), this, SLOT(tabChanged(int)));

  QLabel *labelCopyright = new QLabel(tr("© 2014-2026 Volodymyr Vovchenko"));

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(dataLay);
  mainLayout->addSpacing(15);
  mainLayout->addWidget(tabWidget);
  mainLayout->addWidget(labelCopyright, 0, Qt::AlignRight);
  centralWidget()->setLayout(mainLayout);

  createMenus();

  QString title = "Thermal-FIST " + QString::number(ThermalFIST_VERSION_MAJOR) + "." + QString::number(ThermalFIST_VERSION_MINOR);
  if (ThermalFIST_VERSION_DEVEL != 0) title += "." + QString::number(ThermalFIST_VERSION_DEVEL);
  setWindowTitle(title);

  tab1->resetTPS();
  tab2->resetTPS();
  tabEoS->resetTPS();
  tab5->resetTPS();
  tabEditor->resetTPS();
  tabCosmicEoS->resetTPS();
}

MainWindow::~MainWindow()
{
  delete model;
  delete TPS;
}

void MainWindow::createMenus()
{
  QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
  QAction *loadAct = new QAction(tr("Load list from file..."), this);
  connect(loadAct, &QAction::triggered, this, &MainWindow::loadList);
  fileMenu->addAction(loadAct);
  QAction *loadDecaysAct = new QAction(tr("Load decays..."), this);
  connect(loadDecaysAct, &QAction::triggered, this, &MainWindow::loadDecays);
  fileMenu->addAction(loadDecaysAct);
  fileMenu->addSeparator();
  QAction *exitAct = new QAction(tr("E&xit"), this);
  exitAct->setShortcuts(QKeySequence::Quit);
  exitAct->setStatusTip(tr("Exit the application"));
  connect(exitAct, &QAction::triggered, this, &QWidget::close);
  fileMenu->addAction(exitAct);

  QMenu *viewMenu = menuBar()->addMenu(tr("&View"));
  QAction *incFontAct = new QAction(tr("Increase font size"), this);
  incFontAct->setShortcuts(QKeySequence::ZoomIn);
  connect(incFontAct, &QAction::triggered, this, &MainWindow::increaseFontSize);
  viewMenu->addAction(incFontAct);
  QAction *decFontAct = new QAction(tr("Decrease font size"), this);
  decFontAct->setShortcuts(QKeySequence::ZoomOut);
  connect(decFontAct, &QAction::triggered, this, &MainWindow::decreaseFontSize);
  viewMenu->addAction(decFontAct);
#ifdef Q_OS_WASM
  viewMenu->addSeparator();
  QAction *fullscreenAct = new QAction(tr("Toggle fullscreen"), this);
  fullscreenAct->setShortcut(QKeySequence(Qt::Key_F11));
  connect(fullscreenAct, &QAction::triggered, this, &MainWindow::toggleFullscreen);
  viewMenu->addAction(fullscreenAct);
#endif

  QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));
  QAction *aboutAct = new QAction(tr("&About Thermal-FIST"), this);
  //aboutAct->setStatusTip(tr("Show the application's About box"));
  connect(aboutAct, &QAction::triggered, this, &MainWindow::about);
  helpMenu->addAction(aboutAct);

  QAction *guideAct = new QAction(tr("Quick start guide"), this);
  connect(guideAct, &QAction::triggered, this, &MainWindow::quickstartguide);
  helpMenu->addAction(guideAct);

  QAction *docAct = new QAction(tr("Documentation"), this);
  connect(docAct, &QAction::triggered, this, &MainWindow::documentation);
  helpMenu->addAction(docAct);
}

void MainWindow::loadDecays()
{
#ifdef Q_OS_WASM
  // WASM: Use getOpenFileContent for browser file picker
  WasmFileIO::openFile(this, tr("Open file with decays"), "*.dat *.txt *.*",
    [this](const QString& sandboxPath) {
      if (sandboxPath.isEmpty())
        return;

      std::vector<std::string> decayspaths;
      decayspaths.push_back(sandboxPath.toStdString());
      TPS->LoadDecays(decayspaths);
      model->ChangeTPS(TPS);
      tab1->resetTPS();
      tab2->resetTPS();
      tab5->resetTPS();
      tabEditor->resetTPS();
      tabCosmicEoS->resetTPS();

      leList->setText(clists);
      leList->setText(leList->text() + " + " + QFileInfo(sandboxPath).fileName());
    }
  );
#else
  // Native: Use standard file dialog
  QString listpathprefix = QString(ThermalFIST_INPUT_FOLDER) + "/list";
  if (cpath.size() != 0)
    listpathprefix = cpath;
  QStringList pathdecays = QFileDialog::getOpenFileNames(this, tr("Open file with decays"), listpathprefix);
  if (pathdecays.length() > 0)
  {
    std::vector<std::string> decayspaths;
    for (int i = 0; i < pathdecays.length(); ++i)
      decayspaths.push_back(pathdecays[i].toStdString());

    m_lastDecayPaths = decayspaths;

    // Rebuild the base display name (list files + new decay files).
    QString disp;
    for (size_t i = 0; i < m_lastListPaths.size(); ++i) {
      QString fn = (i == 0) ? shortListDisplayName(QString::fromStdString(m_lastListPaths[0]))
                            : QFileInfo(QString::fromStdString(m_lastListPaths[i])).fileName();
      disp += (i == 0 ? QString() : QString(" + ")) + fn;
    }
    for (size_t idec = 0; idec < decayspaths.size(); ++idec)
      disp += " + " + QFileInfo(QString::fromStdString(decayspaths[idec])).fileName();
    clists = disp;

    // Rebuild from base with the new decays, re-applying phase shifts if enabled
    // (a plain TPS->LoadDecays would wipe the phase-shift cluster decays).
    rebuildCurrentList();
  }
#endif
}

void MainWindow::tabChanged(int newIndex)
{
  if (currentTab == tabWidget->indexOf(tabEditor) && tabEditor->haveChangesToList()) {
    tabEditor->applyChanges();
    tab1->resetTPS();
    tab2->resetTPS();
    tab5->resetTPS();
    tabEditor->resetTPS();
    tabCosmicEoS->resetTPS();
  }

  currentTab = newIndex;
}

void MainWindow::about()
{
#ifdef Q_OS_WASM
  // WASM: Use heap-allocated dialog with show() to avoid Asyncify issues
  AboutDialog *dialog = new AboutDialog(this);
  dialog->setAttribute(Qt::WA_DeleteOnClose);
  dialog->setModal(true);
  dialog->show();
#else
  AboutDialog dialog(this);
  dialog.setWindowFlags(Qt::Window);
  dialog.exec();
#endif
}

void MainWindow::quickstartguide()
{
  QDesktopServices::openUrl(QUrl("https://github.com/vlvovch/Thermal-FIST/blob/master/docs/quickstart.md"));
}

void MainWindow::documentation()
{
  QDesktopServices::openUrl(QUrl("https://fias.uni-frankfurt.de/~vovchenko/project/thermal-fist/doc/"));
}

void MainWindow::increaseFontSize()
{
  QFont font = QApplication::font();
  int sz = font.pointSize();
  font.setPointSize(sz + 1);
  QApplication::setFont(font);

  tab1->updateFontSizes();
  tab2->updateFontSizes();
}

void MainWindow::decreaseFontSize()
{
  QFont font = QApplication::font();
  int sz = font.pointSize();
  if (sz < 2)
    sz = 2;
  font.setPointSize(sz - 1);
  QApplication::setFont(font);

  tab1->updateFontSizes();
  tab2->updateFontSizes();
}

#ifdef Q_OS_WASM
void MainWindow::toggleFullscreen()
{
  // Use JavaScript Fullscreen API directly for better compatibility
  EM_ASM({
    if (document.fullscreenElement) {
      document.exitFullscreen();
    } else {
      // Try to fullscreen the canvas or body
      var elem = document.querySelector('canvas') || document.body;
      if (elem.requestFullscreen) {
        elem.requestFullscreen();
      } else if (elem.webkitRequestFullscreen) {
        elem.webkitRequestFullscreen();
      } else if (elem.mozRequestFullScreen) {
        elem.mozRequestFullScreen();
      }
    }
  });
}
#endif

void MainWindow::loadList()
{
#ifdef Q_OS_WASM
  // WASM: Use getOpenFileContent for browser file picker
  // Note: Multi-file selection not supported in WASM, load one file at a time
  WasmFileIO::openFile(this, tr("Open file with particle list"), "*.dat *.txt *.*",
    [this](const QString& listSandboxPath) {
      if (listSandboxPath.isEmpty())
        return;

      std::vector<std::string> paths = { listSandboxPath.toStdString() };

      // Try default decays next to list (in sandbox FS)
      QString decaysPath = QFileInfo(listSandboxPath).absolutePath() + "/decays.dat";
      std::vector<std::string> decays;

      if (QFileInfo(decaysPath).exists()) {
        decays.push_back(decaysPath.toStdString());
      }
      // Note: In WASM, we skip the dialog asking for decays file for simplicity
      // User can load decays separately via the Load Decays menu

      QString disp = QFileInfo(listSandboxPath).fileName();
      for (size_t idec = 0; idec < decays.size(); ++idec)
        disp += " + " + QFileInfo(QString::fromStdString(decays[idec])).fileName();
      clists = disp;

      m_lastListPaths  = paths;
      m_lastDecayPaths = decays;
      cpath = listSandboxPath;

      // Mark combos as custom (user-uploaded file)
      comboPDGEdition->setCurrentIndex(-1);
      comboListVariant->setCurrentIndex(-1);

      rebuildCurrentList();   // build base + (re)apply phase shifts + display + reset
    }
  );
#else
  // Native: Use standard file dialog
  QString listpathprefix = QString(ThermalFIST_INPUT_FOLDER) + "/list";
  if (cpath.size() != 0)
    listpathprefix = cpath;
  //QString path = QFileDialog::getOpenFileName(this, tr("Open file with particle list"), listpathprefix);
  QStringList pathlist = QFileDialog::getOpenFileNames(this, tr("Open file(s) with particle list"), listpathprefix);
  if (pathlist.length() > 0)
  {
    //*TPS = ThermalParticleSystem(path.toStdString());
    std::vector<std::string> paths(0);
    for (int i = 0; i < pathlist.length(); ++i)
      paths.push_back(pathlist[i].toStdString());
    //*TPS = ThermalParticleSystem(paths);
    //*TPS = ThermalParticleSystem(
    //  { path.toStdString() },
    //  { "" },
    //  { ThermalParticleSystem::flag_noexcitednuclei, ThermalParticleSystem::flag_nonuclei }
    //);

    std::vector<std::string> decays;
    QStringList decpath;
    decpath.push_back(QFileInfo(pathlist[0]).absolutePath() + "/decays.dat");
    //if (!TPS->CheckDecayChannelsAreSpecified() &&  !QFileInfo(decpath).exists()) {
    if (!QFileInfo(decpath[0]).exists()) {

      QMessageBox::StandardButton reply;
      reply = QMessageBox::question(this,
        "Decays",
        "Decays file was not found at `decays.dat`. Would you like to load decays from another file?",
        QMessageBox::Yes | QMessageBox::No);
      if (reply == QMessageBox::Yes) {
        decpath = QFileDialog::getOpenFileNames(this, tr("Open file with decays"), decpath[0]);
        if (decpath.length() > 0)
        {
          for (int i = 0; i < decpath.length(); ++i)
            decays.push_back(decpath[i].toStdString());
        }
      }
    }
    else {
      decays.push_back(decpath[0].toStdString());
    }
    // Build the base display name (list + decay file names, no particle count).
    QString disp = shortListDisplayName(pathlist[0]);
    for (int il = 1; il < pathlist.size(); ++il)
      disp += " + " + QFileInfo(pathlist[il]).fileName();
    for (size_t idec = 0; idec < decays.size(); ++idec)
      disp += " + " + QFileInfo(QString::fromStdString(decays[idec])).fileName();
    clists = disp;

    m_lastListPaths  = paths;
    m_lastDecayPaths = decays;
    cpath = pathlist[0];

    // Mark combos as custom (user-selected file)
    comboPDGEdition->setCurrentIndex(-1);
    comboListVariant->setCurrentIndex(-1);

    rebuildCurrentList();   // build base + (re)apply phase shifts + display + reset
  }
#endif
}

QString MainWindow::shortListDisplayName(const QString& fullPath)
{
  // Extract "PDG2025/list-withnuclei.dat" from a full path like
  // "/path/to/input/list/PDG2025/list-withnuclei.dat"
  QFileInfo fi(fullPath);
  QString dirName = fi.dir().dirName();   // e.g. "PDG2025" or "default"
  QString fileName = fi.fileName();       // e.g. "list-withnuclei.dat"
  if (dirName.isEmpty() || dirName == ".")
    return fileName;
  return dirName + "/" + fileName;
}

void MainWindow::applyParticleList(const std::vector<std::string>& listPaths,
                                   const std::vector<std::string>& decayPaths,
                                   const QString& displayName)
{
  m_lastListPaths  = listPaths;
  m_lastDecayPaths = decayPaths;
  clists = displayName;
  if (!listPaths.empty())
    cpath = QString::fromStdString(listPaths[0]);
  rebuildCurrentList();
}

void MainWindow::resetAllTabs()
{
  // Each tab's resetTPS() recomputes its EV/vdW interaction matrix over the
  // current particle list (QvdWParameters::GetParameters), so the clusters
  // (B=0 mesons) always carry meson parameters. Enabling/disabling phase shifts
  // only flips the cluster densities (the list is unchanged), so no manual
  // interaction re-apply is ever needed.
  tab1->resetTPS();
  tab2->resetTPS();
  tabEoS->resetTPS();
  tab5->resetTPS();
  tabEditor->resetTPS();
  if (!cpath.isEmpty())
    tabEditor->setListPath(cpath);
  tabCosmicEoS->resetTPS();
}

void MainWindow::refreshListDisplay()
{
  QString t = clists;
  const bool present = TPS && PhaseShifts::CountPhaseShiftDensities(*TPS) > 0;
  if (present && chkPhaseShifts && chkPhaseShifts->isChecked() && !m_phaseShiftConf.isEmpty())
    t += " + PS[" + QFileInfo(m_phaseShiftConf).fileName() + "]";
  else if (present)
    t += " (phase shifts off)";   // clusters still in the list but disabled
  if (TPS)
    t += " [" + QString::number(TPS->Particles().size()) + " particles]";
  leList->setText(t);
}

void MainWindow::applyPhaseShiftsIfEnabled()
{
  if (!chkPhaseShifts || !chkPhaseShifts->isChecked()) return;
  if (m_phaseShiftConf.isEmpty()) return;
  try {
    PhaseShifts::AddPhaseShiftChannelsFromFile(*TPS, m_phaseShiftConf.toStdString());
  } catch (const std::exception& e) {
    QMessageBox::warning(this, tr("Phase shifts"),
      tr("Failed to apply phase-shift config:\n%1\n\nDisabling phase shifts.").arg(e.what()));
    chkPhaseShifts->blockSignals(true);
    chkPhaseShifts->setChecked(false);
    chkPhaseShifts->blockSignals(false);
  }
}

void MainWindow::rebuildCurrentList()
{
  // Rebuild the base list from the stored selection, then (re)apply the
  // phase-shift channels on top if enabled. Rebuilding from scratch is how
  // disabling/swapping works, since AddPhaseShiftChannel mutates the list.
  *TPS = ThermalParticleSystem(m_lastListPaths, m_lastDecayPaths);
  model->ChangeTPS(TPS);
  applyPhaseShiftsIfEnabled();
  refreshListDisplay();
  resetAllTabs();
}

void MainWindow::onPhaseShiftToggled()
{
  const bool on = chkPhaseShifts->isChecked();
  if (PhaseShifts::CountPhaseShiftDensities(*TPS) > 0) {
    // Clusters already in the list: just switch their density on/off. No rebuild,
    // so the EV/vdW matrix (which already includes them) is untouched.
    PhaseShifts::SetPhaseShiftsEnabled(*TPS, on);
    refreshListDisplay();
    resetAllTabs();
  } else if (on) {
    // First time on: add the channels to the list (one rebuild).
    rebuildCurrentList();
  }
}

void MainWindow::loadPhaseShiftConf()
{
  QString prefix = QString(ThermalFIST_INPUT_FOLDER) + "/list/phaseshifts";
  if (!m_phaseShiftConf.isEmpty())
    prefix = m_phaseShiftConf;
  QString path = QFileDialog::getOpenFileName(this, tr("Open phase-shift config file"), prefix,
                   tr("Config files (*.conf *.txt *.dat);;All files (*)"));
  if (path.isEmpty())
    return;
  m_phaseShiftConf = path;
  // A new config means different channels: rebuild from the base list so the old
  // clusters are dropped and the new ones added (with the box checked).
  chkPhaseShifts->blockSignals(true);
  chkPhaseShifts->setChecked(true);
  chkPhaseShifts->blockSignals(false);
  rebuildCurrentList();
}

void MainWindow::updateListVariants()
{
  QString edition = comboPDGEdition->currentData().toString();

  comboListVariant->blockSignals(true);
  comboListVariant->clear();

  // Common variants available in all editions
  comboListVariant->addItem(tr("Default (with nuclei)"),  QString("list-withnuclei.dat"));
  comboListVariant->addItem(tr("No nuclei"),              QString("list.dat"));
  comboListVariant->addItem(tr("With charm"),             QString("list-withcharm.dat"));
  comboListVariant->addItem(tr("With excited nuclei"),    QString("list-withexcitednuclei.dat"));

  // Isospin symmetric only available for PDG2020 and PDG2025
  if (edition == "PDG2020" || edition == "PDG2025") {
    comboListVariant->addItem(tr("Isospin symmetric"), QString("list-isospin-symmetric.dat"));
  }

  comboListVariant->setCurrentIndex(0);
  comboListVariant->blockSignals(false);
}

void MainWindow::switchParticleList()
{
  QString edition = comboPDGEdition->currentData().toString();
  QString listFile = comboListVariant->currentData().toString();
  if (edition.isEmpty() || listFile.isEmpty())
    return;

  QString displayName = edition + "/" + listFile + " + decays.dat";

#ifdef Q_OS_WASM
  // Extract the selected list from embedded resources to sandbox
  QString sandboxDir = WasmFileIO::getSandboxTempDir() + "/default";
  QDir().mkpath(sandboxDir);

  // Extract list file
  QString resourcePath = ":/data/lists/" + edition + "/" + listFile;
  QFile listRes(resourcePath);
  QString listPath = sandboxDir + "/" + listFile;
  if (listRes.open(QIODevice::ReadOnly)) {
    QFile listOut(listPath);
    if (listOut.open(QIODevice::WriteOnly)) {
      listOut.write(listRes.readAll());
      listOut.close();
    }
    listRes.close();
  } else {
    qWarning() << "Could not open resource:" << resourcePath;
    return;
  }

  // Extract decays file for this edition
  QString decaysResPath = ":/data/lists/" + edition + "/decays.dat";
  QFile decaysRes(decaysResPath);
  QString decaysPath = sandboxDir + "/decays.dat";
  if (decaysRes.open(QIODevice::ReadOnly)) {
    QFile decaysOut(decaysPath);
    if (decaysOut.open(QIODevice::WriteOnly)) {
      decaysOut.write(decaysRes.readAll());
      decaysOut.close();
    }
    decaysRes.close();
  }

  std::vector<std::string> paths = { listPath.toStdString() };
  std::vector<std::string> decays;
  if (QFileInfo(decaysPath).exists())
    decays.push_back(decaysPath.toStdString());
#else
  // Native: construct path from input folder
  QString listDir = QString(ThermalFIST_INPUT_FOLDER) + "/list/" + edition + "/";
  QString listPath = listDir + listFile;
  QString decaysPath = listDir + "decays.dat";

  std::vector<std::string> paths = { listPath.toStdString() };
  std::vector<std::string> decays;
  if (QFileInfo(decaysPath).exists())
    decays.push_back(decaysPath.toStdString());
#endif

  applyParticleList(paths, decays, displayName);
}
