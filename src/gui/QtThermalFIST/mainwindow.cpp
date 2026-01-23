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
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonicalCharm.h"

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
  //cpath = QString(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list-withnuclei.dat";
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
  QLabel *labelData = new QLabel(tr("Particle list file:"));
  dataLay->setAlignment(Qt::AlignLeft);
  leList = new QLineEdit("");//QApplication::applicationDirPath());
  leList->setReadOnly(true);
  if (TPS->Particles().size() > 0)
    leList->setText(listpath + " + decays.dat");

  buttonLoad = new QPushButton(tr("Load particle list..."));
  connect(buttonLoad, SIGNAL(clicked()), this, SLOT(loadList()));

  buttonLoadDecays = new QPushButton(tr("Load decays..."));
  connect(buttonLoadDecays, SIGNAL(clicked()), this, SLOT(loadDecays()));

  dataLay->addWidget(labelData);
  dataLay->addWidget(leList, 1);
  dataLay->addWidget(buttonLoad);
  dataLay->addWidget(buttonLoadDecays);

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

  QLabel *labelCopyright = new QLabel(tr("© 2014-2025 Volodymyr Vovchenko"));

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
  QAction *loadAct = new QAction(tr("Load particle list..."), this);
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
  if (leList->text().size() != 0)
    listpathprefix = leList->text();
  QStringList pathdecays = QFileDialog::getOpenFileNames(this, tr("Open file with decays"), listpathprefix);
  if (pathdecays.length() > 0)
  {
    std::vector<std::string> decayspaths;
    for (int i = 0; i < pathdecays.length(); ++i)
      decayspaths.push_back(pathdecays[i].toStdString());
    TPS->LoadDecays(decayspaths);
    model->ChangeTPS(TPS);
    tab1->resetTPS();
    tab2->resetTPS();
    tab5->resetTPS();
    tabEditor->resetTPS();
    tabCosmicEoS->resetTPS();

    leList->setText(clists);

    if (QFileInfo(QString::fromStdString(decayspaths[0])).dir().absolutePath() !=
      QFileInfo(clists).dir().absolutePath())
      leList->setText(leList->text() + " + " + QString::fromStdString(decayspaths[0]));
    else
      leList->setText(leList->text() + " + " + QFileInfo(QString::fromStdString(decayspaths[0])).fileName());

    for (int idec = 1; idec < decayspaths.size(); ++idec) {
      leList->setText(leList->text() + " + " + QFileInfo(QString::fromStdString(decayspaths[idec])).fileName());
    }
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

      *TPS = ThermalParticleSystem(paths, decays);
      model->ChangeTPS(TPS);

      leList->setText(QFileInfo(listSandboxPath).fileName());
      clists = leList->text();

      if (decays.size() > 0) {
        leList->setText(leList->text() + " + " + QFileInfo(QString::fromStdString(decays[0])).fileName());
      }

      tab1->resetTPS();
      tab2->resetTPS();
      tabEoS->resetTPS();
      tab5->resetTPS();
      tabEditor->resetTPS();
      tabEditor->setListPath(listSandboxPath);
      tabCosmicEoS->resetTPS();

      cpath = listSandboxPath;
    }
  );
#else
  // Native: Use standard file dialog
  QString listpathprefix = QString(ThermalFIST_INPUT_FOLDER) + "/list";
  if (leList->text().size() != 0)
    listpathprefix = leList->text();
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
    *TPS = ThermalParticleSystem(paths, decays);

    //TPS->SetSortMode(ThermalParticleSystem::SortByBaryonAndMassAndPDG);
    model->ChangeTPS(TPS);
    leList->setText(pathlist[0]);
    if (pathlist.size() > 1) {
      for (int il = 1; il < pathlist.size(); ++il) {
        leList->setText(leList->text() + " + " + QFileInfo(pathlist[il]).fileName());
      }
    }
    clists = leList->text();

    if (decays.size() > 0) {

      if (QFileInfo(QString::fromStdString(decays[0])).dir().absolutePath() !=
        QFileInfo(clists).dir().absolutePath())
        leList->setText(leList->text() + " + " + QString::fromStdString(decays[0]));
      else
        leList->setText(leList->text() + " + " + QFileInfo(QString::fromStdString(decays[0])).fileName());

      for (int idec = 1; idec < decays.size(); ++idec) {
        leList->setText(leList->text() + " + " + QFileInfo(QString::fromStdString(decays[idec])).fileName());
      }
    }

    tab1->resetTPS();
    tab2->resetTPS();

    tabEoS->resetTPS();

    tab5->resetTPS();

    tabEditor->resetTPS();
    tabEditor->setListPath(pathlist[0]);

    tabCosmicEoS->resetTPS();

    cpath = pathlist[0];
  }
#endif
}
