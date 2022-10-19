/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "mainwindow.h"

#include <QLayout>
#include <QFileDialog>
#include <QLabel>
#include <QApplication>
#include <QMessageBox>
#include <QElapsedTimer>
#include <QDebug>


#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonicalCharm.h"

#include "aboutdialog.h"

using namespace thermalfist;

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
{
  //cpath = QString(ThermalFIST_INPUT_FOLDER) + "/list/PDG2020/list-withnuclei.dat";
  cpath = QString(ThermalFIST_DEFAULT_LIST_FILE);

  QString listpath = cpath;
  TPS = new ThermalParticleSystem(listpath.toStdString());

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
  tab1->resetTPS();

  tab2 = new FitToExperimentTab(NULL, model);
  tab1->setFitTab(tab2);

  tabEoS = new EquationOfStateTab(NULL, model);

  tab5 = new EventGeneratorTab(NULL, model);
  tabEditor = new ListEditorTab(NULL, model);
  tabEditor->setListPath(cpath);

  

  tabWidget = new QTabWidget();
  tabWidget->addTab(tab1, QString(tr("Thermal model")));

  tabWidget->addTab(tab2, QString(tr("Thermal fits")));

  tabWidget->addTab(tabEoS, QString(tr("Equation of state")));

  tabWidget->addTab(tab5, QString(tr("Event generator")));

  tabWidget->addTab(tabEditor, QString(tr("Particle list editor")));

  currentTab = tabWidget->currentIndex();
  connect(tabWidget, SIGNAL(currentChanged(int)), this, SLOT(tabChanged(int)));

  QLabel *labelCopyright = new QLabel(tr("© 2014-2020 Volodymyr Vovchenko"));

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
}

void MainWindow::tabChanged(int newIndex)
{
  if (currentTab == tabWidget->indexOf(tabEditor) && tabEditor->haveChangesToList()) {
    tabEditor->applyChanges();
    tab1->resetTPS();
    tab2->resetTPS();
    tab5->resetTPS();
    tabEditor->resetTPS();
  }

  currentTab = newIndex;
}

void MainWindow::about()
{
  AboutDialog dialog(this);
  dialog.setWindowFlags(Qt::Window);
  dialog.exec();
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

void MainWindow::loadList()
{
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

    cpath = pathlist[0];
  }
}
