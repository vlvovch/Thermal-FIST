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
	cpath = QString(INPUT_FOLDER) + "/list/PDG2014/list-withnuclei.dat";

	QString listpath = cpath;
	TPS = new ThermalParticleSystem(listpath.toStdString());
	model = new ThermalModelIdeal(TPS);

	QWidget *centW = new QWidget;
	setCentralWidget(centW);

	//QVBoxLayout *dirTabLay1 = new QVBoxLayout();
	QHBoxLayout *dataLay = new QHBoxLayout();
	QLabel *labelData = new QLabel(tr("Particle list file:"));
	dataLay->setAlignment(Qt::AlignLeft);
	leDatabase = new QLineEdit("");//QApplication::applicationDirPath());
	leDatabase->setReadOnly(true);
	if (TPS->Particles().size() > 0)
		leDatabase->setText(listpath);

	buttonLoad = new QPushButton(tr("Load particle list..."));
	connect(buttonLoad, SIGNAL(clicked()), this, SLOT(loadDatabase()));

	buttonLoadDecays = new QPushButton(tr("Load decays..."));
	connect(buttonLoadDecays, SIGNAL(clicked()), this, SLOT(loadDecays()));

	dataLay->addWidget(labelData);
	dataLay->addWidget(leDatabase, 1);
	dataLay->addWidget(buttonLoad);
	dataLay->addWidget(buttonLoadDecays);

	tab1 = new ModelTab(NULL, model);
	tab1->resetTPS();

	tab2 = new FitToExperimentTab(NULL, model);
  tab1->setFitTab(tab2);

  tabEoS = new EquationOfStateTab(NULL, model);

	tab5 = new EventGeneratorTab(NULL, model);
  tabEditor = new ListEditorTab(NULL, model);

  

	tabWidget = new QTabWidget();
	tabWidget->addTab(tab1, QString(tr("Thermal model")));

	tabWidget->addTab(tab2, QString(tr("Thermal fits")));

  tabWidget->addTab(tabEoS, QString(tr("Equation of state")));

	tabWidget->addTab(tab5, QString(tr("Event generator")));

  tabWidget->addTab(tabEditor, QString(tr("Particle list editor")));

  currentTab = tabWidget->currentIndex();
  connect(tabWidget, SIGNAL(currentChanged(int)), this, SLOT(tabChanged(int)));

	QLabel *labelCopyright = new QLabel(tr("© 2014-2019 Volodymyr Vovchenko"));

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

//#ifndef QT_NO_CONTEXTMENU
//void MainWindow::contextMenuEvent(QContextMenuEvent *event)
//{
//  QMenu menu(this);
//  menu.exec(event->globalPos());
//}
//#endif // QT_NO_CONTEXTMENU

void MainWindow::createMenus()
{
  QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
  QAction *loadAct = new QAction(tr("Load particle list..."), this);
  connect(loadAct, &QAction::triggered, this, &MainWindow::loadDatabase);
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
}

void MainWindow::loadDecays()
{
	QString listpathprefix = QString(INPUT_FOLDER) + "/list";
	if (leDatabase->text().size() != 0)
		listpathprefix = leDatabase->text();
	QString path = QFileDialog::getOpenFileName(this, tr("Open file with decays"), listpathprefix);
	if (path.length() > 0)
	{
		TPS->LoadDecays(path.toStdString());
		model->ChangeTPS(TPS);
		tab1->resetTPS();
		tab2->resetTPS();
		tab5->resetTPS();
    tabEditor->resetTPS();
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

void MainWindow::loadDatabase()
{
	QString listpathprefix = QString(INPUT_FOLDER) + "/list";
	if (leDatabase->text().size() != 0)
		listpathprefix = leDatabase->text();
	QString path = QFileDialog::getOpenFileName(this, tr("Open file with particle list"), listpathprefix);
	if (path.length() > 0)
	{
		*TPS = ThermalParticleSystem(path.toStdString());
		model->ChangeTPS(TPS);
		leDatabase->setText(path);
		tab1->resetTPS();
		tab2->resetTPS();

    tabEoS->resetTPS();

		tab5->resetTPS();

    tabEditor->resetTPS();

		cpath = path;
	}
}
