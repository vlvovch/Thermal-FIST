#include "mainwindow.h"

#include <QLayout>
#include <QFileDialog>
#include <QLabel>
#include <QApplication>
#include <QMessageBox>
#include <QElapsedTimer>
#include <QDebug>

//#include "modeltab.h"
#include "ThermalFISTConfig.h"
#include "HRGBase/ThermalModelIdeal.h"
#include "HRGEV/ThermalModelEVDiagonal.h"
#include "HRGBase/ThermalModelCanonicalStrangeness.h"
#include "HRGBase/ThermalModelCanonicalCharm.h"


MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	cpath = QString(INPUT_FOLDER) + "/list/PDG2014/list-withnuclei.dat";

	//QString listpath = QString(INPUT_FOLDER) + "/list/PDG2014/list-multibaryons.dat";
	//QString listpath = QDir::currentPath() + "/input/list/PDG2014/list.dat";
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

#ifdef USE_MINUIT
	tab2 = new FitToExperimentTab(NULL, model);
#endif

	tab3 = new EnergyDependenceTab(NULL, model);
	tab4 = new ContourPlotTab(NULL, model);
	tab5 = new EventGeneratorTab(NULL, model);

	tabWidget = new QTabWidget();
	tabWidget->addTab(tab1, QString(tr("Thermal model")));

#ifdef USE_MINUIT
	tabWidget->addTab(tab2, QString(tr("Thermal fits")));
#endif

	tabWidget->addTab(tab5, QString(tr("Event generator")));

	QLabel *labelCopyright = new QLabel(tr("© 2014-2018 Volodymyr Vovchenko"));

	QVBoxLayout *mainLayout = new QVBoxLayout;
	mainLayout->addLayout(dataLay);
	mainLayout->addSpacing(15);
	mainLayout->addWidget(tabWidget);
	mainLayout->addWidget(labelCopyright, 0, Qt::AlignRight);
	centralWidget()->setLayout(mainLayout);

	QString title = "Thermal-FIST " + QString::number(ThermalFIST_VERSION_MAJOR) + "." + QString::number(ThermalFIST_VERSION_MINOR);
	if (ThermalFIST_VERSION_DEVEL != 0) title += "." + QString::number(ThermalFIST_VERSION_DEVEL);
	//title = "Thermal-FIST 0.3";
	setWindowTitle(title);
}

MainWindow::~MainWindow()
{
	delete model;
	delete TPS;
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
#ifdef USE_MINUIT
		tab2->resetTPS();
#endif
		tab5->resetTPS();
	}
}

void MainWindow::loadDatabase()
{
	QString listpathprefix = QString(INPUT_FOLDER) + "/list";
	if (leDatabase->text().size() != 0)
		listpathprefix = leDatabase->text();
	QString path = QFileDialog::getOpenFileName(this, tr("Open file with particle database"), listpathprefix);
	if (path.length() > 0)
	{
		*TPS = ThermalParticleSystem(path.toStdString());
		model->ChangeTPS(TPS);
		leDatabase->setText(path);
		tab1->resetTPS();
#ifdef USE_MINUIT
		tab2->resetTPS();
#endif
		tab5->resetTPS();

		cpath = path;
	}
}
