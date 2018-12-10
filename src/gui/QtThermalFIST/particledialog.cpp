/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "particledialog.h"

#include <QLayout>
#include <QHeaderView>
#include <QTimer>
#include <QLabel>
#include <QApplication>
#include <QDebug>

#include "SpectralFunctionDialog.h"

#include "HRGBase/ThermalModelBase.h"

using namespace thermalfist;

ParticleDialog::ParticleDialog(QWidget *parent, ThermalModelBase *mod, int ParticleID) :
    QDialog(parent), model(mod), pid(ParticleID)
{
    QVBoxLayout *layout = new QVBoxLayout();

    QFont font = QApplication::font();
    font.setPointSize(12);
    QLabel *labInfo = new QLabel(tr("Information:"));
    labInfo->setFont(font);
    data = new QTextEdit();
    data->setReadOnly(true);
    data->setPlainText(GetParticleInfo());
    data->setMinimumHeight(350);
    data->setMinimumWidth(300);

    QLabel *labDecays = new QLabel(QString(model->TPS()->Particles()[pid].Name().c_str()) + tr(" decays:"));
    labDecays->setFont(font);
    myModel = new DecayTableModel(0, &model->TPS()->Particle(pid), model->TPS());
    tableDecays = new QTableView();
    tableDecays->setModel(myModel);
    tableDecays->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);

    QHBoxLayout *layoutButtons = new QHBoxLayout();
    layoutButtons->setAlignment(Qt::AlignLeft);

    layout->addWidget(labInfo);
    layout->addWidget(data);
    if (model->TPS()->Particles()[pid].Decays().size()>0) {
        layout->addWidget(labDecays);
        layout->addWidget(tableDecays);
    }

		if (model->TPS()->Particles()[pid].ResonanceWidth() > 0.) {
			buttonSpectralFunction = new QPushButton(tr("Spectral function..."));
			connect(buttonSpectralFunction, SIGNAL(clicked()), this, SLOT(showSpectralFunction()));
			layout->addWidget(buttonSpectralFunction, 0, Qt::AlignRight);
		}

    QLabel *labProd = new QLabel(tr("Particle production"));
    labProd->setFont(font);

    QLabel *labelPrimDensity = new QLabel(tr("Primordial density = "));
    if (model->IsCalculated()) labelPrimDensity->setText(labelPrimDensity->text() +
                                                      QString::number(model->GetParticlePrimordialDensity(pid)) +
                                                      " fm<sup>-3</sup>");

    QLabel *labelPrimMultiplicity = new QLabel(tr("Primordial multiplicity = "));
    if (model->IsCalculated()) labelPrimMultiplicity->setText(labelPrimMultiplicity->text() +
                                                      QString::number(model->GetParticlePrimordialDensity(pid) * model->Volume()));

    QLabel *labelTotalMultiplicity = new QLabel(tr("Total multiplicity = "));
    if (model->IsCalculated()) labelTotalMultiplicity->setText(labelTotalMultiplicity->text() +
                                                      QString::number(model->GetParticleTotalDensity(pid) * model->Volume()));

    tableSources = new QTableWidget();
    tableSources->setColumnCount(3);
    tableSources->setRowCount(model->TPS()->Particles()[pid].DecayContributions().size()+1);
    tableSources->verticalHeader()->hide();
    tableSources->setHorizontalHeaderItem(0, new QTableWidgetItem(tr("Source")));
    tableSources->setHorizontalHeaderItem(1, new QTableWidgetItem(tr("Multiplicity")));
    tableSources->setHorizontalHeaderItem(2, new QTableWidgetItem(tr("Fraction (%)")));

    if (model->IsCalculated()) {

        std::vector< std::pair<double, int> > sources = model->TPS()->Particles()[pid].DecayContributions();
        for(int i=0;i<sources.size();++i) sources[i].first *= model->GetParticlePrimordialDensity(sources[i].second);

        qSort(sources.begin(), sources.end());

        tableSources->setItem(0, 0, new QTableWidgetItem(tr("Direct thermal")));
        tableSources->setItem(0, 1, new QTableWidgetItem(QString::number(model->GetParticlePrimordialDensity(pid) * model->Volume())));
        tableSources->setItem(0, 2, new QTableWidgetItem(QString::number(model->GetParticlePrimordialDensity(pid) / model->GetParticleTotalDensity(pid) * 100.)));

				
				for (int j = 0; j < 3; ++j) {
					QTableWidgetItem *item = tableSources->item(0, j);
					item->setFlags(item->flags() &  ~Qt::ItemIsEditable);
				}

        for(int i=0;i<sources.size();++i)
          {
                int tindex = sources.size() - i - 1;
                int decaypartid = sources[tindex].second;
                tableSources->setItem(i+1, 0, new QTableWidgetItem(tr("Decays from thermal ") + QString(model->TPS()->Particles()[decaypartid].Name().c_str())));
                tableSources->setItem(i+1, 1, new QTableWidgetItem(QString::number(sources[tindex].first * model->Volume())));
                tableSources->setItem(i+1, 2, new QTableWidgetItem(QString::number(sources[tindex].first / model->GetParticleTotalDensity(pid) * 100.)));
								
								for (int j = 0; j < 3; ++j) {
									QTableWidgetItem *item = tableSources->item(i+1, j);
									item->setFlags(item->flags() &  ~Qt::ItemIsEditable);
								}
					}
    }

    tableSources->resizeColumnsToContents();

    QVBoxLayout *layoutProd = new QVBoxLayout();
    layoutProd->setAlignment(Qt::AlignTop);
    layoutProd->addWidget(labProd);
    layoutProd->addWidget(labelPrimDensity);
    layoutProd->addWidget(labelPrimMultiplicity);
    layoutProd->addWidget(labelTotalMultiplicity);
    layoutProd->addWidget(tableSources);

    QHBoxLayout *mainLayout = new QHBoxLayout();

    mainLayout->addLayout(layout);
    mainLayout->addLayout(layoutProd);

    setLayout(mainLayout);

    this->setWindowTitle(QString(model->TPS()->Particles()[pid].Name().c_str()));

    QTimer::singleShot(0, this, SLOT(checkFixTableSize()));
}

QString ParticleDialog::GetParticleInfo() {
    QString ret = "";

    ret += tr("Name") + "\t\t= ";
    ret += QString(model->TPS()->Particles()[pid].Name().c_str()) + "\r\n";

    ret += tr("PDG ID") + "\t\t= ";
    ret += QString::number(model->TPS()->Particles()[pid].PdgId()) + "\r\n";

    ret += tr("Mass") + "\t\t= ";
    ret += QString::number(model->TPS()->Particles()[pid].Mass()*1.e3) + " " + tr("MeV") + "\r\n";

    ret += tr("Type") + "\t\t= ";
    if (model->TPS()->Particles()[pid].BaryonCharge()==0) ret += tr("Meson") + "\r\n";
    else if (model->TPS()->Particles()[pid].BaryonCharge()>0) ret += tr("Baryon") + "\r\n";
    else ret += "Anti-Baryon\r\n";

    ret +=  tr("Stable?") + "\t\t= ";
    if (model->TPS()->Particles()[pid].IsStable()) ret += tr("Yes") + "\r\n";
    else ret += "No\r\n";

    ret += tr("Neutral?") + "\t\t= ";
    if (model->TPS()->Particles()[pid].IsNeutral()) ret +=  tr("Yes") + "\r\n";
    else ret += "No\r\n";

    ret += tr("Spin degeneracy") + "\t= ";
    ret += QString::number(model->TPS()->Particles()[pid].Degeneracy()) + "\r\n";

    ret += tr("Electric charge") + "\t= ";
    ret += QString::number(model->TPS()->Particles()[pid].ElectricCharge()) + "\r\n";

    ret += tr("Strangeness") + "\t= ";
    ret += QString::number(model->TPS()->Particles()[pid].Strangeness()) + "\r\n";

    ret += tr("Charm") + "\t\t= ";
    ret += QString::number(model->TPS()->Particles()[pid].Charm()) + "\r\n";

    ret += tr("Abs. strangeness") + "\t= ";
    ret += QString::number(model->TPS()->Particles()[pid].AbsoluteStrangeness()) + "\r\n";

    ret += tr("Charm content") + "\t= ";
    ret += QString::number(model->TPS()->Particles()[pid].AbsoluteCharm()) + "\r\n";

		ret += tr("Isospin") + "\t= ";
		ret += QString::number(2*model->TPS()->Particles()[pid].ElectricCharge() - model->TPS()->Particles()[pid].BaryonCharge() - model->TPS()->Particles()[pid].Strangeness() - model->TPS()->Particles()[pid].Charm()) + "\r\n";

		ret += tr("I3") + "\t= ";
		ret += QString::number(model->TPS()->Particles()[pid].ElectricCharge() - (model->TPS()->Particles()[pid].BaryonCharge() + model->TPS()->Particles()[pid].Strangeness() + model->TPS()->Particles()[pid].Charm())/2.) + "\r\n";

    if (model->TPS()->Particles()[pid].ResonanceWidth()>1e-5) {
        ret += tr("Decay width") + "\t= ";
        ret += QString::number(model->TPS()->Particles()[pid].ResonanceWidth()*1.e3) + " " + tr("MeV") + "\r\n";

        if (model->TPS()->Particles()[pid].DecayThresholdMass()>1e-5) {
            ret += tr("Threshold mass") + "\t= ";
            ret += QString::number(model->TPS()->Particles()[pid].DecayThresholdMass()*1.e3) + " " + tr("MeV");
						ret += "\r\n";
        }
    }

    return ret;
}

void ParticleDialog::checkFixTableSize() {
    int tle = 0;
    for(int i=0;i<tableDecays->horizontalHeader()->count();++i)
        tle += tableDecays->horizontalHeader()->sectionSize(i);
    tableDecays->setMinimumWidth(tle+23);
    tle = 0;
    for(int i=0;i<tableSources->horizontalHeader()->count();++i)
        tle += tableSources->horizontalHeader()->sectionSize(i);
    tableSources->setMinimumWidth(tle+23);
}

void ParticleDialog::addDecay() {
    myModel->addDecay();
}

void ParticleDialog::removeDecay() {
    QModelIndexList selectedList = tableDecays->selectionModel()->selectedRows();
    //qDebug() << selectedList.count() << "\n";
    for(unsigned int i = 0; i < selectedList.count(); ++i)
        myModel->removeDecay(selectedList.at(i).row());
}

void ParticleDialog::addColumn() {
    myModel->addColumn();
    checkFixTableSize();
}

void ParticleDialog::removeColumn() {
    myModel->removeColumn();
}

void ParticleDialog::showSpectralFunction()
{
	SpectralFunctionDialog dialog(this, &model->TPS()->Particle(pid), model->Parameters().T, model->ChemicalPotential(pid), static_cast<int>(model->TPS()->ResonanceWidthIntegrationType() == ThermalParticle::eBW));
	dialog.setWindowFlags(Qt::Window);
	dialog.setMinimumSize(QSize(800, 400));
	dialog.exec();
}
