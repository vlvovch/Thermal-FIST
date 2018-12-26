/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "decayseditor.h"
#include <QLayout>
#include <QHeaderView>
#include <QTimer>
#include <QDebug>

using namespace thermalfist;

decayseditor::decayseditor(QWidget *parent, ThermalParticle *Particle, ThermalParticleSystem *TPS) :
    QDialog(parent)
{
    QVBoxLayout *layout = new QVBoxLayout();

    myModel = new DecayEditorTableModel(0, Particle, TPS);
    tableDecays = new QTableView();
    tableDecays->setModel(myModel);
    //tableDecays->setSelectionBehavior(QAbstractItemView::SelectRows);
    //tableDecays->setSelectionBehavior(QAbstractItemView::SelectRows);
    tableDecays->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    //tableDecays->setEditTriggers(QAbstractItemView::EditTrigger);
    //tableDecays->resizeColumnsToContents();

    QHBoxLayout *layoutButtons = new QHBoxLayout();
    layoutButtons->setAlignment(Qt::AlignLeft);

    buttonAddDecay = new QPushButton(tr("Add decay"));
    connect(buttonAddDecay, SIGNAL(clicked()), this, SLOT(addDecay()));
    buttonRemoveDecay = new QPushButton(tr("Remove decay"));
    connect(buttonRemoveDecay, SIGNAL(clicked()), this, SLOT(removeDecay()));

    buttonAddDaughterColumn = new QPushButton(tr("Add column"));
    connect(buttonAddDaughterColumn, SIGNAL(clicked()), this, SLOT(addColumn()));
    buttonRemoveDaughterColumn = new QPushButton(tr("Remove column"));
    connect(buttonRemoveDaughterColumn, SIGNAL(clicked()), this, SLOT(removeColumn()));

    layoutButtons->addWidget(buttonAddDecay);
    layoutButtons->addWidget(buttonRemoveDecay);
    layoutButtons->addWidget(buttonAddDaughterColumn);
    layoutButtons->addWidget(buttonRemoveDaughterColumn);

    layout->addWidget(tableDecays);
    layout->addLayout(layoutButtons);

    setLayout(layout);

    this->setWindowTitle(QString("%1 decays").arg(Particle->Name().c_str()));

    QTimer::singleShot(0, this, SLOT(checkFixTableSize()));
}

void decayseditor::checkFixTableSize() {
    int tle = 0;
    for(int i=0;i<tableDecays->horizontalHeader()->count();++i)
        tle += tableDecays->horizontalHeader()->sectionSize(i);
    
    if (tableDecays->size().width() - tle < 23)
    {
        int diff = tle - tableDecays->size().width() + 23;
        int cwi = width();
        resize(cwi+diff, height());
        setMinimumWidth(cwi+diff);
    }
}

void decayseditor::addDecay() {
    myModel->addDecay();
}

void decayseditor::removeDecay() {
    QModelIndexList selectedList = tableDecays->selectionModel()->selectedRows();
    qDebug() << selectedList.count() << "\n";
    for(unsigned int i = 0; i < selectedList.count(); ++i)
        myModel->removeDecay(selectedList.at(i).row());
}

void decayseditor::addColumn() {
    myModel->addColumn();
    checkFixTableSize();
}

void decayseditor::removeColumn() {
    myModel->removeColumn();
}
