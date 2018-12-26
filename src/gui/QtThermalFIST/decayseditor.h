/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef DECAYSEDITOR_H
#define DECAYSEDITOR_H

#include <QDialog>
#include <QTableView>
#include <QPushButton>
#include "listtablemodel.h"

#include "HRGBase/ThermalParticleSystem.h"

class decayseditor : public QDialog
{
    Q_OBJECT

    QTableView *tableDecays;
    DecayEditorTableModel *myModel;

    QPushButton *buttonAddDecay;
    QPushButton *buttonRemoveDecay;

    QPushButton *buttonAddDaughterColumn;
    QPushButton *buttonRemoveDaughterColumn;
public:
    explicit decayseditor(QWidget *parent = 0, thermalfist::ThermalParticle *Particle = NULL, thermalfist::ThermalParticleSystem *TPS = NULL);

signals:

public slots:
    void checkFixTableSize();
    void addDecay();
    void removeDecay();
    void addColumn();
    void removeColumn();
};

#endif // DECAYSEDITOR_H
