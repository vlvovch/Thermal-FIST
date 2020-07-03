/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef PARTICLEDIALOG_H
#define PARTICLEDIALOG_H

#include <QDialog>
#include <QTableView>
#include <QTableWidget>
#include <QPushButton>
#include <QTextEdit>
#include "decaytablemodel.h"

#include "HRGBase/ThermalModelBase.h"

class ParticleDialog : public QDialog
{
    Q_OBJECT

    thermalfist::ThermalModelBase *model;
    int pid;


    QTextEdit *data;

    QTableView *tableDecays;
    QTableWidget *tableSources;
    DecayTableModel *myModel;

    QPushButton *buttonAddDecay;
    QPushButton *buttonRemoveDecay;

    QPushButton *buttonAddDaughterColumn;
    QPushButton *buttonRemoveDaughterColumn;

    QPushButton *buttonSpectralFunction;

    QString GetParticleInfo();

public:
    explicit  ParticleDialog(QWidget *parent = 0, thermalfist::ThermalModelBase *mod = NULL, int ParticleID = -1);

signals:

public slots:
    void checkFixTableSize();
    void addDecay();
    void removeDecay();
    void addColumn();
    void removeColumn();
    void showSpectralFunction();
    void contextMenuRequest(QPoint pos);
    void copyFeeddownTable();

protected:
  bool eventFilter(QObject* obj, QEvent* ev);
};

#endif // DECAYSEDITOR_H
