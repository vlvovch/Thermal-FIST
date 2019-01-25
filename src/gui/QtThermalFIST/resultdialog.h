/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef RESULTDIALOG_H
#define RESULTDIALOG_H

#include <QDialog>
#include <QTableView>
#include <QTableWidget>
#include <QPushButton>
#include <QTextEdit>

#include "BaseStructures.h"

#include "HRGBase/ThermalModelBase.h"
#include "tablemodel.h"


class ResultDialog : public QDialog
{
    Q_OBJECT


    thermalfist::ThermalModelBase *model;
    ChargesFluctuations *flucts;


    QTextEdit *parameters;
    QTextEdit *results;

    QString GetParameters();
    QString GetResults();
public:
    explicit  ResultDialog(QWidget *parent = 0, thermalfist::ThermalModelBase *mod = NULL, ChargesFluctuations *flucts_in = NULL);

signals:

public slots:
    void checkFixTableSize();
};


#endif // RESULTDIALOG_H
