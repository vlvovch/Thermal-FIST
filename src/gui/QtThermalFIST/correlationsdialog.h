/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef CORRELATIONSDIALOG_H
#define CORRELATIONSDIALOG_H

#include <QDialog>
#include <QTableView>
#include <QTableWidget>
#include <QPushButton>
#include <QTextEdit>
#include <QComboBox>

#include "HRGBase/ThermalModelBase.h"

class CorrelationsDialog : public QDialog
{
    Q_OBJECT

    thermalfist::ThermalModelBase *model;
    
    QTableWidget *tableCorr;
    
    QComboBox *comboFeeddown;
    QComboBox *comboType;
    QComboBox *comboQuantity;

public:
    explicit  CorrelationsDialog(QWidget *parent = 0, thermalfist::ThermalModelBase *mod = NULL);

signals:

public slots:
    void checkFixTableSize();
    void recalculate();

protected:
    bool eventFilter(QObject* obj, QEvent* ev);
};

#endif // DECAYSEDITOR_H
