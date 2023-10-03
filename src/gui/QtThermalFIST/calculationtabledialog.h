/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef CALCULATIONTABLEDIALOG_H
#define CALCULATIONTABLEDIALOG_H

#include <QDialog>
#include <QTableView>
#include <QTableWidget>
#include <QPushButton>
#include <QTextEdit>
#include <QComboBox>
#include <QLabel>

#include "HRGBase/ThermalModelBase.h"

struct CalculationTable {
  QString parameter_name;
  QVector<double> parameter_values;
  QVector<double> temperature_values;
  QVector<QString> quantities_names;
  QVector<std::vector<double>> quantities_values;
  QVector<QString> densities_names;
  QVector<std::vector<std::vector<double>>> densities_values;
  CalculationTable() { clear(); }
  void clear();
};

class CalculationTableDialog : public QDialog
{
    Q_OBJECT

    CalculationTable m_table;

    QComboBox *comboQuantity;
    QLabel *labelFeeddown;
    QComboBox *comboFeeddown;
    QTableWidget *tableValues;
    QPushButton *buttonWriteToFile;

public:
    explicit CalculationTableDialog(QWidget *parent = 0, const CalculationTable& table = CalculationTable());

signals:

public slots:
    void checkFixTableSize();
    void fillTable();
    void contextMenuRequest(QPoint pos);
    void copyTable();
    void writeTableToFile();

protected:
    bool eventFilter(QObject* obj, QEvent* ev);
};

#endif
