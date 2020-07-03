/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2019 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef HELPERROUTINES_H
#define HELPERROUTINES_H

#include <QTableWidget>

class QTableWidgetCC : public QTableWidget
{
  Q_OBJECT
public:
  explicit  QTableWidgetCC(QWidget *parent = 0)
    : QTableWidget(parent)
  { }

  public slots :
    void keyPressEvent(QKeyEvent *event);
};

void copyTableViewSelectionToClipBoard(QTableView *view);

#endif