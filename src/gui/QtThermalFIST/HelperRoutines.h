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
#include <QTableView>
#include <QHeaderView>
#include <QApplication>
#include <QFontMetrics>

// Configure table view row heights to match current font size
// Call this after creating a QTableView to ensure proper row heights
inline void configureTableRowHeight(QTableView* table) {
    if (!table) return;
    QFontMetrics fm(QApplication::font());
    // Row height = font height + comfortable padding
    int rowHeight = fm.height() + 10;
    table->verticalHeader()->setDefaultSectionSize(rowHeight);
}

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