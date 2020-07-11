/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef ABOUTDIALOG_H
#define ABOUTDIALOG_H

#include <QDialog>
#include <QPushButton>
#include <QTextEdit>

class AboutDialog : public QDialog
{
  Q_OBJECT
public:
  explicit  AboutDialog(QWidget *parent = 0);
protected:
  void showEvent(QShowEvent* ev);
};

#endif // RESULTDIALOG_H
