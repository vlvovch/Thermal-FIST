/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2023 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "trajectoriestab.h"
#include <QHBoxLayout>

TrajectoriesTab::TrajectoriesTab(QWidget *parent, const QVector<QWidget*>& tabs, const QVector<QString>& tabnames) :
    QWidget(parent)
{
  tabWidget = new QTabWidget(this);
  for(int i = 0; i < tabs.size(); i++)
    tabWidget->addTab(tabs[i], tabnames[i]);

  QHBoxLayout *layout = new QHBoxLayout(this);
  layout->addWidget(tabWidget);
  setLayout(layout);
}
