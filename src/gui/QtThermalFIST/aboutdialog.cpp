/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "aboutdialog.h"

#include <QLayout>
#include <QLabel>
#include <QApplication>
#include <QDebug>

#include "ThermalFISTConfig.h"


AboutDialog::AboutDialog(QWidget *parent) :
  QDialog(parent)
{
  QVBoxLayout *layout = new QVBoxLayout();
  layout->setAlignment(Qt::AlignCenter);

  QPixmap logo(":/images/FIST-about.png");
  QLabel *labellogo = new QLabel;
  labellogo->setPixmap(logo);

  


  QFont fontTitle = QApplication::font();
  fontTitle.setPointSize(14);
  fontTitle.setBold(true);

  QString title = "Thermal-FIST";
  QString titleVersion = "Version " + QString::number(ThermalFIST_VERSION_MAJOR) + "." + QString::number(ThermalFIST_VERSION_MINOR);
  if (ThermalFIST_VERSION_DEVEL != 0) titleVersion += "." + QString::number(ThermalFIST_VERSION_DEVEL);

  QLabel *labTitle = new QLabel(title);
  labTitle->setFont(fontTitle);

  layout->addWidget(labellogo, 0, Qt::AlignCenter);
  layout->addWidget(labTitle, 0, Qt::AlignCenter);
  layout->addWidget(new QLabel(titleVersion), 0, Qt::AlignCenter);
  layout->addSpacing(20);
  layout->addWidget(new QLabel(tr("Copyright (c) 2014-2018 Volodymyr Vovchenko")), 0, Qt::AlignCenter);
  layout->addWidget(new QLabel(tr("GNU General Public License (GPLv3 or later)")), 0, Qt::AlignCenter);
  

  QLabel *labelLatest = new QLabel(tr("The latest version of the program is available at <a href=\"https://github.com/vlvovch/Thermal-FIST\">https://github.com/vlvovch/Thermal-FIST</a>"));
  labelLatest->setTextFormat(Qt::RichText);
  labelLatest->setTextInteractionFlags(Qt::TextBrowserInteraction);
  labelLatest->setOpenExternalLinks(true);
  QFont fontSmaller = QApplication::font();
  fontSmaller.setPointSize(9);
  labelLatest->setFont(fontSmaller);

  layout->addSpacing(20);
  layout->addWidget(labelLatest);

  // contact by e-mail
  // some obfuscation (just in case)
  QString email = "";
  email += 'v';
  email += QChar(email[0].toLatin1() - 7);
  email += email[0];
  email += QChar('a' + 2);
  email += QChar(email[email.size() - 1].toLatin1() + 5);
  email += "enk";
  email += "@";
  email += "o";
  QChar ch1 = email[email.size() - 2], ch2 = email[email.size() - 1];
  email[email.size() - 1] = ch1;
  email[email.size() - 2] = ch2;
  email += "fias";
  email += ".";
  email += "uni-frankfurt.de";
  
  QLabel *labelEmail = new QLabel("For questions, suggestions or bug reports please contact <a href='mailto:" + email + "?subject=About Thermal-FIST'>" + email + "</a>");
  labelEmail->setTextFormat(Qt::RichText);
  labelEmail->setTextInteractionFlags(Qt::TextBrowserInteraction);
  labelEmail->setOpenExternalLinks(true);
  labelEmail->setFont(fontSmaller);

  layout->addWidget(labelEmail);

  setLayout(layout);

  this->setWindowTitle(tr("About"));
}