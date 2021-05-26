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
  fontTitle.setPointSize(QApplication::font().pointSize() + 4);
  fontTitle.setBold(true);

  QFont fontDefault = QApplication::font();
  //fontDefault.setPointSize(10);

  QString title = "Thermal-FIST";
  QString titleVersion = "Version " + QString::number(ThermalFIST_VERSION_MAJOR) + "." + QString::number(ThermalFIST_VERSION_MINOR);
  if (ThermalFIST_VERSION_DEVEL != 0) titleVersion += "." + QString::number(ThermalFIST_VERSION_DEVEL);

  QLabel *labTitle = new QLabel(title);
  labTitle->setFont(fontTitle);

  layout->addWidget(labellogo, 0, Qt::AlignCenter);
  layout->addWidget(labTitle, 0, Qt::AlignCenter);
  QLabel *labelVersion = new QLabel(titleVersion);
  labelVersion->setFont(fontDefault);
  layout->addWidget(labelVersion, 0, Qt::AlignCenter);
  layout->addSpacing(20);
  QLabel *labelCC = new QLabel(tr("Copyright (c) 2014-2021 Volodymyr Vovchenko"));
  labelCC->setFont(fontDefault);
  layout->addWidget(labelCC, 0, Qt::AlignCenter);
  QLabel *labelLic = new QLabel(tr("GNU General Public License (GPLv3 or later)"));
  labelLic->setFont(fontDefault);
  layout->addWidget(labelLic, 0, Qt::AlignCenter);

  QFont fontSmaller = QApplication::font();
  fontSmaller.setPointSize(QApplication::font().pointSize() - 1);

  QLabel *labelAttrib = new QLabel(tr("Publication of results obtained using this code should include a reference to:"));
  QLabel *labelAttrib2 = new QLabel(tr("V. Vovchenko, H. Stoecker, <a href=\"https://dx.doi.org/10.1016/j.cpc.2019.06.024\">Comput. Phys. Commun. 244, 295 (2019)</a> [<a href=\"https://arxiv.org/abs/1901.05249\">arXiv:1901.05249</a>]"));
  labelAttrib2->setTextFormat(Qt::RichText);
  labelAttrib2->setTextInteractionFlags(Qt::TextBrowserInteraction);
  labelAttrib2->setOpenExternalLinks(true);
  labelAttrib->setFont(fontSmaller);
  labelAttrib2->setFont(fontSmaller);

  layout->addSpacing(20);
  layout->addWidget(labelAttrib);
  layout->addWidget(labelAttrib2);
  

  QLabel *labelLatest = new QLabel(tr("The latest version of the program is available at <a href=\"https://github.com/vlvovch/Thermal-FIST\">https://github.com/vlvovch/Thermal-FIST</a>"));
  labelLatest->setTextFormat(Qt::RichText);
  labelLatest->setTextInteractionFlags(Qt::TextBrowserInteraction);
  labelLatest->setOpenExternalLinks(true);
  labelLatest->setFont(fontSmaller);
  

  layout->addSpacing(10);
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

  QLabel *labelSupportedBy= new QLabel(tr("Supported by"));
  labelSupportedBy->setFont(fontSmaller);

  QHBoxLayout *laySupported = new QHBoxLayout;
  laySupported->setAlignment(Qt::AlignCenter);
  QPixmap logoFIAS(":/images/fias.png");
  QLabel *labellogoFIAS = new QLabel;
  labellogoFIAS->setPixmap(logoFIAS.scaledToHeight(60, Qt::SmoothTransformation));

  QPixmap logoGoethe(":/images/Goethe-Uni.png");
  QLabel *labellogoGoethe = new QLabel;
  labellogoGoethe->setPixmap(logoGoethe.scaledToHeight(70, Qt::SmoothTransformation));

  QPixmap logoKNU(":/images/knu.png");
  QLabel *labellogoKnu = new QLabel;
  labellogoKnu->setPixmap(logoKNU.scaledToHeight(80, Qt::SmoothTransformation));

  laySupported->addWidget(labellogoGoethe);
  laySupported->addSpacing(70);
  laySupported->addWidget(labellogoFIAS);
  laySupported->addSpacing(70);
  laySupported->addWidget(labellogoKnu);

  layout->addSpacing(20);
  //layout->addWidget(labelSupportedBy, 0, Qt::AlignLeft);
  layout->addLayout(laySupported, 1);

  setLayout(layout);

  this->setWindowTitle(tr("About"));
}