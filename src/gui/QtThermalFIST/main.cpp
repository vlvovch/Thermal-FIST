/*
 * Thermal-FIST package
 * 
 * Copyright (c) 2014-2018 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#include "mainwindow.h"
#include <QApplication>
#include <QTranslator>
#include <QLibraryInfo>

int main(int argc, char *argv[])
{

    QApplication a(argc, argv);

    QTranslator qtTranslator;
    qtTranslator.load("qt_" + QLocale::system().name(),
            QLibraryInfo::location(QLibraryInfo::TranslationsPath));
    a.installTranslator(&qtTranslator);

    QTranslator myappTranslator;
    myappTranslator.load("HadronResonanceGas_" + QLocale::system().name());
    a.installTranslator(&myappTranslator);

    QFont font = QApplication::font();
    font.setPointSize(10);
    QApplication::setFont(font);

    //std::ios_base::sync_with_stdio(false);

    MainWindow w;
    w.showMaximized();

    return a.exec();
}
