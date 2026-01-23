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
#include <QtGlobal>
#include <QScreen>

#ifdef Q_OS_WASM
#include <emscripten.h>
#include <emscripten/html5.h>
#endif

//#define CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#ifdef Q_OS_WASM
// Global pointer for resize callback
static MainWindow* g_mainWindow = nullptr;

// Callback for browser resize events
EM_BOOL onBrowserResize(int eventType, const EmscriptenUiEvent *uiEvent, void *userData)
{
    Q_UNUSED(eventType);
    Q_UNUSED(uiEvent);
    Q_UNUSED(userData);

    if (g_mainWindow) {
        // Get current browser window size
        int width = EM_ASM_INT({ return window.innerWidth; });
        int height = EM_ASM_INT({ return window.innerHeight; });

        // Resize the main window to fill the browser
        g_mainWindow->resize(width, height);
    }
    return EM_TRUE;
}
#endif

int main(int argc, char *argv[])
{
    //_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

    QApplication a(argc, argv);

    QFont font = QApplication::font();
    QFont fontSmall = font;
    font.setPointSize(10);
    QApplication::setFont(font);

    QPixmap pixmap(":/images/FIST.png");
    QSplashScreen splash(pixmap);
    splash.show();
    splash.setFont(fontSmall);
    splash.showMessage(QObject::tr("Initializing QtThermalFIST..."),
      Qt::AlignLeft | Qt::AlignTop, Qt::black);  //This line represents the alignment of text, color and position

    QTranslator qtTranslator;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    qtTranslator.load("qt_" + QLocale::system().name(),
            QLibraryInfo::location(QLibraryInfo::TranslationsPath));
#else
    qtTranslator.load("qtbase_" + QLocale::system().name(),
            QLibraryInfo::path(QLibraryInfo::TranslationsPath));
#endif
    a.installTranslator(&qtTranslator);

    QTranslator myappTranslator;
    myappTranslator.load("HadronResonanceGas_" + QLocale::system().name());
    a.installTranslator(&myappTranslator);



    //std::ios_base::sync_with_stdio(false);


    MainWindow w;

#ifdef Q_OS_WASM
    // WASM: Set up automatic window resizing to fill browser
    g_mainWindow = &w;

    // Register resize callback
    emscripten_set_resize_callback(EMSCRIPTEN_EVENT_TARGET_WINDOW, nullptr, EM_TRUE, onBrowserResize);

    // Set initial size to fill browser window
    int width = EM_ASM_INT({ return window.innerWidth; });
    int height = EM_ASM_INT({ return window.innerHeight; });
    w.resize(width, height);
    w.show();
#else
    w.showMaximized();
#endif

    splash.finish(&w);

    return a.exec();
}
