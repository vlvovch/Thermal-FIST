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

// Resize canvas and Qt window to fill browser
static void resizeToFillBrowser()
{
    // Get browser window size and device pixel ratio
    double dpr = EM_ASM_DOUBLE({ return window.devicePixelRatio || 1; });
    int cssWidth = EM_ASM_INT({ return window.innerWidth; });
    int cssHeight = EM_ASM_INT({ return window.innerHeight; });

    // Canvas internal resolution (accounting for device pixel ratio)
    int canvasWidth = static_cast<int>(cssWidth * dpr);
    int canvasHeight = static_cast<int>(cssHeight * dpr);

    // Set canvas internal size (resolution) and CSS display size
    EM_ASM({
        var canvas = document.querySelector('canvas');
        if (canvas) {
            // Set internal resolution
            canvas.width = $0;
            canvas.height = $1;
            // Set CSS display size
            canvas.style.width = $2 + 'px';
            canvas.style.height = $3 + 'px';
        }
        // Prevent scrollbars
        document.body.style.margin = '0';
        document.body.style.overflow = 'hidden';
        document.documentElement.style.overflow = 'hidden';
    }, canvasWidth, canvasHeight, cssWidth, cssHeight);

    // Resize Qt window to match CSS size and update layout
    if (g_mainWindow) {
        g_mainWindow->setGeometry(0, 0, cssWidth, cssHeight);
        g_mainWindow->updateGeometry();
    }
}

// Callback for browser resize events
EM_BOOL onBrowserResize(int eventType, const EmscriptenUiEvent *uiEvent, void *userData)
{
    Q_UNUSED(eventType);
    Q_UNUSED(uiEvent);
    Q_UNUSED(userData);
    resizeToFillBrowser();
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

    // Register resize callback for browser window resize events
    emscripten_set_resize_callback(EMSCRIPTEN_EVENT_TARGET_WINDOW, nullptr, EM_TRUE, onBrowserResize);

    // Resize canvas and show window (use show() instead of showMaximized() for better WASM sizing)
    resizeToFillBrowser();
    w.show();
#else
    w.showMaximized();
#endif

    splash.finish(&w);

    return a.exec();
}
