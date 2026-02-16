/*
 * Thermal-FIST package
 *
 * Copyright (c) 2014-2026 Volodymyr Vovchenko
 *
 * GNU General Public License (GPLv3 or later)
 */
#ifndef WASMFILEIO_H
#define WASMFILEIO_H

#include <QFile>
#include <QDir>
#include <QFileInfo>
#include <QByteArray>
#include <QString>
#include <QFileDialog>
#include <QBuffer>
#include <QTextStream>
#include <functional>

#ifdef Q_OS_WASM
#include <emscripten.h>
#include <emscripten/threading.h>
#endif

/**
 * @brief Helper functions for WebAssembly file I/O
 *
 * In WebAssembly, browser file pickers don't provide real filesystem paths.
 * Instead, we use Qt's getOpenFileContent/saveFileContent which work with
 * file contents directly. This helper provides utilities to:
 * 1. Write uploaded file contents to a sandbox filesystem path
 * 2. Trigger browser downloads for exported data
 */
namespace WasmFileIO {

/**
 * @brief Write uploaded file bytes to a temporary sandbox path
 *
 * Creates a file in /tmp with the given name and contents.
 * The returned path can be used with existing file-path-based APIs.
 *
 * @param originalName The original filename (used for naming the temp file)
 * @param bytes The file contents
 * @return The sandbox filesystem path, or empty string on failure
 */
inline QString writeUploadedFileToSandbox(const QString& originalName,
                                          const QByteArray& bytes)
{
    // /tmp exists in Qt/WASM builds typically; if not, create it.
    QDir dir("/tmp");
    if (!dir.exists())
        dir.mkpath(".");

    QString base = QFileInfo(originalName).fileName();
    if (base.isEmpty())
        base = "upload.dat";

    QString path = dir.filePath(base);

    QFile f(path);
    if (!f.open(QIODevice::WriteOnly))
        return QString();

    f.write(bytes);
    f.close();

    return path;
}

/**
 * @brief Get the sandbox directory path for temporary files
 * @return The path to the sandbox temp directory
 */
inline QString getSandboxTempDir()
{
    return QStringLiteral("/tmp");
}

/**
 * @brief Check if a file exists in the sandbox filesystem
 * @param path The path to check
 * @return true if file exists
 */
inline bool sandboxFileExists(const QString& path)
{
    return QFileInfo::exists(path);
}

#ifdef Q_OS_WASM

/**
 * @brief Open a file using browser file picker (WASM only)
 *
 * Uses QFileDialog::getOpenFileContent for async file selection.
 * The callback receives the sandbox path where the file was written.
 *
 * @param parent Parent widget for the dialog
 * @param caption Dialog caption
 * @param filter File filter (e.g., "*.dat *.txt")
 * @param callback Function called with (sandboxPath) on success, empty string on cancel
 */
inline void openFile(QWidget* parent,
                     const QString& caption,
                     const QString& filter,
                     std::function<void(const QString& sandboxPath)> callback)
{
    QFileDialog::getOpenFileContent(
        filter,
        [callback](const QString& fileName, const QByteArray& fileContent) {
            if (fileName.isEmpty()) {
                callback(QString());
                return;
            }

            QString sandboxPath = writeUploadedFileToSandbox(fileName, fileContent);
            callback(sandboxPath);
        }
    );
    Q_UNUSED(parent);
    Q_UNUSED(caption);
}

/**
 * @brief Save data to a file via browser download (WASM only)
 *
 * Triggers a browser download with the given data and filename.
 *
 * @param data The data to save
 * @param suggestedName Suggested filename for the download
 */
inline void saveFile(const QByteArray& data, const QString& suggestedName)
{
    QFileDialog::saveFileContent(data, suggestedName);
}

/**
 * @brief Save text to a file via browser download (WASM only)
 *
 * Convenience overload for text data.
 *
 * @param text The text to save
 * @param suggestedName Suggested filename for the download
 */
inline void saveTextFile(const QString& text, const QString& suggestedName)
{
    saveFile(text.toUtf8(), suggestedName);
}

#endif // Q_OS_WASM

/**
 * @brief Check if threading is supported in the current WASM environment
 *
 * WASM threading requires:
 * 1. Build with -pthread flag
 * 2. Browser with SharedArrayBuffer support
 * 3. Server serving COOP/COEP headers
 *
 * Note: Safari has unreliable threading support even when it claims availability,
 * so we force synchronous execution on Safari.
 *
 * @return true if threading is available, false otherwise
 */
inline bool isThreadingAvailable()
{
#ifdef Q_OS_WASM
    // Check if threading support is enabled at runtime
    if (!emscripten_has_threading_support())
        return false;

    // Safari has unreliable WASM threading - force synchronous execution
    // Check user agent for Safari (but not Chrome which also contains "Safari")
    bool isSafari = EM_ASM_INT({
        var ua = navigator.userAgent;
        return (ua.indexOf('Safari') !== -1 && ua.indexOf('Chrome') === -1) ? 1 : 0;
    });

    if (isSafari)
        return false;

    return true;
#else
    // Native builds always have threading
    return true;
#endif
}

} // namespace WasmFileIO

#endif // WASMFILEIO_H
