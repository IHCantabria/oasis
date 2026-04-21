// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LOGGER_HPP__
#define LOGGER_HPP__

#include <string>
#include <fstream>
#include <mutex>

// =============================================================================
// Logger — centralised output for OASIS
//
// Verbosity levels (least → most verbose):
//   FAST  — timestep progress + errors only
//   INFO  — FAST + section headers, lifecycle messages, summaries  [default]
//   DEBUG — INFO + all property dumps and internal diagnostics
//
// Usage (from any translation unit):
//   #include "Logger.hpp"
//   Logger::info("Reading body properties...");
//   Logger::debug("  mass = " + std::to_string(m));
//   Logger::progress("  t = 12.50 / 600.00 s");   // overwrites previous line
//   Logger::error("NaN detected in system matrix");
// =============================================================================

class Logger
{
public:
    enum class Level { FAST = 0, INFO = 1, DEBUG = 2 };

    // ------------------------------------------------------------------
    // Initialise the logger.  Call once from main() before any output.
    //   level       — verbosity level (default: INFO)
    //   logFilePath — if non-empty, also write to this file
    // ------------------------------------------------------------------
    static void init(Level level = Level::INFO, const std::string& logFilePath = "");

    // ------------------------------------------------------------------
    // GPL banner + version.  Always printed regardless of level.
    // ------------------------------------------------------------------
    static void banner();

    // ------------------------------------------------------------------
    // Logging methods
    // ------------------------------------------------------------------
    static void debug(const std::string& msg);    // printed at DEBUG only
    static void info(const std::string& msg);     // printed at INFO and DEBUG
    static void fast(const std::string& msg);     // always printed
    static void warning(const std::string& msg);  // printed at INFO and DEBUG, prefix [WARNING]
    static void error(const std::string& msg);    // always printed to stderr, prefix [ERROR]

    // ------------------------------------------------------------------
    // In-place progress line.
    // On an interactive terminal: uses \r to overwrite the current line.
    // When piped / redirected or writing to log file: uses \n.
    // ------------------------------------------------------------------
    static void progress(const std::string& msg);

    // ------------------------------------------------------------------
    // Close the log file (called at the end of main()).
    // ------------------------------------------------------------------
    static void close();

    // Prevent instantiation
    Logger() = delete;

private:
    static Level          s_level;
    static std::ofstream  s_file;
    static bool           s_fileOpen;
    static bool           s_tty;         // true when stdout is an interactive terminal
    static bool           s_progressActive; // true after a \r progress line was written
    static std::mutex     s_mutex;

    // Write to terminal and/or log file
    static void write(const std::string& msg, bool toStderr = false);
    // Write a log-file-only line (with newline, no \r)
    static void writeFile(const std::string& msg);
};

#endif // LOGGER_HPP__
