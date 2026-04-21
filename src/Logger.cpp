// SPDX-License-Identifier: GPL-3.0-or-later
#include "Logger.hpp"
#include "version.hpp"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <sstream>

#ifdef _WIN32
  #include <io.h>
  #define ISATTY(fd) _isatty(fd)
  #define FILENO(f)  _fileno(f)
#else
  #include <unistd.h>
  #define ISATTY(fd) isatty(fd)
  #define FILENO(f)  fileno(f)
#endif

// ---- static member definitions ----
Logger::Level   Logger::s_level          = Logger::Level::INFO;
std::ofstream   Logger::s_file;
bool            Logger::s_fileOpen       = false;
bool            Logger::s_tty            = false;
bool            Logger::s_progressActive = false;
std::mutex      Logger::s_mutex;

// ---------------------------------------------------------------------------
void Logger::init(Level level, const std::string& logFilePath)
{
    std::lock_guard<std::mutex> lock(s_mutex);
    s_level = level;
    s_tty   = (ISATTY(FILENO(stdout)) != 0);

    if (!logFilePath.empty())
    {
        s_file.open(logFilePath, std::ios::out | std::ios::trunc);
        s_fileOpen = s_file.is_open();
    }
}

// ---------------------------------------------------------------------------
static std::string timestamp()
{
    auto now  = std::chrono::system_clock::now();
    auto tt   = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf{};
#ifdef _WIN32
    localtime_s(&tm_buf, &tt);
#else
    localtime_r(&tt, &tm_buf);
#endif
    std::ostringstream ss;
    ss << std::put_time(&tm_buf, "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

// ---------------------------------------------------------------------------
void Logger::banner()
{
    std::lock_guard<std::mutex> lock(s_mutex);

    const std::string line = "----------------------------------------------";
    const std::string b =
        "\n" + line + "\n"
        "OASIS " + gVERSION + " -- Offshore Advanced Simulation Software\n"
        "Copyright (C) 2026  IHCantabria\n"
        "Licensed under GNU GPL v3 -- see LICENSE file.\n"
        "This program comes with ABSOLUTELY NO WARRANTY.\n"
        + line + "\n";

    std::cout << b << std::flush;
    if (s_fileOpen)
        s_file << "[" << timestamp() << "] " << b;
}

// ---------------------------------------------------------------------------
void Logger::write(const std::string& msg, bool toStderr)
{
    // If the last thing printed was a \r progress line, close it first
    if (s_progressActive && !toStderr)
    {
        std::cout << "\n";
        s_progressActive = false;
    }

    if (toStderr)
        std::cerr << msg << "\n" << std::flush;
    else
        std::cout << msg << "\n" << std::flush;

    if (s_fileOpen)
        s_file << "[" << timestamp() << "] " << msg << "\n" << std::flush;
}

// ---------------------------------------------------------------------------
void Logger::debug(const std::string& msg)
{
    if (s_level < Level::DEBUG) return;
    std::lock_guard<std::mutex> lock(s_mutex);
    write("[DEBUG] " + msg);
}

// ---------------------------------------------------------------------------
void Logger::info(const std::string& msg)
{
    if (s_level < Level::INFO) return;
    std::lock_guard<std::mutex> lock(s_mutex);
    write(msg);
}

// ---------------------------------------------------------------------------
void Logger::fast(const std::string& msg)
{
    std::lock_guard<std::mutex> lock(s_mutex);
    write(msg);
}

// ---------------------------------------------------------------------------
void Logger::warning(const std::string& msg)
{
    if (s_level < Level::INFO) return;
    std::lock_guard<std::mutex> lock(s_mutex);
    write("[WARNING] " + msg);
}

// ---------------------------------------------------------------------------
void Logger::error(const std::string& msg)
{
    std::lock_guard<std::mutex> lock(s_mutex);
    // Close any active progress line on stdout first
    if (s_progressActive)
    {
        std::cout << "\n" << std::flush;
        s_progressActive = false;
    }
    std::cerr << "[ERROR] " << msg << "\n" << std::flush;
    if (s_fileOpen)
        s_file << "[" << timestamp() << "] [ERROR] " << msg << "\n" << std::flush;
}

// ---------------------------------------------------------------------------
void Logger::progress(const std::string& msg)
{
    std::lock_guard<std::mutex> lock(s_mutex);

    if (s_tty)
    {
        // Overwrite current line on interactive terminal
        std::cout << "\r" << msg << std::flush;
        s_progressActive = true;
    }
    else
    {
        // Piped/redirected: print normally
        std::cout << msg << "\n" << std::flush;
    }

    // Always write full line to log file
    if (s_fileOpen)
        s_file << "[" << timestamp() << "] " << msg << "\n" << std::flush;
}

// ---------------------------------------------------------------------------
void Logger::close()
{
    std::lock_guard<std::mutex> lock(s_mutex);
    if (s_progressActive)
    {
        std::cout << "\n" << std::flush;
        s_progressActive = false;
    }
    if (s_fileOpen)
    {
        s_file.close();
        s_fileOpen = false;
    }
}
