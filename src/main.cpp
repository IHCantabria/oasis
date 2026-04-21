// SPDX-License-Identifier: GPL-3.0-or-later

#include <string>
#include <filesystem>

#include "Logger.hpp"
#include "Simulations/Simulation.hpp"
#include "Exceptions/Exception.hpp"

int main(int argc, char* argv[])
{
    // ------------------------------------------------------------------
    // Parse command-line arguments
    //   argv[1]              project root path (required)
    //   --verbosity LEVEL    DEBUG | INFO | FAST  (default: INFO)
    //   --log <file>         also write output to this file
    // ------------------------------------------------------------------
    std::string project_path;
    Logger::Level verbosity = Logger::Level::INFO;
    std::string log_file;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--verbosity" && i + 1 < argc)
        {
            std::string v = argv[++i];
            if (v == "DEBUG")      verbosity = Logger::Level::DEBUG;
            else if (v == "FAST")  verbosity = Logger::Level::FAST;
            else                   verbosity = Logger::Level::INFO;
        }
        else if (arg == "--log" && i + 1 < argc)
        {
            log_file = argv[++i];
        }
        else if (project_path.empty() && arg[0] != '-')
        {
            project_path = arg;
        }
    }

    Logger::init(verbosity, log_file);
    Logger::banner();

    if (project_path.empty())
    {
        Logger::error("No project path provided. Usage: oasis <project_path> [--verbosity DEBUG|INFO|FAST] [--log <file>]");
        return 1;
    }

    Logger::info("Simulation directory: " + (project_path.empty() ? "./  (pwd)" : project_path));

    try
    {
        // Auto-detect input format
        std::string data_format;
        std::filesystem::path yaml_path = std::filesystem::path(project_path) / "input" / "dataProblem.yaml";
        std::filesystem::path dat_path  = std::filesystem::path(project_path) / "input" / "dataProblem.dat";

        if (std::filesystem::exists(yaml_path))
        {
            data_format = "YAML";
            Logger::info("Input format: YAML (dataProblem.yaml)");
        }
        else if (std::filesystem::exists(dat_path))
        {
            data_format = "ASCII";
            Logger::info("Input format: ASCII (dataProblem.dat)");
        }
        else
        {
            Logger::error("No input file found. Expected input/dataProblem.yaml or input/dataProblem.dat");
            return 1;
        }

        Logger::info("Creating simulation...");
        Simulation* mySim = new Simulation(project_path, data_format);

        Logger::info("Loading...");
        mySim->LoadCase();

        Logger::info("Initializing...");
        mySim->Initialize();

        Logger::info("Running...");
        mySim->Run();

        Logger::info("Closing simulation...");
        mySim->CloseCase();

        delete mySim;
    }
    catch (Exception& error)
    {
        error.PrintDebug();
    }

    Logger::info("End of OASIS.");
    Logger::close();

    return 0;
}