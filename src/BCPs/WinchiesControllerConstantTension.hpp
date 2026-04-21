// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef WINCHIE_CONTROLLER_CONSTANT_TENSION_FLAG
#define WINCHIE_CONTROLLER_CONSTANT_TENSION_FLAG
#include "WinchiesController.hpp"

class WinchieControllerConstantTension : public WinchieController
{
public:
    double targetTension; // Target tension per winch [N]

    WinchieControllerConstantTension(int n, Winchie** Ws, Simulation* pIncSim);
    ~WinchieControllerConstantTension() override;

    void ReadPropertiesASCII(FILE* pFile) override;
    void ReadPropertiesYAML(YAML::Node node) override;
    void SetUpWinchiesController(void) override;
    void controlWinchies(double time) override;
};

#endif
