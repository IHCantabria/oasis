// SPDX-License-Identifier: GPL-3.0-or-later
#include <armadillo>
#include <string>
#include <math.h>
#include "WinchiesControllerHorizontal.hpp"
#include "../Simulations/Simulation.hpp"
#include "../CommonTools.hpp"
#include "../os_tools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../Logger.hpp"
#include "../MathTools.hpp"

WinchieControllerHorizontal::WinchieControllerHorizontal(int n, Winchie** Ws, Simulation* pIncSim)
    : WinchieController(n, Ws, pIncSim)
{
    pfile_FC = NULL;
    pfile_RP = NULL;
}

WinchieControllerHorizontal::~WinchieControllerHorizontal()
{
}

void WinchieControllerHorizontal::ReadPropertiesASCII(FILE* pFile)
{

    // Declare variables
    char buffer_line[1000];
    double dtemp;
    int itemp, itemp2;

    // NOTE: The 3 header lines and controller type have already been read by the factory.
    // We start reading directly from the body index.

    fscanf(pFile, "%d %[^\n]\n", &indBody, buffer_line);
    indBody = indBody - 1;

    // Read all parameters
    if (fscanf(pFile, "%lf", &Ac) != 1)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the controller state-space model parameters \n";
        throw ValueError(ss.str());
    }
    if (fscanf(pFile, "%lf", &Bc) != 1)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the controller state-space model parameters \n";
        throw ValueError(ss.str());
    }
    if (fscanf(pFile, "%lf", &Cc) != 1)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the controller state-space model parameters \n";
        throw ValueError(ss.str());
    }
    if (fscanf(pFile, "%lf", &Dc) != 1)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the controller state-space model parameters \n";
        throw ValueError(ss.str());
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    for (int ii = 0; ii < 3; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the lead controller gains \n";
            throw ValueError(ss.str());
        }
        Kc(ii, 0) = dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    for (int ii = 0; ii < 3; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the integral time gains \n";
            throw ValueError(ss.str());
        }
        Ki(ii, 0) = dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    if (fscanf(pFile, "%lf", &Af) != 1)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the filter state-space model parameters \n";
        throw ValueError(ss.str());
    }
    if (fscanf(pFile, "%lf", &Bf) != 1)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the filter state-space model parameters \n";
        throw ValueError(ss.str());
    }
    if (fscanf(pFile, "%lf", &Cf) != 1)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the filter state-space model parameters \n";
        throw ValueError(ss.str());
    }
    if (fscanf(pFile, "%lf", &Df) != 1)
    {
        std::stringstream ss;
        ss << "An error occurred when trying to read the filter state-space model parameters \n";
        throw ValueError(ss.str());
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    fscanf(pFile, "%lf %[^\n]\n", &Kw, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &time_ini, buffer_line);
    fscanf(pFile, "%ld %[^\n]\n", &reference_flag, buffer_line);

    if (reference_flag < 1)
    {
        // Skip the "State space reference position" header line
        fgets(buffer_line, sizeof(buffer_line), pFile);

        if (fscanf(pFile, "%lf", &Ar) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the reference state-space model parameters \n";
            throw ValueError(ss.str());
        }
        if (fscanf(pFile, "%lf", &Br) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the reference state-space model parameters \n";
            throw ValueError(ss.str());
        }
        if (fscanf(pFile, "%lf", &Cr) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the reference state-space model parameters \n";
            throw ValueError(ss.str());
        }
        if (fscanf(pFile, "%lf", &Dr) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the reference state-space model parameters \n";
            throw ValueError(ss.str());
        }
        fscanf(pFile, "%[^\n]\n", buffer_line);

        for (int ii = 0; ii < 3; ii++)
        {
            if (fscanf(pFile, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read the vector of reference positions \n";
                throw ValueError(ss.str());
            }
            ur(ii, 0) = dtemp;
        }
        fscanf(pFile, "%[^\n]\n", buffer_line);

        // Skip the "Discrete points reference position" header lines
        for (int ii = 0; ii < 5; ii++)
        {
            fgets(buffer_line, sizeof(buffer_line), pFile);
        }
    }
    else
    {
        // Skip the "State space reference position" header lines
        for (int ii = 0; ii < 3; ii++)
        {
            fgets(buffer_line, sizeof(buffer_line), pFile);
        }

        // Skip the "Discrete points reference position" header line
        fgets(buffer_line, sizeof(buffer_line), pFile);

        // Initialise reference signal data arrays
        t_ref = arma::zeros(reference_flag, 1);
        x_ref = arma::zeros(reference_flag, 1);
        y_ref = arma::zeros(reference_flag, 1);
        yaw_ref = arma::zeros(reference_flag, 1);

        for (int ii = 0; ii < reference_flag; ii++)
        {
            if (fscanf(pFile, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read the times for reference positions \n";
                throw ValueError(ss.str());
            }
            t_ref(ii, 0) = dtemp;
        }
        fscanf(pFile, "%[^\n]\n", buffer_line);
        for (int ii = 0; ii < reference_flag; ii++)
        {
            if (fscanf(pFile, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read the x for reference positions \n";
                throw ValueError(ss.str());
            }
            x_ref(ii, 0) = dtemp;
        }
        fscanf(pFile, "%[^\n]\n", buffer_line);
        for (int ii = 0; ii < reference_flag; ii++)
        {
            if (fscanf(pFile, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read the y for reference positions \n";
                throw ValueError(ss.str());
            }
            y_ref(ii, 0) = dtemp;
        }
        fscanf(pFile, "%[^\n]\n", buffer_line);
        for (int ii = 0; ii < reference_flag; ii++)
        {
            if (fscanf(pFile, "%lf", &dtemp) != 1)
            {
                std::stringstream ss;
                ss << "An error occurred when trying to read the yaw for reference positions \n";
                throw ValueError(ss.str());
            }
            yaw_ref(ii, 0) = dtemp;
        }
        fscanf(pFile, "%[^\n]\n", buffer_line);
    }

    // Skip the three "Inversor block inputs" header lines
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(buffer_line, sizeof(buffer_line), pFile);
    }

    fscanf(pFile, "%d %[^\n]\n", &inversor_flag, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &T_max, buffer_line, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &T_min, buffer_line, buffer_line);

    // Skip the "Straight lines inversor inputs" header line
    fgets(buffer_line, sizeof(buffer_line), pFile);

    fscanf(pFile, "%d %[^\n]\n", &nIterMax, buffer_line);
    fscanf(pFile, "%lf %[^\n]\n", &atol, buffer_line);

    // Skip the "Coefficients inversor inputs" header line
    fgets(buffer_line, sizeof(buffer_line), pFile);

    fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
    ind_x_pos = arma::zeros<arma::uvec>(itemp, 1);
    coef_x_pos = arma::zeros(itemp, 1);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%d", &itemp2) != 1)
        {
            Logger::debug("itemp2 = " + std::to_string(itemp2));
            std::stringstream ss;
            ss << "An error occurred when trying to read the indices of lines for positive X force \n";
            throw ValueError(ss.str());
        }
        ind_x_pos(ii) = itemp2 - 1;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the coefficients of lines for positive X force \n";
            throw ValueError(ss.str());
        }
        coef_x_pos(ii) = dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
    ind_x_neg = arma::zeros<arma::uvec>(itemp, 1);
    coef_x_neg = arma::zeros(itemp, 1);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%d", &itemp2) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the indices of lines for negative X force \n";
            throw ValueError(ss.str());
        }
        ind_x_neg(ii) = itemp2 - 1;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the coefficients of lines for negative X force \n";
            throw ValueError(ss.str());
        }
        coef_x_neg(ii) = dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
    ind_y_pos = arma::zeros<arma::uvec>(itemp, 1);
    coef_y_pos = arma::zeros(itemp, 1);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%d", &itemp2) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the indices of lines for positive Y force \n";
            throw ValueError(ss.str());
        }
        ind_y_pos(ii) = itemp2 - 1;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the coefficients of lines for positive Y force \n";
            throw ValueError(ss.str());
        }
        coef_y_pos(ii) = dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
    ind_y_neg = arma::zeros<arma::uvec>(itemp, 1);
    coef_y_neg = arma::zeros(itemp, 1);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%d", &itemp2) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the indices of lines for negative Y force \n";
            throw ValueError(ss.str());
        }
        ind_y_neg(ii) = itemp2 - 1;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the coefficients of lines for negative Y force \n";
            throw ValueError(ss.str());
        }
        coef_y_neg(ii) = dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
    ind_g_pos = arma::zeros<arma::uvec>(itemp, 1);
    coef_g_pos = arma::zeros(itemp, 1);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%d", &itemp2) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the indices of lines for positive YAW force \n";
            throw ValueError(ss.str());
        }
        ind_g_pos(ii) = itemp2 - 1;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the coefficients of lines for positive YAW force \n";
            throw ValueError(ss.str());
        }
        coef_g_pos(ii) = dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);

    fscanf(pFile, "%d %[^\n]\n", &itemp, buffer_line);
    ind_g_neg = arma::zeros<arma::uvec>(itemp, 1);
    coef_g_neg = arma::zeros(itemp, 1);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%d", &itemp2) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the indices of lines for negative YAW force \n";
            throw ValueError(ss.str());
        }
        ind_g_neg(ii) = itemp2 - 1;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
    for (int ii = 0; ii < itemp; ii++)
    {
        if (fscanf(pFile, "%lf", &dtemp) != 1)
        {
            std::stringstream ss;
            ss << "An error occurred when trying to read the coefficients of lines for negative YAW force \n";
            throw ValueError(ss.str());
        }
        coef_g_neg(ii) = dtemp;
    }
    fscanf(pFile, "%[^\n]\n", buffer_line);
}

void WinchieControllerHorizontal::ReadPropertiesYAML(YAML::Node node)
{
    indBody = node["body_index"].as<int>() - 1;

    Ac = node["Ac"].as<double>();
    Bc = node["Bc"].as<double>();
    Cc = node["Cc"].as<double>();
    Dc = node["Dc"].as<double>();

    YAML::Node kcNode = node["Kc"];
    for (int ii = 0; ii < 3; ii++)
        Kc(ii, 0) = kcNode[ii].as<double>();

    YAML::Node kiNode = node["Ki"];
    for (int ii = 0; ii < 3; ii++)
        Ki(ii, 0) = kiNode[ii].as<double>();

    Af = node["Af"].as<double>();
    Bf = node["Bf"].as<double>();
    Cf = node["Cf"].as<double>();
    Df = node["Df"].as<double>();

    Kw = node["Kw"].as<double>();
    time_ini = node["time_ini"].as<double>();
    reference_flag = node["reference_flag"].as<long>();

    if (reference_flag < 1)
    {
        Ar = node["Ar"].as<double>();
        Br = node["Br"].as<double>();
        Cr = node["Cr"].as<double>();
        Dr = node["Dr"].as<double>();

        YAML::Node urNode = node["ur"];
        for (int ii = 0; ii < 3; ii++)
            ur(ii, 0) = urNode[ii].as<double>();
    }
    else
    {
        t_ref = arma::zeros(reference_flag, 1);
        x_ref = arma::zeros(reference_flag, 1);
        y_ref = arma::zeros(reference_flag, 1);
        yaw_ref = arma::zeros(reference_flag, 1);

        YAML::Node tNode = node["t_ref"];
        for (int ii = 0; ii < reference_flag; ii++)
            t_ref(ii, 0) = tNode[ii].as<double>();

        YAML::Node xNode = node["x_ref"];
        for (int ii = 0; ii < reference_flag; ii++)
            x_ref(ii, 0) = xNode[ii].as<double>();

        YAML::Node yNode = node["y_ref"];
        for (int ii = 0; ii < reference_flag; ii++)
            y_ref(ii, 0) = yNode[ii].as<double>();

        YAML::Node yawNode = node["yaw_ref"];
        for (int ii = 0; ii < reference_flag; ii++)
            yaw_ref(ii, 0) = yawNode[ii].as<double>();
    }

    inversor_flag = node["inversor_flag"].as<int>();
    T_max = node["T_max"].as<double>();
    T_min = node["T_min"].as<double>();
    nIterMax = node["nIterMax"].as<int>();
    atol = node["atol"].as<double>();

    auto readInversorGroup =
        [&](const std::string& indKey, const std::string& coefKey, arma::uvec& ind, arma::mat& coef)
    {
        YAML::Node indNode = node[indKey];
        YAML::Node coefNode = node[coefKey];
        int n = (int)indNode.size();
        ind = arma::zeros<arma::uvec>(n, 1);
        coef = arma::zeros(n, 1);
        for (int ii = 0; ii < n; ii++)
            ind(ii) = indNode[ii].as<int>() - 1;
        for (int ii = 0; ii < n; ii++)
            coef(ii) = coefNode[ii].as<double>();
    };

    readInversorGroup("ind_x_pos", "coef_x_pos", ind_x_pos, coef_x_pos);
    readInversorGroup("ind_x_neg", "coef_x_neg", ind_x_neg, coef_x_neg);
    readInversorGroup("ind_y_pos", "coef_y_pos", ind_y_pos, coef_y_pos);
    readInversorGroup("ind_y_neg", "coef_y_neg", ind_y_neg, coef_y_neg);
    readInversorGroup("ind_g_pos", "coef_g_pos", ind_g_pos, coef_g_pos);
    readInversorGroup("ind_g_neg", "coef_g_neg", ind_g_neg, coef_g_neg);
}

void WinchieControllerHorizontal::SetUpWinchiesController(void)
{

    int ne = ceil(pSim->simulationTime / pSim->maxTimeStep);
    error = arma::zeros(ne, 3);

    T_min = T_min * pSim->gravity * 1000;
    T_max = T_max * pSim->gravity * 1000;
    Logger::info("Winches controller T_min = " + std::to_string(T_min) + " N");
    Logger::info("Winches controller T_max = " + std::to_string(T_max) + " N");

    T = arma::ones(nWinchies, 1) * T_min;

    coef_x_pos = arma::pow(coef_x_pos, -1);
    coef_x_pos = coef_x_pos / arma::accu(coef_x_pos);
    coef_x_neg = arma::pow(coef_x_neg, -1);
    coef_x_neg = coef_x_neg / arma::accu(coef_x_neg);
    coef_y_pos = arma::pow(coef_y_pos, -1);
    coef_y_pos = coef_y_pos / arma::accu(coef_y_pos);
    coef_y_neg = arma::pow(coef_y_neg, -1);
    coef_y_neg = coef_y_neg / arma::accu(coef_y_neg);
    coef_g_pos = arma::pow(coef_g_pos, -1);
    coef_g_pos = coef_g_pos / arma::accu(coef_g_pos);
    coef_g_neg = arma::pow(coef_g_neg, -1);
    coef_g_neg = coef_g_neg / arma::accu(coef_g_neg);

    posAnchG = arma::zeros(3, nWinchies);
    posFairL = arma::zeros(3, nWinchies);
    int indAnch, indFair;

    for (int ii = 0; ii < nWinchies; ii = ii + 1)
    {

        Winchies[ii]->tau = T_min * Winchies[ii]->radius;

        indFair = Winchies[ii]->LineBCP - 1;
        indAnch = 0;
        if (indFair == 0)
            indAnch = 1;
        posAnchG.col(ii) = Winchies[ii]->LineW->pLineBcps[indAnch]->pos;
        posFairL.col(ii) = Winchies[ii]->LineW->pLineBcps[indFair]->pos;
    }
}

void WinchieControllerHorizontal::controlWinchies(double time)
{

    arma::mat pos;

    if (time >= time_ini)
    {
        arma::mat e1;
        // Reference signal
        if (reference_flag < 1)
        {
            xr = Ar * xr + Br * ur;
            yr = Cr * xr + Dr * ur;
        }
        else
        {
            arma::mat temp_input, temp_output;
            temp_input = arma::ones(1) * time;
            arma::interp1(t_ref, x_ref, temp_input, temp_output, "linear",
                          arma::as_scalar(x_ref(reference_flag - 1, 0)));
            yr(0, 0) = arma::as_scalar(temp_output);
            arma::interp1(t_ref, y_ref, temp_input, temp_output, "linear",
                          arma::as_scalar(y_ref(reference_flag - 1, 0)));
            yr(1, 0) = arma::as_scalar(temp_output);
            arma::interp1(t_ref, yaw_ref, temp_input, temp_output, "linear",
                          arma::as_scalar(yaw_ref(reference_flag - 1, 0)));
            yr(2, 0) = arma::as_scalar(temp_output);
        }
        // Extract position from body
        // Aqui habria que meter ruido gausiano para el ruido de los sensores
        pos = pSim->pBodies[indBody]->pos;
        yb(0, 0) = pos(0, 0);
        yb(1, 0) = pos(1, 0);
        yb(2, 0) = pos(5, 0);
        // First order filter
        xf = Af * xf + Bf * yb;
        yf = Cf * xf + Df * yb;
        // Compute error
        error.row(k) = (yr - yf).t();
        // Integrate error
        e1 = Ki % arma::trapz(error.rows(0, k)).t() - yf;
        xc = Ac * xc + Bc * e1;
        yc = Kc % (Cc * xc + Dc * e1);
        k = k + 1;
        inversorBlock();
    }

    applyTensions();
}

void WinchieControllerHorizontal::inversorBlock(void)
{

    if (inversor_flag == 1)
    {

        arma::mat rotMat = pSim->pBodies[indBody]->rotMat;
        arma::mat posBody = pSim->pBodies[indBody]->pos.rows(0, 2);
        arma::mat posFairG_temp, rG;
        arma::mat Aeq = arma::zeros(3, nWinchies);
        double alpha, beta, dx, dy, dz, rx, ry, rz;

        for (int ii = 0; ii < nWinchies; ii = ii + 1)
        {
            rG = rotMat * posFairL.col(ii);
            posFairG_temp = posBody + rG;
            dx = arma::as_scalar(posAnchG(0, ii) - posFairG_temp(0, 0));
            dy = arma::as_scalar(posAnchG(1, ii) - posFairG_temp(1, 0));
            dz = arma::as_scalar(posAnchG(2, ii) - posFairG_temp(2, 0));
            rx = arma::as_scalar(rG(0, 0));
            ry = arma::as_scalar(rG(1, 0));
            rz = arma::as_scalar(rG(2, 0));
            alpha = atan2(dy, dx);
            beta = atan2(-dz, sqrt(dx * dx + dy * dy));
            Aeq(0, ii) = cos(alpha) * cos(beta);
            Aeq(1, ii) = sin(alpha) * cos(beta);
            Aeq(2, ii) = (ry * sin(beta) - rz * sin(alpha) * cos(beta)) * rotMat(0, 2) +
                         (rz * cos(alpha) * cos(beta) - rx * sin(beta)) * rotMat(1, 2) +
                         (rx * sin(alpha) * cos(beta) - ry * cos(alpha) * cos(beta)) * rotMat(2, 2);
        }

        arma::mat beq = Kw * yc;

        arma::mat A = arma::join_vert(arma::eye(nWinchies, nWinchies), arma::ones(1, nWinchies) / nWinchies);
        arma::mat b = arma::join_vert(T, arma::zeros(1, 1));
        arma::mat C = Aeq;
        arma::mat d = beq;
        arma::mat AA = A.t() * A;
        arma::mat invAA = arma::inv(AA);

        arma::mat x = invAA * (A.t() * b - C.t() * arma::inv(C * invAA * C.t()) * (C * invAA * A.t() * b - d));

        int kk = 1;
        arma::mat CC = C.t() * arma::inv(C * C.t());
        while ((arma::any(arma::any(x > T_max + atol)) || arma::any(arma::any(x < T_min - atol))) && (kk <= nIterMax))
        {
            b = arma::clamp(x, T_min, T_max);
            x = b - CC * (C * b - d);
            kk = kk + 1;
        }

        if (kk > nIterMax)
        {
            Logger::warning("Tensions clamped on winches controller!");
        }

        T = arma::clamp(x, T_min, T_max);
    }
    else if (inversor_flag == 2)
    {

        double Fx, Fy, Fg;
        Fx = arma::as_scalar(yc(0, 0));
        Fy = arma::as_scalar(yc(1, 0));
        Fg = arma::as_scalar(yc(2, 0));
        arma::mat f = arma::zeros(nWinchies, 3);
        arma::uvec ind_0 = arma::zeros<arma::uvec>(1);
        if (Fx > 0)
            f.submat(ind_x_pos, ind_0) = coef_x_pos;
        if (Fx < 0)
            f.submat(ind_x_neg, ind_0) = coef_x_neg;
        if (Fy > 0)
            f.submat(ind_y_pos, ind_0 + 1) = coef_y_pos;
        if (Fy < 0)
            f.submat(ind_y_neg, ind_0 + 1) = coef_y_neg;
        if (Fg > 0)
            f.submat(ind_g_pos, ind_0 + 2) = coef_g_pos;
        if (Fg < 0)
            f.submat(ind_g_neg, ind_0 + 2) = coef_g_neg;
        T = arma::clamp(T_min + Kw * f * arma::abs(yc), T_min, T_max);
    }
    else
    {
        std::stringstream ss;
        ss << "Inversor flag not available. \n";
        throw ValueError(ss.str());
    }
}

void WinchieControllerHorizontal::OpenOutputFilesASCII(std::string path)
{
    // Open base output files (WinchesTensions.txt, WinchedLinesLengths.txt)
    WinchieController::OpenOutputFilesASCII(path);

    // Open horizontal-controller-specific output files
    std::string file_path;

    file_path = JoinPath(path, "ControlForce.txt");
    pfile_FC = fopen(file_path.c_str(), "w");
    if (pfile_FC == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: ControlForce.txt\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }

    file_path = JoinPath(path, "ReferencePosition.txt");
    pfile_RP = fopen(file_path.c_str(), "w");
    if (pfile_RP == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: ReferencePosition.txt\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
}

void WinchieControllerHorizontal::CloseOutputFilesASCII(void)
{
    if (pfile_FC)
        fclose(pfile_FC);
    if (pfile_RP)
        fclose(pfile_RP);
    WinchieController::CloseOutputFilesASCII();
}

void WinchieControllerHorizontal::OpenOutputFilesCSV(std::string path)
{
    // Open base CSV output files (winches_tensions.csv, winches_lengths.csv)
    WinchieController::OpenOutputFilesCSV(path);

    // Open horizontal-controller-specific CSV output files
    std::string file_path;

    file_path = JoinPath(path, "control_force.csv");
    pfile_FC_csv = fopen(file_path.c_str(), "w");
    if (pfile_FC_csv == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: control_force.csv\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    fprintf(pfile_FC_csv, "time,ctrl_fx,ctrl_fy,ctrl_fz\n");

    file_path = JoinPath(path, "reference_position.csv");
    pfile_RP_csv = fopen(file_path.c_str(), "w");
    if (pfile_RP_csv == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: reference_position.csv\n    ->Dir: " << path << std::endl;
        throw IOError(ss.str());
    }
    fprintf(pfile_RP_csv, "time,ref_x,ref_y,ref_z\n");
}

void WinchieControllerHorizontal::CloseOutputFilesCSV(void)
{
    if (pfile_FC_csv)
        fclose(pfile_FC_csv);
    if (pfile_RP_csv)
        fclose(pfile_RP_csv);
    WinchieController::CloseOutputFilesCSV();
}

void WinchieControllerHorizontal::WriteOut(double t)
{
    // Write base outputs (tensions and line lengths)
    WinchieController::WriteOut(t);

    if (pSim->outputFormat == 1)
    {
        // CSV format
        fprintf(pfile_FC_csv, "%f", t);
        for (int ii = 0; ii < 3; ii = ii + 1)
            fprintf(pfile_FC_csv, ",%f", Kw * yc(ii, 0));
        fprintf(pfile_FC_csv, "\n");

        fprintf(pfile_RP_csv, "%f", t);
        for (int ii = 0; ii < 3; ii = ii + 1)
            fprintf(pfile_RP_csv, ",%f", yr(ii, 0));
        fprintf(pfile_RP_csv, "\n");
    }
    else
    {
        // ASCII format (original)
        fprintf(pfile_FC, "%f    ", t);
        for (int ii = 0; ii < 3; ii = ii + 1)
            fprintf(pfile_FC, "%f    ", Kw * yc(ii, 0));
        fprintf(pfile_FC, "\n");

        fprintf(pfile_RP, "%f    ", t);
        for (int ii = 0; ii < 3; ii = ii + 1)
            fprintf(pfile_RP, "%f    ", yr(ii, 0));
        fprintf(pfile_RP, "\n");
    }
}
