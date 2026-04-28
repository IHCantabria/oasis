// SPDX-License-Identifier: GPL-3.0-or-later

#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>

#include "Simulation.hpp"
#include "../Logger.hpp"
#include "../CommonTools.hpp"
#include "../Exceptions/Exception.hpp"
#include "../os_tools.hpp"
#include "../Bodies/Bodies.hpp"
#include "../BCPs/BCPs.hpp"
#include "../BCPs/Winchies.hpp"
#include "../BCPs/WinchiesController.hpp"
#include "../Waves/Wave.hpp"
#include "../ODE_solvers/ODE_solvers.hpp"
#include "../WindTurbine/WindTurbine.hpp"

#ifndef __has_include
static_assert(false, "__has_include not supported");
#else
#if (__cplusplus >= 201703L || (defined(_MSVC_LANG) && _MSVC_LANG >= 201703L)) && __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#elif __has_include(<boost/filesystem.hpp>)
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#endif
#endif

arma::mat Simulation::CalculateSystemDynamics(double time, arma::mat y)
{
    // TODO: Implement a logger with different levels of verbosity
    // std::cout << "Time: " << time << " s\n";
    numCallsSysFun++;
    // std::cout << "Simulation::CalculateSystemDynamics - At first" << std::endl;
    arma::mat yprime = arma::zeros(size(y));
    int i0;
    int ini;
    // Copy info from y to the objects.
    // std::cout << "Simulation::CalculateSystemDynamics - Copy info from y to objects" << std::endl;
    ini = 0;
    // Load variables from second order systems
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        for (int jj = 0; jj < this->pBodiesFree[ii]->numDofs; jj = jj + 1)
        {
            int itemp = this->pBodiesFree[ii]->pDofs[jj];
            this->pBodiesFree[ii]->pos(itemp, 0) = y(ini);
            ini = ini + 1;
        }
        // std::cout << "Simulation::CalculateSystemDynamics - Body " << this->pBodiesFree[ii]->GetId() + 1 << " pos: "
        // << this->pBodiesFree[ii]->pos.t() << std::endl;
    }
    for (int ii = 0; ii < numLines; ii++)
    {
        for (int jj = pLines[ii]->first_node; jj < pLines[ii]->last_node; jj = jj + 1)
        {
            pLines[ii]->pos.row(jj) = y.rows(ini, ini + 2).t();
            ini = ini + 3;
        }
    }
    for (int ii = 0; ii < numWinches; ii++)
    {
        pWinches[ii]->theta = arma::as_scalar(y.row(ini));
        ini = ini + 1;
    }
    for (int ii = 0; ii < numWindTurbines; ii++)
    {
        pWindTurbines[ii]->rotPos = arma::as_scalar(y.row(ini));
        ini = ini + 1;
    }
    // Load variable derivatives from second order systems
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        for (int jj = 0; jj < this->pBodiesFree[ii]->numDofs; jj = jj + 1)
        {
            int itemp = this->pBodiesFree[ii]->pDofs[jj];
            this->pBodiesFree[ii]->vel(itemp, 0) = y(ini);
            ini = ini + 1;
        }
        // std::cout << "Simulation::CalculateSystemDynamics - Body " << this->pBodiesFree[ii]->GetId() + 1 << " vel: "
        // << this->pBodiesFree[ii]->vel.t() << std::endl;
    }
    for (int ii = 0; ii < numLines; ii++)
    {
        for (int jj = pLines[ii]->first_node; jj < pLines[ii]->last_node; jj = jj + 1)
        {
            pLines[ii]->vel.row(jj) = y.rows(ini, ini + 2).t();
            ini = ini + 3;
        }
    }
    for (int ii = 0; ii < numWinches; ii++)
    {
        pWinches[ii]->omega = arma::as_scalar(y.row(ini));
        ini = ini + 1;
    }
    for (int ii = 0; ii < numWindTurbines; ii++)
    {
        pWindTurbines[ii]->rotSpeed = arma::as_scalar(y.row(ini));
        ini = ini + 1;
    }
    // Load variables from first order systems
    for (int ii = 0; ii < numOWCs; ii++)
    {
        if (this->pOWCs[ii]->turbine_type >= 0)
        {
            pOWCs[ii]->rel_pressure = arma::as_scalar(y.row(ini));
            ini = ini + 1;
        }
        if (this->pOWCs[ii]->turbine_type > 0)
        {
            int turb_valve_type = this->pOWCs[ii]->pOWCTurbine->pOWCTurbineType->valve_type;
            if (turb_valve_type == 0 || turb_valve_type == 1)
            {
                pOWCs[ii]->pOWCTurbine->angular_velocity = arma::as_scalar(y.row(ini));
                ini = ini + 1;
            }
            else if (turb_valve_type == 2 || turb_valve_type == 3)
            {
                pOWCs[ii]->pOWCTurbine->gap_rel_pressure = arma::as_scalar(y.row(ini));
                ini = ini + 1;
                pOWCs[ii]->pOWCTurbine->angular_velocity = arma::as_scalar(y.row(ini));
                ini = ini + 1;
            }
        }
    }

    for (int ii = 0; ii < numBodiesLock; ii++)
    {
        pBodiesLock[ii]->UpdateLockBody(time);
    }

    // Update BodyBCP positions and velocities
    // std::cout << "Simulation::CalculateSystemDynamics - Update BCP positions and velocities" << std::endl;
    for (int ii = 0; ii < numBodies; ii++)
    {
        pBodies[ii]->UpdateBcps();
        pBodies[ii]->ResetBcps();
    }
    // Set boundary conditions on pos and vel of Lines if the BCP is not a joint
    // std::cout << "Simulation::CalculateSystemDynamics - Set Boundary conditios" << std::endl;
    for (int ii = 0; ii < numLines; ii++)
    {
        if (pLines[ii]->pLineBcps[0]->GetType() != 3 && pLines[ii]->pLineBcps[0]->GetType() != 5)
        {
            pLines[ii]->pLineBcps[0]->GetValues(time);
            pLines[ii]->pos.row(0) = pLines[ii]->pLineBcps[0]->pos.t();
            pLines[ii]->vel.row(0) = pLines[ii]->pLineBcps[0]->vel.t();
        }
        if (pLines[ii]->pLineBcps[1]->GetType() != 3 && pLines[ii]->pLineBcps[1]->GetType() != 5)
        {
            pLines[ii]->pLineBcps[1]->GetValues(time);
            pLines[ii]->pos.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->pos.t();
            pLines[ii]->vel.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->vel.t();
        }
    }
    // Set boundary conditions on pos and vel of Lines if the BCP is a joint or elastic anchor
    // std::cout << "Simulation::CalculateSystemDynamics - Set Boundary conditios if BCP is a Joint or ElasticAnchor" <<
    // std::endl;
    for (int ii = 0; ii < numLines; ii++)
    {
        if (pLines[ii]->pLineBcps[0]->GetType() == 3 || pLines[ii]->pLineBcps[0]->GetType() == 5)
        {
            pLines[ii]->pLineBcps[0]->posLines.row(pLines[ii]->pLineBcps[0]->iLJ) = pLines[ii]->pos.row(0);
            pLines[ii]->pLineBcps[0]->velLines.row(pLines[ii]->pLineBcps[0]->iLJ) = pLines[ii]->vel.row(0);
            pLines[ii]->pLineBcps[0]->iLJ = pLines[ii]->pLineBcps[0]->iLJ + 1;
        }
        if (pLines[ii]->pLineBcps[1]->GetType() == 3 || pLines[ii]->pLineBcps[1]->GetType() == 5)
        {
            pLines[ii]->pLineBcps[1]->posLines.row(pLines[ii]->pLineBcps[1]->iLJ) =
                pLines[ii]->pos.row(pLines[ii]->N - 1);
            pLines[ii]->pLineBcps[1]->velLines.row(pLines[ii]->pLineBcps[1]->iLJ) =
                pLines[ii]->vel.row(pLines[ii]->N - 1);
            pLines[ii]->pLineBcps[1]->iLJ = pLines[ii]->pLineBcps[1]->iLJ + 1;
        }
    }
    for (int ii = 0; ii < numLines; ii++)
    {
        if (pLines[ii]->pLineBcps[0]->GetType() == 3 || pLines[ii]->pLineBcps[0]->GetType() == 5)
        {
            pLines[ii]->pLineBcps[0]->GetValues(time);
            pLines[ii]->pos.row(0) = pLines[ii]->pLineBcps[0]->pos.t();
            pLines[ii]->vel.row(0) = pLines[ii]->pLineBcps[0]->vel.t();
        }
        if (pLines[ii]->pLineBcps[1]->GetType() == 3 || pLines[ii]->pLineBcps[1]->GetType() == 5)
        {
            pLines[ii]->pLineBcps[1]->GetValues(time);
            pLines[ii]->pos.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->pos.t();
            pLines[ii]->vel.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->vel.t();
        }
    }

    // Compute forces vector for the different Lines
    // std::cout << "Simulation::CalculateSystemDynamics - Compute forces vector for different lines" << std::endl;
    for (int ii = 0; ii < numLines; ii++)
    {
        pLines[ii]->SEM_computeF(time);
    }

    // Compute forces of Springs
    // std::cout << "Simulation::CalculateSystemDynamics - Compute spring" << std::endl;
    for (int ii = 0; ii < numSprings; ii++)
    {
        pSprings[ii]->computeSpringForces();
    }

    // Update hydrostatic parameters if there is sinking
    // std::cout << "Simulation::CalculateSystemDynamics - Update Sinking Hydrostatics" << std::endl;
    for (int ii = 0; ii < numSinking; ii = ii + 1)
    {
        pSinking[ii]->UpdateSinkingHydrostatics(time);
    }

    // Initialize forces vector
    arma::mat Fb = arma::zeros(6 * numBodiesFree, 1);
    double coefHydro = ((time - lastHydroTime_old) * (time - lastHydroTime_old2)) /
                       ((lastHydroTime - lastHydroTime_old) * (lastHydroTime - lastHydroTime_old2));
    double coefHydro_old = ((time - lastHydroTime) * (time - lastHydroTime_old2)) /
                           ((lastHydroTime_old - lastHydroTime) * (lastHydroTime_old - lastHydroTime_old2));
    double coefHydro_old2 = ((time - lastHydroTime) * (time - lastHydroTime_old)) /
                            ((lastHydroTime_old2 - lastHydroTime) * (lastHydroTime_old2 - lastHydroTime_old));

    // Compute hydrostatic and hydrodynamic forces
    // std::cout << "Simulation::CalculateSystemDynamics - Compute hydrodynamic and hydrostatic forces" << std::endl;
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        // Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) = pBodies[ii]->pHydro->CalculateHydrodynamicForces(time) +
        // pBodies[ii]->pHydro->CalculateHydrostaticForces(time);
        Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) = coefHydro * pBodiesFree[ii]->Fb +
                                                      coefHydro_old * pBodiesFree[ii]->Fb_old +
                                                      coefHydro_old2 * pBodiesFree[ii]->Fb_old2;
    }
    if (Fb.has_nan() | Fb.has_inf())
    {
        Logger::error("NaN or Inf in Fb for hydrostatic/hydrodynamic forces (t=" + std::to_string(time) + ")");
        throw std::exception();
    }

    // Compute forces on BCPs
    // std::cout << "Simulation::CalculateSystemDynamics - Compute forces on BCPs" << std::endl;
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        pBodiesFree[ii]->ComputeBcpForces();
        Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) =
            Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) + pBodiesFree[ii]->bcpForces;
    }
    if (Fb.has_nan() | Fb.has_inf())
    {
        Logger::error("NaN or Inf in Fb for BCP forces (t=" + std::to_string(time) + ")");
        throw std::exception();
    }

    // Add wind turbine forces
    // std::cout << "Simulation::CalculateSystemDynamics - Add wind turbine forces" << std::endl;
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        pBodiesFree[ii]->ComputeWindTurbForces();
        Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) =
            Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) + pBodiesFree[ii]->windTurbForces;
    }
    if (Fb.has_nan() | Fb.has_inf())
    {
        Logger::error("NaN or Inf in Fb for wind turbine forces (t=" + std::to_string(time) + ")");
        throw std::exception();
    }

    // Compute OWCs dynamics
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        pBodiesFree[ii]->owcForces = arma::zeros(6, 1);
    }
    for (int ii = 0; ii < numOWCs; ii++)
    {
        pOWCs[ii]->ComputeForces(time);
    }
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) =
            Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) + pBodiesFree[ii]->owcForces;
    }

    // Compute also everything for locked bodies so it can be displayed on the output files
    for (int ii = 0; ii < numBodiesLock; ii++)
    {
        pBodiesLock[ii]->ComputeWindTurbForces();
        pBodiesLock[ii]->ComputeBcpForces();
    }

    // Compute body acceleration
    // TODO: Regarding OWCs, this could be more efficient if we remove unused DOFs from system matrix
    // std::cout << "Simulation::CalculateSystemDynamics - Compute Bodies accelerations" << std::endl;
    arma::mat accB;
    if (numBodiesFree > 0)
    {
        if (numBodiesLock > 0)
        { // If locked bodies, reduce system matrix and include effect acc from other bodies

            // Assemble a vector with all accelerations
            arma::mat accAll = arma::zeros(6 * numBodies, 1);
            for (int ii = 0; ii < numBodies; ii++)
            {
                accAll(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) = pBodies[ii]->acc;
            }

            if (rotSimpFlag)
            { // If simplification, always use same matrix

                accB =
                    (*pSystemMatrixFFInv) * (Fb.rows(sysMatIndFree) - (*pSystemMatrixFL) * accAll.rows(sysMatIndLock));
            }
            else
            { // If no simplification, update matrix for new rotation states

                arma::mat tmpMat = (*pSystemMatrix);
                arma::mat tmpVec = Fb;
                arma::mat auxM, auxV, aMat, aMat_dot, invRotMat, phi_dot;
                for (int ii = 0; ii < numBodiesFree; ii++)
                {
                    aMat = pBodiesFree[ii]->aMat;
                    aMat_dot = pBodiesFree[ii]->aMat_dot;
                    invRotMat = pBodiesFree[ii]->invRotMat;
                    phi_dot = pBodiesFree[ii]->vel.rows(arma::span(3, 5));

                    auxM = tmpMat(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), arma::span(6 * ii + 3, 6 * (ii + 1) - 1));
                    tmpMat(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), arma::span(6 * ii + 3, 6 * (ii + 1) - 1)) =
                        auxM * invRotMat * aMat;

                    auxV = tmpVec(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), 0);
                    tmpVec(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), 0) =
                        auxV - auxM * invRotMat * aMat_dot * phi_dot -
                        arma::cross(invRotMat * aMat * phi_dot, auxM * invRotMat * aMat * phi_dot);
                }

                arma::mat tmpMat_FF = tmpMat(sysMatIndFree, sysMatIndFree);
                arma::mat tmpMat_FL = tmpMat(sysMatIndFree, sysMatIndLock);
                arma::mat tmpVec_FF = tmpVec.rows(sysMatIndFree);

                accB = arma::solve(tmpMat_FF, tmpVec_FF - tmpMat_FL * accAll.rows(sysMatIndLock));
            }
        }
        else
        { // If no locked bodies, keep it simple

            if (rotSimpFlag)
            { // If simplification, always use same matrix

                accB = (*pSystemMatrixInv) * Fb;
            }
            else
            { // If no simplification, update matrix for new rotation states

                arma::mat tmpMat = (*pSystemMatrix);
                arma::mat tmpVec = Fb;
                arma::mat auxM, auxV, aMat, aMat_dot, invRotMat, phi_dot;
                for (int ii = 0; ii < numBodies; ii++)
                {
                    aMat = pBodies[ii]->aMat;
                    aMat_dot = pBodies[ii]->aMat_dot;
                    invRotMat = pBodies[ii]->invRotMat;
                    phi_dot = pBodies[ii]->vel.rows(arma::span(3, 5));

                    auxM = tmpMat(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), arma::span(6 * ii + 3, 6 * (ii + 1) - 1));
                    tmpMat(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), arma::span(6 * ii + 3, 6 * (ii + 1) - 1)) =
                        auxM * invRotMat * aMat;

                    auxV = tmpVec(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), 0);
                    tmpVec(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), 0) =
                        auxV - auxM * invRotMat * aMat_dot * phi_dot -
                        arma::cross(invRotMat * aMat * phi_dot, auxM * invRotMat * aMat * phi_dot);
                }

                accB = arma::solve(tmpMat, tmpVec);
            }
        }
    }
    if (accB.has_nan() | accB.has_inf())
    {
        Logger::error("NaN or Inf in bodies accelerations (t=" + std::to_string(time) + ")");
        throw std::exception();
    }
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        pBodiesFree[ii]->acc = accB(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) % pBodiesFree[ii]->isDofActive;
    }

    // Update BodyBCP accelerations
    // std::cout << "Simulation::CalculateSystemDynamics - Compute BCP accelerations" << std::endl;
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        pBodiesFree[ii]->UpdateBcps();
    }

    // Obtain Lines accelerations, imposing boundary conditions if the BCP is not a joint or elastic anchor
    // std::cout << "Simulation::CalculateSystemDynamics - Compute lines accelerations" << std::endl;
    for (int ii = 0; ii < numLines; ii++)
    {
        // For body BCPs, the accelertion has changed, so GetValues() routine is called again
        if (pLines[ii]->pLineBcps[0]->GetType() == 4)
        {
            pLines[ii]->pLineBcps[0]->GetValues(time);
        }
        if (pLines[ii]->pLineBcps[1]->GetType() == 4)
        {
            pLines[ii]->pLineBcps[1]->GetValues(time);
        }
        // Impose BCP accelerations on lines forces vectors (not for joints or elastic anchors)
        if (pLines[ii]->pLineBcps[0]->GetType() != 3 && pLines[ii]->pLineBcps[0]->GetType() != 5)
        {
            pLines[ii]->F.row(0) = pLines[ii]->pLineBcps[0]->acc.t();
        }
        if (pLines[ii]->pLineBcps[1]->GetType() != 3 && pLines[ii]->pLineBcps[1]->GetType() != 5)
        {
            pLines[ii]->F.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->acc.t();
        }
    }

    // Assemble lines forces vectors
    // std::cout << "Simulation::CalculateSystemDynamics - Assemble lines forces vectors" << std::endl;
    arma::mat LinesCouplingVector = arma::zeros(numAllLinesNodes, 3);
    for (int ii = 0; ii < numLines; ii++)
    {
        if (pLines[ii]->F.has_nan() | pLines[ii]->F.has_inf())
        {
            Logger::error("NaN or Inf in force vector for line " + std::to_string(ii + 1) + " (t=" + std::to_string(time) + ")");
            throw std::exception();
        }
        LinesCouplingVector.rows(pLines[ii]->ind4CouplingMat) += pLines[ii]->F;
    }
    for (int ii = 0; ii < numBcps; ii++)
    {
        if ((pBcps[ii]->GetType() == 3 || pBcps[ii]->GetType() == 5) && (pBcps[ii]->flag_assigned == 1))
        {
            LinesCouplingVector.row(pBcps[ii]->couplingMatIndex) += pBcps[ii]->JointForce;
        }
    }
    if (LinesCouplingVector.has_nan() | LinesCouplingVector.has_inf())
    {
        Logger::error("NaN or Inf in lines force coupling vector (t=" + std::to_string(time) + ")");
        throw std::exception();
    }

    // Solve lines accelerations
    // std::cout << "Simulation::CalculateSystemDynamics - Solve lines accelerations" << std::endl;
    arma::mat LinesAccelerations;
    // TODO: 100 should be a parameter
    if (numAllLinesNodes >= 100)
    {
        arma::superlu_opts opts;
        opts.allow_ugly = false;
        arma::spsolve(LinesAccelerations, *pLinesCouplingMatrix_sp, LinesCouplingVector, "superlu", opts);
    }
    else
    {
        if (useWinches)
        {
            LinesAccelerations = arma::solve(*pLinesCouplingMatrix, LinesCouplingVector);
        }
        else
        {
            LinesAccelerations = (*pLinesCouplingMatrixInv) * LinesCouplingVector;
        }
    }
    if (LinesAccelerations.has_nan() | LinesAccelerations.has_inf())
    {
        Logger::error("NaN or Inf in lines accelerations (t=" + std::to_string(time) + ")");
        throw std::exception();
    }
    for (int ii = 0; ii < numLines; ii++)
    {
        pLines[ii]->acc = LinesAccelerations.rows(pLines[ii]->ind4CouplingMat);
    }

    // Compute Winches
    // std::cout << "Simulation::CalculateSystemDynamics - Compute Winchies" << std::endl;
    for (int ii = 0; ii < numWinches; ii = ii + 1)
    {
        pWinches[ii]->computeWinchie();
    }

    // Recompute Lines coupling matrix for new dL values
    // std::cout << "Simulation::CalculateSystemDynamics - Recompute Lines coupling matrix for new dL values" <<
    // std::endl;
    if (useWinches && numJointBcps > 0)
    {
        ComputeLinesCouplingMatrix();
    }

    // Compute Wind Turbines rotor acceleration
    // std::cout << "Simulation::CalculateSystemDynamics - Compute Wind Turbines" << std::endl;
    for (int ii = 0; ii < numWindTurbines; ii = ii + 1)
    {
        pWindTurbines[ii]->ComputeRotorAcc();
    }

    // Copy info from the objects to yprime
    // std::cout << "Simulation::CalculateSystemDynamics - Copy info from objects to yprime" << std::endl;
    ini = 0;
    // Return variable derivatives from second order systems
    for (int ii = 0; ii < numBodiesFree; ii = ii + 1)
    {
        for (int jj = 0; jj < this->pBodiesFree[ii]->numDofs; jj = jj + 1)
        {
            int itemp = this->pBodiesFree[ii]->pDofs[jj];
            yprime(ini) = this->pBodiesFree[ii]->vel(itemp, 0);
            ini = ini + 1;
        }
    }
    for (int ii = 0; ii < numLines; ii = ii + 1)
    {
        for (int jj = pLines[ii]->first_node; jj < pLines[ii]->last_node; jj = jj + 1)
        {
            yprime.rows(ini, ini + 2) = pLines[ii]->vel.row(jj).t();
            ini = ini + 3;
        }
    }
    for (int ii = 0; ii < numWinches; ii = ii + 1)
    {
        yprime.row(ini) = pWinches[ii]->omega;
        ini = ini + 1;
    }
    for (int ii = 0; ii < numWindTurbines; ii = ii + 1)
    {
        yprime.row(ini) = pWindTurbines[ii]->rotSpeed;
        ini = ini + 1;
    }
    // Return variable second derivatives from second order systems
    for (int ii = 0; ii < numBodiesFree; ii = ii + 1)
    {
        for (int jj = 0; jj < this->pBodiesFree[ii]->numDofs; jj = jj + 1)
        {
            int itemp = this->pBodiesFree[ii]->pDofs[jj];
            yprime(ini) = this->pBodiesFree[ii]->acc(itemp, 0);
            ini = ini + 1;
        }
    }
    for (int ii = 0; ii < numLines; ii = ii + 1)
    {
        for (int jj = pLines[ii]->first_node; jj < pLines[ii]->last_node; jj = jj + 1)
        {
            yprime.rows(ini, ini + 2) = pLines[ii]->acc.row(jj).t();
            ini = ini + 3;
        }
    }
    for (int ii = 0; ii < numWinches; ii = ii + 1)
    {
        yprime.row(ini) = pWinches[ii]->alpha;
        ini = ini + 1;
    }
    for (int ii = 0; ii < numWindTurbines; ii = ii + 1)
    {
        yprime.row(ini) = pWindTurbines[ii]->rotAcc;
        ini = ini + 1;
    }
    // Return variable derivatives from first order systems
    for (int ii = 0; ii < numOWCs; ii++)
    {
        if (this->pOWCs[ii]->turbine_type >= 0)
        {
            yprime.row(ini) = pOWCs[ii]->rel_pressure_dot;
            ini = ini + 1;
        }
        if (this->pOWCs[ii]->turbine_type > 0)
        {
            int turb_valve_type = this->pOWCs[ii]->pOWCTurbine->pOWCTurbineType->valve_type;
            if (turb_valve_type == 0 || turb_valve_type == 1)
            {
                yprime.row(ini) = pOWCs[ii]->pOWCTurbine->angular_acceleration;
                ini = ini + 1;
            }
            else if (turb_valve_type == 2 || turb_valve_type == 3)
            {
                yprime.row(ini) = pOWCs[ii]->pOWCTurbine->gap_rel_pressure_dot;
                ini = ini + 1;
                yprime.row(ini) = pOWCs[ii]->pOWCTurbine->angular_acceleration;
                ini = ini + 1;
            }
        }
    }

    // Check if yprime has a NaN
    // std::cout << "Simulation::CalculateSystemDynamics - Check if yprime has a NaN" << std::endl;
    if (yprime.has_nan() | yprime.has_inf())
    {
        Logger::error("NaN or Inf in yprime (t=" + std::to_string(time) + ")");
        throw std::exception();
    }

    // std::cout << "Simulation::CalculateSystemDynamics - End of fcn" << std::endl;
    return yprime;
}

void Simulation::CalculateSystemDynamicsStatic(double time)
{
    // This function is called when numSystem == 0 (no degrees of freedom)
    // It computes forces without solving differential equations
    // Similar to CalculateSystemDynamics but without y vector input and without computing accelerations
    // This is useful for analyzing wave loads on fixed structures or with prescribed motions
    numCallsSysFun++;

    // Update locked bodies (those with prescribed motion)
    for (int ii = 0; ii < numBodiesLock; ii++)
    {
        pBodiesLock[ii]->UpdateLockBody(time);
    }

    // Update BodyBCP positions and velocities
    for (int ii = 0; ii < numBodies; ii++)
    {
        pBodies[ii]->UpdateBcps();
        pBodies[ii]->ResetBcps();
    }

    // Set boundary conditions on pos and vel of Lines if the BCP is not a joint
    for (int ii = 0; ii < numLines; ii++)
    {
        if (pLines[ii]->pLineBcps[0]->GetType() != 3 && pLines[ii]->pLineBcps[0]->GetType() != 5)
        {
            pLines[ii]->pLineBcps[0]->GetValues(time);
            pLines[ii]->pos.row(0) = pLines[ii]->pLineBcps[0]->pos.t();
            pLines[ii]->vel.row(0) = pLines[ii]->pLineBcps[0]->vel.t();
        }
        if (pLines[ii]->pLineBcps[1]->GetType() != 3 && pLines[ii]->pLineBcps[1]->GetType() != 5)
        {
            pLines[ii]->pLineBcps[1]->GetValues(time);
            pLines[ii]->pos.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->pos.t();
            pLines[ii]->vel.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->vel.t();
        }
    }

    // Set boundary conditions on pos and vel of Lines if the BCP is a joint or elastic anchor
    for (int ii = 0; ii < numLines; ii++)
    {
        if (pLines[ii]->pLineBcps[0]->GetType() == 3 || pLines[ii]->pLineBcps[0]->GetType() == 5)
        {
            pLines[ii]->pLineBcps[0]->posLines.row(pLines[ii]->pLineBcps[0]->iLJ) = pLines[ii]->pos.row(0);
            pLines[ii]->pLineBcps[0]->velLines.row(pLines[ii]->pLineBcps[0]->iLJ) = pLines[ii]->vel.row(0);
            pLines[ii]->pLineBcps[0]->iLJ = pLines[ii]->pLineBcps[0]->iLJ + 1;
        }
        if (pLines[ii]->pLineBcps[1]->GetType() == 3 || pLines[ii]->pLineBcps[1]->GetType() == 5)
        {
            pLines[ii]->pLineBcps[1]->posLines.row(pLines[ii]->pLineBcps[1]->iLJ) =
                pLines[ii]->pos.row(pLines[ii]->N - 1);
            pLines[ii]->pLineBcps[1]->velLines.row(pLines[ii]->pLineBcps[1]->iLJ) =
                pLines[ii]->vel.row(pLines[ii]->N - 1);
            pLines[ii]->pLineBcps[1]->iLJ = pLines[ii]->pLineBcps[1]->iLJ + 1;
        }
    }

    for (int ii = 0; ii < numLines; ii++)
    {
        if (pLines[ii]->pLineBcps[0]->GetType() == 3 || pLines[ii]->pLineBcps[0]->GetType() == 5)
        {
            pLines[ii]->pLineBcps[0]->GetValues(time);
            pLines[ii]->pos.row(0) = pLines[ii]->pLineBcps[0]->pos.t();
            pLines[ii]->vel.row(0) = pLines[ii]->pLineBcps[0]->vel.t();
        }
        if (pLines[ii]->pLineBcps[1]->GetType() == 3 || pLines[ii]->pLineBcps[1]->GetType() == 5)
        {
            pLines[ii]->pLineBcps[1]->GetValues(time);
            pLines[ii]->pos.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->pos.t();
            pLines[ii]->vel.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->vel.t();
        }
    }

    // Compute forces vector for the different Lines
    for (int ii = 0; ii < numLines; ii++)
    {
        pLines[ii]->SEM_computeF(time);
    }

    // Compute forces of Springs
    for (int ii = 0; ii < numSprings; ii++)
    {
        pSprings[ii]->computeSpringForces();
    }

    // Update hydrostatic parameters if there is sinking
    for (int ii = 0; ii < numSinking; ii = ii + 1)
    {
        pSinking[ii]->UpdateSinkingHydrostatics(time);
    }

    // Compute hydrostatic and hydrodynamic forces for bodies
    for (int ii = 0; ii < numBodies; ii++)
    {
        arma::mat tmp_hs = pBodies[ii]->pHydro->CalculateHydrostaticForces(time);
        arma::mat tmp_hd = pBodies[ii]->Fb;
        // Store combined forces in body for output purposes
        pBodies[ii]->Fb = tmp_hs + tmp_hd;
    }

    // Compute forces on BCPs for all bodies
    for (int ii = 0; ii < numBodies; ii++)
    {
        pBodies[ii]->ComputeBcpForces();
    }

    // Add wind turbine forces
    for (int ii = 0; ii < numBodies; ii++)
    {
        pBodies[ii]->ComputeWindTurbForces();
    }

    // Compute OWCs dynamics
    for (int ii = 0; ii < numBodies; ii++)
    {
        pBodies[ii]->owcForces = arma::zeros(6, 1);
    }
    for (int ii = 0; ii < numOWCs; ii++)
    {
        pOWCs[ii]->ComputeForces(time);
    }
}

void Simulation::CloseCase()
{
    for (int ii = 0; ii < numLines; ii++)
    {
        pLines[ii]->CloseOutputFiles();
    }

    for (int ii = 0; ii < numBodies; ii++)
    {
        pBodies[ii]->CloseOutputFiles();
    }

    for (int ii = 0; ii < numSprings; ii++)
    {
        pSprings[ii]->CloseOutputFiles();
    }

    for (int ii = 0; ii < numSinking; ii++)
    {
        pSinking[ii]->CloseOutputFiles();
    }

    if (useWinches && WinchesController != nullptr)
    {
        WinchesController->CloseOutputFiles();
        delete WinchesController;
        WinchesController = nullptr;
    }

    for (int ii = 0; ii < numWindTurbines; ii++)
    {
        pWindTurbines[ii]->Finalize();
    }
}

void Simulation::Initialize()
{
    double start_time = 0.0;

    // Initialize system vector
    Logger::debug("Num. Bodies Free    : " + std::to_string(this->numBodiesFree));
    Logger::debug("Num. Mooring DOFs   : " + std::to_string(this->numDofLinesTotal));
    Logger::debug("Num. Winches        : " + std::to_string(this->numWinches));
    Logger::debug("Num. Wind Turbines  : " + std::to_string(this->numWindTurbines));
    Logger::debug("Num. OWCs           : " + std::to_string(this->numOWCs));

    // Get the number of DOFs of second order systems
    numSystem2 = 0;
    // Add DOFs of free bodies
    numDofBodiesFree = 0;
    for (int ii = 0; ii < this->numBodiesFree; ii++)
    {
        numDofBodiesFree = numDofBodiesFree + this->pBodiesFree[ii]->numDofs;
    }
    numSystem2 = numSystem2 + numDofBodiesFree;
    // Add DOFs of lines
    numSystem2 = numSystem2 + 3 * this->numDofLinesTotal;
    // Add DOFs of winches
    numSystem2 = numSystem2 + this->numWinches;
    // Add DOFs of wind turbines
    numSystem2 = numSystem2 + this->numWindTurbines;

    // Compute the total number of DOFs for the second order system multiplying by 2
    numSystem = 2 * numSystem2;

    // Add DOFs of first order systems at the end of the system vector
    // Add DOFs of OWCs
    for (int ii = 0; ii < this->numOWCs; ii++)
    {
        if (this->pOWCs[ii]->turbine_type >= 0)
        {
            numSystem = numSystem + 1; // Chamber relative pressure
        }
        if (this->pOWCs[ii]->turbine_type > 0)
        {
            int valve_type = this->pOWCs[ii]->pOWCTurbine->pOWCTurbineType->valve_type;
            if (valve_type == 0 || valve_type == 1) // Cases with no throttle valve
            {
                numSystem = numSystem + 1; // Turbine angular velocity
            }
            else if (valve_type == 2 || valve_type == 3) // Cases with throttle valve
            {
                numSystem = numSystem + 2; // Turbine angular velocity and gap relative pressure
            }
            else
            {
                Logger::error("OWC Turbine valve type not implemented!");
                throw std::exception();
            }
        }
    }

    Logger::debug("System size: " + std::to_string(numSystem));

    arma::mat y = arma::zeros(numSystem, 1);
    arma::mat yprime = arma::zeros(numSystem, 1);

    int ini = 0;
    std::string file_path;
    Logger::debug("Initialising system vector...");
    if (this->readEquilibrium == 0)
    {
        for (int ii = 0; ii < this->numBodiesFree; ii = ii + 1)
        {
            Logger::debug("  including body " + std::to_string(ii + 1));
            for (int jj = 0; jj < this->pBodiesFree[ii]->numDofs; jj = jj + 1)
            {
                int itemp = this->pBodiesFree[ii]->pDofs[jj];
                y(ini) = this->pBodiesFree[ii]->pos(itemp, 0);
                ini = ini + 1;
            }
        }
        for (int ii = 0; ii < this->numLines; ii = ii + 1)
        {
            Logger::debug("  including line " + std::to_string(ii + 1));
            for (int jj = pLines[ii]->first_node; jj < pLines[ii]->last_node; jj = jj + 1)
            {
                y.rows(ini, ini + 2) = this->pLines[ii]->pos.row(jj).t();
                ini = ini + 3;
            }
        }
        for (int ii = 0; ii < this->numWinches; ii = ii + 1)
        {
            Logger::debug("  including winch " + std::to_string(ii + 1));
            ini = ini + 1;
        }
        for (int ii = 0; ii < this->numWindTurbines; ii = ii + 1)
        {
            Logger::debug("  including wind turbine " + std::to_string(ii + 1));
            y(ini) = this->pWindTurbines[ii]->rotPos;
            y(numSystem2 + ini) = this->pWindTurbines[ii]->rotSpeed;
            ini = ini + 1;
        }
        ini = 0;
        for (int ii = 0; ii < this->numOWCs; ii = ii + 1)
        {
            Logger::debug("  including OWC " + std::to_string(ii + 1));
            if (this->pOWCs[ii]->turbine_type >= 0)
            {
                y(2 * numSystem2 + ini) = this->pOWCs[ii]->rel_pressure;
                ini = ini + 1;
            }
            if (this->pOWCs[ii]->turbine_type > 0)
            {
                int valve_type = this->pOWCs[ii]->pOWCTurbine->pOWCTurbineType->valve_type;
                if (valve_type == 0 || valve_type == 1) // Cases with no throttle valve
                {
                    y(2 * numSystem2 + ini) = this->pOWCs[ii]->pOWCTurbine->angular_velocity;
                    ini = ini + 1;
                }
                else if (valve_type == 2 || valve_type == 3) // Cases with throttle valve
                {
                    y(2 * numSystem2 + ini) = this->pOWCs[ii]->pOWCTurbine->gap_rel_pressure;
                    ini = ini + 1;
                    y(2 * numSystem2 + ini) = this->pOWCs[ii]->pOWCTurbine->angular_velocity;
                    ini = ini + 1;
                }
            }
        }
    }
    else
    {
        Logger::debug("Reading initial condition from dataStaticIC.dat...");
        file_path = JoinPath(inputFolderPath, "dataStaticIC.dat");
        y.load(file_path, arma::arma_ascii);

        ini = 6 * this->numBodiesFree;
        for (int ii = 0; ii < this->numLines; ii = ii + 1)
        {
            for (int jj = this->pLines[ii]->first_node; jj < this->pLines[ii]->last_node; jj = jj + 1)
            {
                this->pLines[ii]->pos.row(jj) = y.rows(ini, ini + 2).t();
                ini = ini + 3;
            }
            if (this->pLines[ii]->first_node == 0)
            {
                this->pLines[ii]->pLineBcps[0]->pos = this->pLines[ii]->pos.row(0).t();
            }
            else
            {
                this->pLines[ii]->pos.row(0) = this->pLines[ii]->pLineBcps[0]->pos.t();
            }
            if (this->pLines[ii]->last_node == this->pLines[ii]->N)
            {
                this->pLines[ii]->pLineBcps[1]->pos = this->pLines[ii]->pos.row(this->pLines[ii]->N - 1).t();
            }
            else
            {
                this->pLines[ii]->pos.row(this->pLines[ii]->N - 1) = this->pLines[ii]->pLineBcps[1]->pos.t();
            }
        }
    }

    // Initialize Temporal Solver
    if (numSystem > 0)
    {
        if (this->timeIntMethod == 1)
        {
            Logger::info("Initializing BDF2 temporal solver...");
            pTimeSolver = new BDF2(start_time, this->simulationTime, this->maxTimeStep, this->writeTimeStep, y, this);
            Logger::debug("  BDF2 constructed");
            pTimeSolver->init();
            pTimeSolver->atol = this->timeIntAbsTol;
            pTimeSolver->rtol = this->timeIntRelTol;
            pTimeSolver->nIterMax = this->maxIterStep;
            Logger::debug("  BDF2 initialised");
        }
        else if (this->timeIntMethod == 2)
        {
            Logger::info("Initializing BDF" + std::to_string(timeIntOrder) + " temporal solver...");
            pTimeSolver = new BDFN(this->timeIntOrder, this->timeIntAdaptivity, start_time, this->simulationTime,
                                   this->maxTimeStep, this->writeTimeStep, y, this);
            Logger::debug("  BDF" + std::to_string(timeIntOrder) + " constructed");
            pTimeSolver->init();
            pTimeSolver->atol = this->timeIntAbsTol;
            pTimeSolver->rtol = this->timeIntRelTol;
            pTimeSolver->nIterMax = this->maxIterStep;
            Logger::debug("  BDF" + std::to_string(timeIntOrder) + " initialised");
        }
        else if (this->timeIntMethod == 3)
        {
            Logger::info("Initializing ESDIRK temporal solver...");
            pTimeSolver = new ESDIRK(this->timeIntAdaptivity, start_time, this->simulationTime, this->maxTimeStep,
                                     this->writeTimeStep, this->timeIntJacNumStepsMax, y, this);
            Logger::debug("  ESDIRK constructed");
            pTimeSolver->atol = this->timeIntAbsTol;
            pTimeSolver->rtol = this->timeIntRelTol;
            pTimeSolver->nIterMax = this->maxIterStep;
            Logger::debug("  ESDIRK initialised");
        }
        else
        {
            Logger::error("Time integration method not recognised (timeIntMethod=" + std::to_string(timeIntMethod) + ")");
            throw std::exception();
        }
    }
    else
    {
        Logger::warning("numSystem == 0: no ODE solver needed. Using simple time-stepping.");
        pTimeSolver = nullptr;
    }

    // Write initial condition to files
    for (int ii = 0; ii < this->numLines; ii = ii + 1)
        this->pLines[ii]->WriteOut(start_time);
    for (int ii = 0; ii < this->numBodies; ii = ii + 1)
        this->pBodies[ii]->WriteOut(start_time);

    Logger::debug("Updating system...");
    this->UpdateSystem(start_time);
    Logger::debug("System updated.");
}

void Simulation::LoadCase()
{
    fs::remove_all(outputFolderPath);
    fs::create_directories(outputFolderPath);

    // Read Simulation Properties
    this->ReadProperties();

    // Read Components Data
    // TODO: For unnecessary components, if the file is not found, it should not throw an error
    this->ReadWaves();
    this->ReadBodies();
    this->ReadLines();
    this->ReadBcps();
    if (useWinches)
    {
        this->ReadWinches();
    }
    this->ReadSprings();
    this->ReadSinking();
    this->ReadSeaFloor();
    this->ReadWindTurbines();
    this->ReadOWCs();

    // Setup case
    this->SetupCase();

    // TODO: Implement a logger with different levels of verbosity
    // this->PrintSetup();
}

void Simulation::PrintSetup(void)
{
    // Print BCP properties
    for (int i = 0; i < this->numBcps; i++)
    {
        pBcps[i]->Print();
    }
}

void Simulation::ReadBcps()
{
    (this->*pReadBcps)();
}

void Simulation::ReadBcpsASCII()
{
    Logger::info("--> Reading BCPs Properties");

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataBCPs.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        Logger::warning("dataBCPs.dat not found. Setting numBCPs = 0.");
        numFairBcps = 0;
        numAnchorBcps = 0;
        numJointBcps = 0;
        numBodyBcps = 0;
        numBcps = 0;
        pBcps = new BCP*[0];
        pFairleadBcps = new FairleadBCP*[0];
        pAnchorBcps = new AnchorBCP*[0];
        pJointBcps = new JointBCP*[0];
        pBodyBcps = new BodyBCP*[0];
        return;
    }

    // Read data
    fscanf(file_pointer, "%d %[^\n]\n", &numFairBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numAnchorBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numJointBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numBodyBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numElasticAnchorBcps, bufferLine);
    numBcps = numFairBcps + numAnchorBcps + numJointBcps + numBodyBcps + numElasticAnchorBcps;

    // printf("Number of fairleads: %d\n", numFairBcps);
    // printf("Number of anchor: %d\n", numAnchorBcps);
    // printf("Number of joint: %d\n", numJointBcps);
    // printf("Number of body: %d\n", numBodyBcps);

    // Allocate a vector of BCP pointers
    pBcps = new BCP*[numBcps];
    int bcp_count = 0;

    // Read each BCP
    pFairleadBcps = new FairleadBCP*[numFairBcps];
    for (int ii = 0; ii < numFairBcps; ii++)
    {
        pFairleadBcps[ii] = new FairleadBCP(bcp_count);
        pFairleadBcps[ii]->ReadPropertiesASCII(file_pointer, inputFolderPath);
        dynamic_cast<FairleadBCP*>(pFairleadBcps[ii])->Initialize(inputFolderPath);
        pBcps[bcp_count] = pFairleadBcps[ii];
        bcp_count++;
    }
    pAnchorBcps = new AnchorBCP*[numAnchorBcps];
    for (int ii = 0; ii < numAnchorBcps; ii++)
    {
        pAnchorBcps[ii] = new AnchorBCP(bcp_count);
        pAnchorBcps[ii]->ReadPropertiesASCII(file_pointer);
        pBcps[bcp_count] = pAnchorBcps[ii];
        bcp_count++;
    }
    pJointBcps = new JointBCP*[numJointBcps];
    for (int ii = 0; ii < numJointBcps; ii++)
    {
        pJointBcps[ii] = new JointBCP(bcp_count);
        pJointBcps[ii]->ReadPropertiesASCII(file_pointer);
        dynamic_cast<JointBCP*>(pJointBcps[ii])->Initialize(this->gravity, this->waterDensity, this->waterDepth);
        pBcps[bcp_count] = pJointBcps[ii];
        bcp_count++;
    }
    pBodyBcps = new BodyBCP*[numBodyBcps];
    for (int ii = 0; ii < numBodyBcps; ii++)
    {
        pBodyBcps[ii] = new BodyBCP(bcp_count);
        pBodyBcps[ii]->ReadPropertiesASCII(file_pointer);
        pBcps[bcp_count] = pBodyBcps[ii];
        bcp_count++;
    }
    pElasticAnchorBcps = new ElasticAnchorBCP*[numElasticAnchorBcps];
    for (int ii = 0; ii < numElasticAnchorBcps; ii++)
    {
        pElasticAnchorBcps[ii] = new ElasticAnchorBCP(bcp_count);
        pElasticAnchorBcps[ii]->ReadPropertiesASCII(file_pointer);
        dynamic_cast<ElasticAnchorBCP*>(pElasticAnchorBcps[ii])->Initialize(this->gravity, this->waterDensity);
        pBcps[bcp_count] = pElasticAnchorBcps[ii];
        bcp_count++;
    }

    // Loop over BCPs to check if it is necessary to read the winches file
    for (int ii = 0; ii < numBcps; ii++)
    {
        if (pBcps[ii]->winchId != 0)
        {
            useWinches = true;
            break;
        }
    }

    // Close file
    fclose(file_pointer);
    Logger::info("--> BCPs Properties Read");
}

void Simulation::ReadBcpsYAML()
{
    Logger::info("--> Reading BCPs Properties");

    if (!yamlRoot["bcps"])
    {
        Logger::warning("'bcps' section not found in YAML. Setting numBCPs = 0.");
        numFairBcps = 0;
        numAnchorBcps = 0;
        numJointBcps = 0;
        numBodyBcps = 0;
        numElasticAnchorBcps = 0;
        numBcps = 0;
        pBcps = new BCP*[0];
        pFairleadBcps = new FairleadBCP*[0];
        pAnchorBcps = new AnchorBCP*[0];
        pJointBcps = new JointBCP*[0];
        pBodyBcps = new BodyBCP*[0];
        pElasticAnchorBcps = new ElasticAnchorBCP*[0];
        return;
    }

    YAML::Node bcpsNode = yamlRoot["bcps"];

    // Count BCPs by type from YAML arrays
    numFairBcps = bcpsNode["fairleads"] ? (int)bcpsNode["fairleads"].size() : 0;
    numAnchorBcps = bcpsNode["anchors"] ? (int)bcpsNode["anchors"].size() : 0;
    numJointBcps = bcpsNode["joints"] ? (int)bcpsNode["joints"].size() : 0;
    numBodyBcps = bcpsNode["body_bcps"] ? (int)bcpsNode["body_bcps"].size() : 0;
    numElasticAnchorBcps = bcpsNode["elastic_anchors"] ? (int)bcpsNode["elastic_anchors"].size() : 0;
    numBcps = numFairBcps + numAnchorBcps + numJointBcps + numBodyBcps + numElasticAnchorBcps;

    pBcps = new BCP*[numBcps];
    int bcp_count = 0;

    pFairleadBcps = new FairleadBCP*[numFairBcps];
    for (int ii = 0; ii < numFairBcps; ii++)
    {
        pFairleadBcps[ii] = new FairleadBCP(bcp_count);
        pFairleadBcps[ii]->ReadPropertiesYAML(bcpsNode["fairleads"][ii], inputFolderPath);
        dynamic_cast<FairleadBCP*>(pFairleadBcps[ii])->Initialize(inputFolderPath);
        pBcps[bcp_count] = pFairleadBcps[ii];
        bcp_count++;
    }
    pAnchorBcps = new AnchorBCP*[numAnchorBcps];
    for (int ii = 0; ii < numAnchorBcps; ii++)
    {
        pAnchorBcps[ii] = new AnchorBCP(bcp_count);
        pAnchorBcps[ii]->ReadPropertiesYAML(bcpsNode["anchors"][ii]);
        pBcps[bcp_count] = pAnchorBcps[ii];
        bcp_count++;
    }
    pJointBcps = new JointBCP*[numJointBcps];
    for (int ii = 0; ii < numJointBcps; ii++)
    {
        pJointBcps[ii] = new JointBCP(bcp_count);
        pJointBcps[ii]->ReadPropertiesYAML(bcpsNode["joints"][ii]);
        dynamic_cast<JointBCP*>(pJointBcps[ii])->Initialize(this->gravity, this->waterDensity, this->waterDepth);
        pBcps[bcp_count] = pJointBcps[ii];
        bcp_count++;
    }
    pBodyBcps = new BodyBCP*[numBodyBcps];
    for (int ii = 0; ii < numBodyBcps; ii++)
    {
        pBodyBcps[ii] = new BodyBCP(bcp_count);
        pBodyBcps[ii]->ReadPropertiesYAML(bcpsNode["body_bcps"][ii]);
        pBcps[bcp_count] = pBodyBcps[ii];
        bcp_count++;
    }
    pElasticAnchorBcps = new ElasticAnchorBCP*[numElasticAnchorBcps];
    for (int ii = 0; ii < numElasticAnchorBcps; ii++)
    {
        pElasticAnchorBcps[ii] = new ElasticAnchorBCP(bcp_count);
        pElasticAnchorBcps[ii]->ReadPropertiesYAML(bcpsNode["elastic_anchors"][ii]);
        dynamic_cast<ElasticAnchorBCP*>(pElasticAnchorBcps[ii])->Initialize(this->gravity, this->waterDensity);
        pBcps[bcp_count] = pElasticAnchorBcps[ii];
        bcp_count++;
    }

    // Check if any BCP needs winches
    for (int ii = 0; ii < numBcps; ii++)
    {
        if (pBcps[ii]->winchId != 0)
        {
            useWinches = true;
            break;
        }
    }

    Logger::info("--> BCPs Properties Read");
}

void Simulation::ReadBodies()
{
    (this->*pReadBodies)();
}

void Simulation::ReadBodiesASCII()
{
    Logger::info("--> Reading Bodies Properties");

    // Declare local variables
    int body_count = 0;
    char buffer_line[1000];
    int diff_count = 0;
    int hydro_database_count = 0;
    std::string hydro_databases_name[300];
    int pos_database = 0;

    // Parse file in order to guess the number of bodies
    Logger::debug("Parsing dataBodies.dat...");
    std::string file_path = JoinPath(inputFolderPath, "dataBodies.dat");
    // Check if file exists first
    FILE* test_file = fopen(file_path.c_str(), "r");
    if (test_file == NULL)
    {
        Logger::warning("dataBodies.dat not found. Setting numBodies = 0.");
        this->numBodies = 0;
    }
    else
    {
        fclose(test_file);
        this->numBodies = parse_file(file_path);
    }
    Logger::debug("Number of bodies: " + std::to_string(this->numBodies));

    if (numBodies > 0)
    {
        // Open file
        Logger::debug("Opening dataBodies.dat...");
        FILE* pFile = fopen(file_path.c_str(), "r");
        if (pFile == NULL)
        {
            std::stringstream ss;
            ss << "Not possible to open the file: dataBodies.dat\n    ->Dir: " << inputFolderPath << std::endl;
            throw IOError(ss.str());
        }

        // Read all bodies
        pBodies = new Body*[numBodies];
        for (int ii = 0; ii < numBodies; ii++)
        {
            Logger::debug("  body " + std::to_string(ii + 1));
            for (int ii = 0; ii < 3; ii++)
            {
                fgets(buffer_line, sizeof(buffer_line), pFile);
            }

            // Get body type line
            fgets(buffer_line, sizeof(buffer_line), pFile);

            // Read body
            if (strncmp(buffer_line, "RAD_DIFF", 8) == 0)
            {
                pBodies[ii] = new Body(ii, this);
                pBodies[ii]->ReadPropertiesASCII(pFile);
                pBodies[ii]->OpenOutputFiles(outputFolderPath);
            }
            else
            {
                std::stringstream ss;
                ss << "Error while parsing file: dataBodies.dat\n --> Expected body: " << ii << " type definition\n";
                throw ValueError(ss.str());
            }
        }

        // Close file
        fclose(pFile);
        Logger::debug("--> Bodies loaded");

        // Loop over bodies in order to get the number of hydrodynamic databases
        hydro_databases_name[hydro_database_count] = pBodies[0]->hydroDatabaseName;
        hydro_database_count++;
        for (int ii = 1; ii < this->numBodies; ii++)
        {
            diff_count = 0;
            for (int jj = 0; jj < hydro_database_count; jj++)
            {
                if (pBodies[ii]->hydroDatabaseName.compare(hydro_databases_name[jj]) != 0)
                {
                    diff_count++;
                }
            }

            if (diff_count == hydro_database_count)
            {
                hydro_databases_name[hydro_database_count] = pBodies[ii]->hydroDatabaseName;
                hydro_database_count++;
            }
        }

        Logger::debug("--> HDB files:");
        for (int ii = 0; ii < hydro_database_count; ii++)
        {
            Logger::debug("    " + hydro_databases_name[ii]);
        }

        // Pre-determine the number of bodies in each HDB file
        int* numBodiesPerHDB = new int[hydro_database_count];
        for (int ii = 0; ii < hydro_database_count; ii++)
        {
            std::string hdb_path = JoinPath(inputFolderPath, hydro_databases_name[ii]);
            numBodiesPerHDB[ii] = HydroDatabase::GetNumBodiesFromFile(hdb_path);
            Logger::debug("    HDB: " + hydro_databases_name[ii] + " — " + std::to_string(numBodiesPerHDB[ii]) + " bodies");
        }

        // Arrange all the bodies by database (different logic for single vs multi-body HDBs)
        Body** pBodiesSort = new Body*[numBodies];
        int* pBody_found = new int[numBodies];
        for (int ii = 0; ii < numBodies; ii++)
        {
            pBody_found[ii] = 0;
        }
        for (int ii = 0; ii < hydro_database_count; ii++)
        {
            if (numBodiesPerHDB[ii] == 1)
            {
                // Single-body HDB: multiple OASIS bodies can share it, add sequentially
                for (int jj = 0; jj < numBodies; jj++)
                {
                    if ((hydro_databases_name[ii].compare(pBodies[jj]->hydroDatabaseName) == 0) &&
                        (pBody_found[jj] == 0))
                    {
                        pBody_found[jj] = 1;
                        pBodiesSort[body_count] = pBodies[jj];
                        body_count++;
                    }
                }
            }
            else
            {
                // Multi-body HDB: bodies must be placed by their hydroDatabaseIndex
                int temp_nB_hdb = numBodiesPerHDB[ii];
                int temp_nB_found = 0;
                arma::uvec temp_hdb_ind(temp_nB_hdb);
                for (int jj = 0; jj < numBodies; jj++)
                {
                    if ((hydro_databases_name[ii].compare(pBodies[jj]->hydroDatabaseName) == 0) &&
                        (pBody_found[jj] == 0))
                    {
                        pBody_found[jj] = 1;
                        pBodiesSort[body_count + pBodies[jj]->hydroDatabaseIndex] = pBodies[jj];
                        if (temp_nB_found < temp_nB_hdb)
                        {
                            temp_hdb_ind(temp_nB_found) = pBodies[jj]->hydroDatabaseIndex;
                        }
                        temp_nB_found++;
                    }
                }
                // Validate: all body indices in the multi-body HDB must be used
                if (temp_nB_found != temp_nB_hdb)
                {
                    std::stringstream ss;
                    ss << "ERROR: Not all bodies in multi-body HDB " << hydro_databases_name[ii]
                       << " are used. Expected " << temp_nB_hdb << " but found " << temp_nB_found << ".\n";
                    throw ValueError(ss.str());
                }
                // Validate: no repeated body indices in multi-body HDB
                arma::uvec temp_hdb_ind_unique = arma::unique(temp_hdb_ind);
                if ((int)temp_hdb_ind_unique.n_elem != temp_nB_hdb)
                {
                    std::stringstream ss;
                    ss << "ERROR: Repeated body indices in multi-body HDB: " << hydro_databases_name[ii] << "\n";
                    throw ValueError(ss.str());
                }
                body_count += temp_nB_hdb;
            }
        }
        delete[] pBody_found;
        for (int ii = 0; ii < numBodies; ii++)
        {
            pBodies[ii] = pBodiesSort[ii];
        }
        delete[] pBodiesSort;

        // Checking free bodies
        Logger::debug("Checking for free bodies (ASCII)...");
        for (int ii = 0; ii < numBodies; ii++)
        {
            if (pBodies[ii]->flag_blocked > 0)
            {
                numBodiesLock++;
            }
            else
            {
                numBodiesFree++;
            }
        }
        pBodiesLock = new Body*[numBodiesLock];
        pBodiesFree = new Body*[numBodiesFree];
        int indL = 0;
        int indF = 0;
        for (int ii = 0; ii < numBodies; ii++)
        {
            if (pBodies[ii]->flag_blocked > 0)
            {
                pBodiesLock[indL] = pBodies[ii];
                indL++;
            }
            else
            {
                pBodiesFree[indF] = pBodies[ii];
                indF++;
            }
        }

        // Set hydrodynamic database to each body
        std::string hydro_file_path;
        for (int ii = 0; ii < this->numBodies; ii++)
        {
            // Look for position of the database
            pos_database = 0;
            while (true)
            {
                if (this->pBodies[ii]->hydroDatabaseName.compare(hydro_databases_name[pos_database]) == 0)
                {
                    break;
                }
                pos_database++;
            }

            int nBodiesInHDB = numBodiesPerHDB[pos_database];
            if (nBodiesInHDB == 1)
            {
                // Single-body HDB: create a private 1-element body array for this body
                Body** pSingleBodyArray = new Body*[1];
                pSingleBodyArray[0] = pBodies[ii];
                this->pBodies[ii]->LoadHydrodynamicDatabase(pSingleBodyArray, 0);
            }
            else
            {
                // Multi-body HDB: build shared body array indexed by hydroDatabaseIndex
                Body** pMultiBodyArray = new Body*[nBodiesInHDB];
                for (int jj = 0; jj < nBodiesInHDB; jj++)
                {
                    pMultiBodyArray[jj] = pBodies[ii - pBodies[ii]->hydroDatabaseIndex + jj];
                }
                this->pBodies[ii]->LoadHydrodynamicDatabase(pMultiBodyArray, pBodies[ii]->hydroDatabaseIndex);
            }
        }

        // Fill System Matrix
        Logger::debug("Filling system matrix...");
        arma::span a1;
        arma::span a2;
        this->pSystemMatrix = new arma::mat(6 * this->numBodies, 6 * this->numBodies, arma::fill::zeros);
        this->pSystemMatrixInv = new arma::mat(6 * this->numBodies, 6 * this->numBodies, arma::fill::zeros);
        for (int ii = 0; ii < this->numBodies; ii++)
        {
            // Body rows
            a1 = arma::span(6 * ii, 6 * ii + 5);

            // Get number of bodies in body's HDB
            int temp_nB_hdb = pBodies[ii]->pHydro->GetNumBodies();
            if (temp_nB_hdb == 1)
            {
                // Single-body HDB: diagonal block only, no cross-coupling
                a2 = a1;
            }
            else
            {
                // Multi-body HDB: full cross-coupling span
                int temp_indHDB = pBodies[ii]->hydroDatabaseIndex;
                a2 = arma::span(6 * (ii - temp_indHDB), 6 * (ii - temp_indHDB + temp_nB_hdb) - 1);
            }

            // Fill system matrix
            Logger::debug("  body " + std::to_string(ii + 1) + " added");
            (*pSystemMatrix)(a1, a2) += pBodies[ii]->pHydro->GetTotalMass();
            pBodies[ii]->sysMatSpan1 = a1;
            pBodies[ii]->sysMatSpan2 = a2;
            pBodies[ii]->sysMatInd1 = arma::regspace<arma::uvec>(6 * ii, 6 * ii + 5);
        }

        // Take free dofs index vectors from the system matrix
        sysMatIndFree = arma::zeros<arma::uvec>(6 * numBodiesFree);
        for (int ii = 0; ii < this->numBodiesFree; ii++)
        {
            sysMatIndFree(arma::span(6 * ii, 6 * (ii + 1) - 1)) = pBodiesFree[ii]->sysMatInd1;
        }
        sysMatIndLock = arma::zeros<arma::uvec>(6 * numBodiesLock);
        for (int ii = 0; ii < this->numBodiesLock; ii++)
        {
            sysMatIndLock(arma::span(6 * ii, 6 * (ii + 1) - 1)) = pBodiesLock[ii]->sysMatInd1;
        }

        // Extract submatrices and invert system matrix
        Logger::debug("Inverting system matrix...");
        if (numBodiesFree > 0)
        {
            pSystemMatrixFF = new arma::mat(6 * numBodiesFree, 6 * numBodiesFree, arma::fill::zeros);
            *pSystemMatrixFF = (*pSystemMatrix)(sysMatIndFree, sysMatIndFree);
            if (numBodiesLock > 0)
            {
                pSystemMatrixFL = new arma::mat(6 * numBodiesFree, 6 * numBodiesLock, arma::fill::zeros);
                *pSystemMatrixFL = (*pSystemMatrix)(sysMatIndFree, sysMatIndLock);
            }
        }
        *pSystemMatrixInv = arma::solve(*pSystemMatrix, eye(size(*pSystemMatrix)));
        if (pSystemMatrixInv->has_nan() || pSystemMatrixInv->has_inf())
        {
            Logger::warning("System matrix inverse has NaN or Inf entries!");
        }
        if (numBodiesFree > 0)
        {
            pSystemMatrixFFInv = new arma::mat(6 * numBodiesFree, 6 * numBodiesFree, arma::fill::zeros);
            *pSystemMatrixFFInv = arma::solve(*pSystemMatrixFF, eye(size(*pSystemMatrixFF)));
            if (pSystemMatrixFFInv->has_nan() || pSystemMatrixFFInv->has_inf())
            {
                Logger::warning("Free-body system matrix inverse has NaN or Inf entries!");
            }
        }
        Logger::debug("System matrix inverted.");

        // Initialize velocity buffers, checking if it is necessary to increase the buffer size
        int time_buffer_size = this->timeBufferSize;
        for (int ii = 0; ii < this->numBodies; ii++)
        {
            if (time_buffer_size < 10 * this->pBodies[ii]->pHydro->GetNumPointsIrf())
            {
                time_buffer_size = 10 * this->pBodies[ii]->pHydro->GetNumPointsIrf();
            }
        }
        for (int ii = 0; ii < this->numBodies; ii++)
        {
            this->pBodies[ii]->velBufferSize = time_buffer_size;
            this->pBodies[ii]->velBuffer = arma::zeros(6, time_buffer_size);
        }
        this->timeBufferSize = time_buffer_size;
        this->timeBuffer = arma::zeros(1, time_buffer_size);

        // Free memory
        delete[] numBodiesPerHDB;

        Logger::info("--> Bodies Properties Read");
    }
}

void Simulation::ReadBodiesYAML()
{
    Logger::info("--> Reading Bodies Properties");

    if (!yamlRoot["bodies"])
    {
        Logger::warning("'bodies' section not found in YAML. Setting numBodies = 0.");
        this->numBodies = 0;
    }
    else
    {
        this->numBodies = (int)yamlRoot["bodies"].size();
    }

    Logger::debug("Number of bodies: " + std::to_string(this->numBodies));

    if (numBodies > 0)
    {
        int body_count = 0;
        int diff_count = 0;
        int hydro_database_count = 0;
        std::string hydro_databases_name[300];
        int pos_database = 0;

        pBodies = new Body*[numBodies];
        for (int ii = 0; ii < numBodies; ii++)
        {
            YAML::Node bodyNode = yamlRoot["bodies"][ii];
            std::string body_type = bodyNode["type"].as<std::string>();

            Logger::debug("  body " + std::to_string(ii + 1));
            if (body_type.compare("RAD_DIFF") == 0)
            {
                pBodies[ii] = new Body(ii, this);
                pBodies[ii]->ReadPropertiesYAML(bodyNode);
                pBodies[ii]->OpenOutputFiles(outputFolderPath);
            }
            else
            {
                std::stringstream ss;
                ss << "Error in YAML bodies section: Unexpected body type '" << body_type << "' for body " << ii
                   << "\n";
                throw ValueError(ss.str());
            }
        }

        Logger::debug("--> Bodies loaded");

        // Loop over bodies to get the number of hydrodynamic databases
        hydro_databases_name[hydro_database_count] = pBodies[0]->hydroDatabaseName;
        hydro_database_count++;
        for (int ii = 1; ii < this->numBodies; ii++)
        {
            diff_count = 0;
            for (int jj = 0; jj < hydro_database_count; jj++)
            {
                if (pBodies[ii]->hydroDatabaseName.compare(hydro_databases_name[jj]) != 0)
                {
                    diff_count++;
                }
            }
            if (diff_count == hydro_database_count)
            {
                hydro_databases_name[hydro_database_count] = pBodies[ii]->hydroDatabaseName;
                hydro_database_count++;
            }
        }

        Logger::debug("--> HDB files:");
        for (int ii = 0; ii < hydro_database_count; ii++)
        {
            Logger::debug("    " + hydro_databases_name[ii]);
        }

        // Pre-determine the number of bodies in each HDB file
        int* numBodiesPerHDB = new int[hydro_database_count];
        for (int ii = 0; ii < hydro_database_count; ii++)
        {
            std::string hdb_path = JoinPath(inputFolderPath, hydro_databases_name[ii]);
            numBodiesPerHDB[ii] = HydroDatabase::GetNumBodiesFromFile(hdb_path);
            Logger::debug("    HDB: " + hydro_databases_name[ii] + " — " + std::to_string(numBodiesPerHDB[ii]) + " bodies");
        }

        // Arrange all the bodies by database (same logic as ASCII)
        Body** pBodiesSort = new Body*[numBodies];
        int* pBody_found = new int[numBodies];
        for (int ii = 0; ii < numBodies; ii++)
        {
            pBody_found[ii] = 0;
        }
        for (int ii = 0; ii < hydro_database_count; ii++)
        {
            if (numBodiesPerHDB[ii] == 1)
            {
                for (int jj = 0; jj < numBodies; jj++)
                {
                    if ((hydro_databases_name[ii].compare(pBodies[jj]->hydroDatabaseName) == 0) &&
                        (pBody_found[jj] == 0))
                    {
                        pBody_found[jj] = 1;
                        pBodiesSort[body_count] = pBodies[jj];
                        body_count++;
                    }
                }
            }
            else
            {
                int temp_nB_hdb = numBodiesPerHDB[ii];
                int temp_nB_found = 0;
                arma::uvec temp_hdb_ind(temp_nB_hdb);
                for (int jj = 0; jj < numBodies; jj++)
                {
                    if ((hydro_databases_name[ii].compare(pBodies[jj]->hydroDatabaseName) == 0) &&
                        (pBody_found[jj] == 0))
                    {
                        pBody_found[jj] = 1;
                        pBodiesSort[body_count + pBodies[jj]->hydroDatabaseIndex] = pBodies[jj];
                        if (temp_nB_found < temp_nB_hdb)
                        {
                            temp_hdb_ind(temp_nB_found) = pBodies[jj]->hydroDatabaseIndex;
                        }
                        temp_nB_found++;
                    }
                }
                if (temp_nB_found != temp_nB_hdb)
                {
                    std::stringstream ss;
                    ss << "ERROR: Not all bodies in multi-body HDB " << hydro_databases_name[ii]
                       << " are used. Expected " << temp_nB_hdb << " but found " << temp_nB_found << ".\n";
                    throw ValueError(ss.str());
                }
                arma::uvec temp_hdb_ind_unique = arma::unique(temp_hdb_ind);
                if ((int)temp_hdb_ind_unique.n_elem != temp_nB_hdb)
                {
                    std::stringstream ss;
                    ss << "ERROR: Repeated body indices in multi-body HDB: " << hydro_databases_name[ii] << "\n";
                    throw ValueError(ss.str());
                }
                body_count += temp_nB_hdb;
            }
        }
        delete[] pBody_found;
        for (int ii = 0; ii < numBodies; ii++)
        {
            pBodies[ii] = pBodiesSort[ii];
        }
        delete[] pBodiesSort;

        // Checking free bodies
        Logger::debug("Checking for free bodies (YAML)...");
        for (int ii = 0; ii < numBodies; ii++)
        {
            if (pBodies[ii]->flag_blocked > 0)
            {
                numBodiesLock++;
            }
            else
            {
                numBodiesFree++;
            }
        }
        pBodiesLock = new Body*[numBodiesLock];
        pBodiesFree = new Body*[numBodiesFree];
        int indL = 0;
        int indF = 0;
        for (int ii = 0; ii < numBodies; ii++)
        {
            if (pBodies[ii]->flag_blocked > 0)
            {
                pBodiesLock[indL] = pBodies[ii];
                indL++;
            }
            else
            {
                pBodiesFree[indF] = pBodies[ii];
                indF++;
            }
        }

        // Set hydrodynamic database to each body
        std::string hydro_file_path;
        for (int ii = 0; ii < this->numBodies; ii++)
        {
            pos_database = 0;
            while (true)
            {
                if (this->pBodies[ii]->hydroDatabaseName.compare(hydro_databases_name[pos_database]) == 0)
                {
                    break;
                }
                pos_database++;
            }

            int nBodiesInHDB = numBodiesPerHDB[pos_database];
            if (nBodiesInHDB == 1)
            {
                Body** pSingleBodyArray = new Body*[1];
                pSingleBodyArray[0] = pBodies[ii];
                this->pBodies[ii]->LoadHydrodynamicDatabase(pSingleBodyArray, 0);
            }
            else
            {
                Body** pMultiBodyArray = new Body*[nBodiesInHDB];
                for (int jj = 0; jj < nBodiesInHDB; jj++)
                {
                    pMultiBodyArray[jj] = pBodies[ii - pBodies[ii]->hydroDatabaseIndex + jj];
                }
                this->pBodies[ii]->LoadHydrodynamicDatabase(pMultiBodyArray, pBodies[ii]->hydroDatabaseIndex);
            }
        }

        // Fill System Matrix
        Logger::debug("Filling system matrix...");
        arma::span a1;
        arma::span a2;
        this->pSystemMatrix = new arma::mat(6 * this->numBodies, 6 * this->numBodies, arma::fill::zeros);
        this->pSystemMatrixInv = new arma::mat(6 * this->numBodies, 6 * this->numBodies, arma::fill::zeros);
        for (int ii = 0; ii < this->numBodies; ii++)
        {
            a1 = arma::span(6 * ii, 6 * ii + 5);
            int temp_nB_hdb = pBodies[ii]->pHydro->GetNumBodies();
            if (temp_nB_hdb == 1)
            {
                a2 = a1;
            }
            else
            {
                int temp_indHDB = pBodies[ii]->hydroDatabaseIndex;
                a2 = arma::span(6 * (ii - temp_indHDB), 6 * (ii - temp_indHDB + temp_nB_hdb) - 1);
            }
            Logger::debug("  body " + std::to_string(ii + 1) + " added");
            (*pSystemMatrix)(a1, a2) += pBodies[ii]->pHydro->GetTotalMass();
            pBodies[ii]->sysMatSpan1 = a1;
            pBodies[ii]->sysMatSpan2 = a2;
            pBodies[ii]->sysMatInd1 = arma::regspace<arma::uvec>(6 * ii, 6 * ii + 5);
        }

        // Take free dofs index vectors from the system matrix
        sysMatIndFree = arma::zeros<arma::uvec>(6 * numBodiesFree);
        for (int ii = 0; ii < this->numBodiesFree; ii++)
        {
            sysMatIndFree(arma::span(6 * ii, 6 * (ii + 1) - 1)) = pBodiesFree[ii]->sysMatInd1;
        }
        sysMatIndLock = arma::zeros<arma::uvec>(6 * numBodiesLock);
        for (int ii = 0; ii < this->numBodiesLock; ii++)
        {
            sysMatIndLock(arma::span(6 * ii, 6 * (ii + 1) - 1)) = pBodiesLock[ii]->sysMatInd1;
        }

        // Extract submatrices and invert system matrix
        Logger::debug("Inverting system matrix...");
        if (numBodiesFree > 0)
        {
            pSystemMatrixFF = new arma::mat(6 * numBodiesFree, 6 * numBodiesFree, arma::fill::zeros);
            *pSystemMatrixFF = (*pSystemMatrix)(sysMatIndFree, sysMatIndFree);
            if (numBodiesLock > 0)
            {
                pSystemMatrixFL = new arma::mat(6 * numBodiesFree, 6 * numBodiesLock, arma::fill::zeros);
                *pSystemMatrixFL = (*pSystemMatrix)(sysMatIndFree, sysMatIndLock);
            }
        }
        *pSystemMatrixInv = arma::solve(*pSystemMatrix, eye(size(*pSystemMatrix)));
        if (pSystemMatrixInv->has_nan() || pSystemMatrixInv->has_inf())
        {
            Logger::warning("System matrix inverse has NaN or Inf entries!");
        }
        if (numBodiesFree > 0)
        {
            pSystemMatrixFFInv = new arma::mat(6 * numBodiesFree, 6 * numBodiesFree, arma::fill::zeros);
            *pSystemMatrixFFInv = arma::solve(*pSystemMatrixFF, eye(size(*pSystemMatrixFF)));
            if (pSystemMatrixFFInv->has_nan() || pSystemMatrixFFInv->has_inf())
            {
                Logger::warning("Free-body system matrix inverse has NaN or Inf entries!");
            }
        }
        Logger::debug("System matrix inverted.");

        // Initialize velocity buffers
        int time_buffer_size = this->timeBufferSize;
        for (int ii = 0; ii < this->numBodies; ii++)
        {
            if (time_buffer_size < 10 * this->pBodies[ii]->pHydro->GetNumPointsIrf())
            {
                time_buffer_size = 10 * this->pBodies[ii]->pHydro->GetNumPointsIrf();
            }
        }
        for (int ii = 0; ii < this->numBodies; ii++)
        {
            this->pBodies[ii]->velBufferSize = time_buffer_size;
            this->pBodies[ii]->velBuffer = arma::zeros(6, time_buffer_size);
        }
        this->timeBufferSize = time_buffer_size;
        this->timeBuffer = arma::zeros(1, time_buffer_size);

        delete[] numBodiesPerHDB;

        Logger::info("--> Bodies Properties Read");
    }
}

void Simulation::ReadLines()
{
    (this->*pReadLines)();
}

void Simulation::ReadLinesASCII()
{
    Logger::info("--> Reading Lines Properties");

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataLines.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        Logger::warning("dataLines.dat not found. Setting numLines = 0.");
        numLines = 0;
        pLines = new Line*[0];
        return;
    }

    // Read line types
    fscanf(file_pointer, "%d %[^\n]\n", &numLineTypes, bufferLine);
    pLineTypes = new LineType*[numLineTypes];
    for (int ii = 0; ii < numLineTypes; ii++)
    {
        pLineTypes[ii] = new LineType();
        pLineTypes[ii]->ReadPropertiesASCII(file_pointer);
    }

    // Read line instances
    fscanf(file_pointer, "%d %[^\n]\n", &numLines, bufferLine);
    pLines = new Line*[numLines];

    for (int ii = 0; ii < numLines; ii++)
    {
        pLines[ii] = new Line(ii, gravity, waterDensity, waterDepth);
        try
        {
            pLines[ii]->ReadPropertiesASCII(file_pointer, pLineTypes, numLineTypes);
            pLines[ii]->OpenOutputFiles(outputFolderPath, outputFormat);
        }
        catch (int e)
        {
            if (e == 0)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is under the floor level.");
            if (e == 1)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " touches the seafloor and it shouldn't.");
            if (e == 2)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is not tense and laying on the seafloor. It should be pretensed.");
            if (e == 3)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is not tense and vertical. It should be pretensed.");
            if (e == 4)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is not tense and it should.");
            if (e == 5)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " initial shape can't be computed with QS method.");
            if (e == 6)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " touches the seafloor although none of its ends are there.");
        }
    }

    // Close the file
    fclose(file_pointer);

    Logger::info("--> Lines Properties Read");
}

void Simulation::ReadLinesYAML()
{
    Logger::info("--> Reading Lines Properties");

    // Read line types
    if (yamlRoot["line_types"])
    {
        numLineTypes = (int)yamlRoot["line_types"].size();
        pLineTypes = new LineType*[numLineTypes];
        for (int ii = 0; ii < numLineTypes; ii++)
        {
            pLineTypes[ii] = new LineType();
            pLineTypes[ii]->ReadPropertiesYAML(yamlRoot["line_types"][ii]);
        }
    }
    else
    {
        numLineTypes = 0;
        pLineTypes = new LineType*[0];
    }

    if (!yamlRoot["lines"])
    {
        Logger::warning("'lines' section not found in YAML. Setting numLines = 0.");
        numLines = 0;
        pLines = new Line*[0];
        return;
    }

    numLines = (int)yamlRoot["lines"].size();
    pLines = new Line*[numLines];

    for (int ii = 0; ii < numLines; ii++)
    {
        pLines[ii] = new Line(ii, gravity, waterDensity, waterDepth);
        try
        {
            pLines[ii]->ReadPropertiesYAML(yamlRoot["lines"][ii], pLineTypes, numLineTypes);
            pLines[ii]->OpenOutputFiles(outputFolderPath, outputFormat);
        }
        catch (int e)
        {
            if (e == 0)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is under the floor level.");
            if (e == 1)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " touches the seafloor and it shouldn't.");
            if (e == 2)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is not tense and laying on the seafloor. It should be pretensed.");
            if (e == 3)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is not tense and vertical. It should be pretensed.");
            if (e == 4)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is not tense and it should.");
            if (e == 5)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " initial shape can't be computed with QS method.");
            if (e == 6)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " touches the seafloor although none of its ends are there.");
        }
    }

    Logger::info("--> Lines Properties Read");
}

void Simulation::ReadSinking()
{
    (this->*pReadSinking)();
}

void Simulation::ReadSinkingASCII()
{
    Logger::info("--> Reading Sinking Properties");

    // Declare local variables
    char bufferLine[1000];
    int sinkingBodyIndex;

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataSinking.dat");
    FILE* pFile = fopen(file_path.c_str(), "r");
    if (pFile == NULL)
    {
        Logger::warning("dataSinking.dat not found. Setting numSinking = 0.");
        numSinking = 0;
        pSinking = new Sinking*[0];
        return;
    }

    // Read total number of sinking bodies to read
    fscanf(pFile, "%d %[^\n]\n", &numSinking, bufferLine);
    fscanf(pFile, "%d %[^\n]\n", &numSinking, bufferLine);

    // Read sinking bodies if there are any
    if (numSinking > 0)
    {
        // Check if there are more than one sinking body
        if (numSinking > 1)
        {
            std::stringstream ss;
            ss << "Multiple bodies sinking is not implemented yet." << std::endl;
            throw IOError(ss.str());
        }

        // Read all sinking bodies
        pSinking = new Sinking*[numSinking];
        for (int ii = 0; ii < numSinking; ii++)
        {
            // Discard header lines
            for (int ii = 0; ii < 3; ii++)
            {
                fgets(bufferLine, sizeof(bufferLine), pFile);
            }

            // Get sinking body id
            fscanf(pFile, "%d %[^\n]\n", &sinkingBodyIndex, bufferLine);

            // Initiallice Sinking object
            pSinking[ii] = new Sinking(sinkingBodyIndex, this);

            // Read sinking body properties
            pSinking[ii]->ReadPropertiesASCII(pFile, inputFolderPath);

            // Open output files
            pSinking[ii]->OpenOutputFiles(outputFolderPath);
        }
    }
    else
    {
        pSinking = new Sinking*[0];
    }

    // Close file
    fclose(pFile);

    Logger::info("--> Sinking Properties Read");
}

void Simulation::ReadSinkingYAML()
{
    Logger::info("--> Reading Sinking Properties");

    if (!yamlRoot["sinking"])
    {
        Logger::warning("'sinking' section not found in YAML. Setting numSinking = 0.");
        numSinking = 0;
        pSinking = new Sinking*[0];
        return;
    }

    numSinking = (int)yamlRoot["sinking"].size();

    if (numSinking > 0)
    {
        if (numSinking > 1)
        {
            std::stringstream ss;
            ss << "Multiple bodies sinking is not implemented yet." << std::endl;
            throw IOError(ss.str());
        }

        pSinking = new Sinking*[numSinking];
        for (int ii = 0; ii < numSinking; ii++)
        {
            YAML::Node sinkNode = yamlRoot["sinking"][ii];
            int sinkingBodyIndex = sinkNode["body_index"].as<int>();
            pSinking[ii] = new Sinking(sinkingBodyIndex, this);
            pSinking[ii]->ReadPropertiesYAML(sinkNode, inputFolderPath);
            pSinking[ii]->OpenOutputFiles(outputFolderPath);
        }
    }
    else
    {
        pSinking = new Sinking*[0];
    }

    Logger::info("--> Sinking Properties Read");
}

void Simulation::ReadSprings()
{
    (this->*pReadSprings)();
}

void Simulation::ReadSpringsASCII()
{
    Logger::info("--> Reading Spring Properties");

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataSprings.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        Logger::warning("dataSprings.dat not found. Setting numSprings = 0.");
        numSprings = 0;
        pSprings = new Spring*[0];
        return;
    }

    // Read spring types
    fscanf(file_pointer, "%d %[^\n]\n", &numSpringTypes, bufferLine);
    pSpringTypes = new SpringType*[numSpringTypes];
    for (int ii = 0; ii < numSpringTypes; ii++)
    {
        pSpringTypes[ii] = new SpringType();
        pSpringTypes[ii]->ReadPropertiesASCII(file_pointer);
    }

    // Read spring instances
    fscanf(file_pointer, "%d %[^\n]\n", &numSprings, bufferLine);
    pSprings = new Spring*[numSprings];
    for (int ii = 0; ii < numSprings; ii++)
    {
        pSprings[ii] = new Spring(ii);
        pSprings[ii]->ReadPropertiesASCII(file_pointer, pSpringTypes, numSpringTypes);
    }

    // Close the file
    fclose(file_pointer);

    Logger::info("--> Spring Properties Read");
}

void Simulation::ReadSpringsYAML()
{
    Logger::info("--> Reading Spring Properties");

    // Read spring types
    if (yamlRoot["spring_types"])
    {
        numSpringTypes = (int)yamlRoot["spring_types"].size();
        pSpringTypes = new SpringType*[numSpringTypes];
        for (int ii = 0; ii < numSpringTypes; ii++)
        {
            pSpringTypes[ii] = new SpringType();
            pSpringTypes[ii]->ReadPropertiesYAML(yamlRoot["spring_types"][ii]);
        }
    }
    else
    {
        numSpringTypes = 0;
        pSpringTypes = new SpringType*[0];
    }

    if (!yamlRoot["springs"])
    {
        Logger::warning("'springs' section not found in YAML. Setting numSprings = 0.");
        numSprings = 0;
        pSprings = new Spring*[0];
        return;
    }

    numSprings = (int)yamlRoot["springs"].size();
    pSprings = new Spring*[numSprings];

    for (int ii = 0; ii < numSprings; ii++)
    {
        pSprings[ii] = new Spring(ii);
        pSprings[ii]->ReadPropertiesYAML(yamlRoot["springs"][ii], pSpringTypes, numSpringTypes);
    }

    Logger::info("--> Spring Properties Read");
}

void Simulation::ReadProperties()
{
    (this->*pReadProperties)();
}

void Simulation::ReadPropertiesASCII()
{
    Logger::info("--> Reading Simulation Properties");

    // Declare local variables
    char bufferLine[1000];
    int dummyBool;

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataProblem.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: dataProblem.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
    }

    // Read data
    fscanf(file_pointer, "%lf %[^\n]\n", &gravity, bufferLine);      // Gravity acceleration [m/s^2]
    fscanf(file_pointer, "%lf %[^\n]\n", &waterDensity, bufferLine); // Water density [kg/m^3]
    fscanf(file_pointer, "%lf %[^\n]\n", &airAtmPresDensity,
           bufferLine);                                            // Air density at atmospheric pressure [kg/m^3]
    fscanf(file_pointer, "%lf %[^\n]\n", &airAtmPres, bufferLine); // Atmospheric pressure [Pa]
    fscanf(file_pointer, "%lf %[^\n]\n", &airAdiabaticDilation, bufferLine); // Air adiabatic dilation [-]
    fscanf(file_pointer, "%lf %[^\n]\n", &waterDepth, bufferLine);           // Seabed vertical coordinate [m]
    fscanf(file_pointer, "%lf %[^\n]\n", &writeTimeStep, bufferLine);        // Output time step [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &maxTimeStep, bufferLine); // Maximum time step for time integration [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &hydroTimeStep,
           bufferLine); // Time step for hydrodynamic forces computation [s]
    lastHydroTime = 0.0;
    lastHydroTime_old = -hydroTimeStep;
    lastHydroTime_old2 = -2.0 * hydroTimeStep;
    fscanf(file_pointer, "%lf %[^\n]\n", &fastTimeStep,
           bufferLine); // Time step for FAST wind turbines forces computation [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &fastControllerTimeStep,
           bufferLine); // Time step for FAST wind turbines controller update [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &timeIRF, bufferLine); // IRF time [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &sinkingTimeStep,
           bufferLine); // Time step for synking hydrodinamic data bases update [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &winchesContTimeStep, bufferLine); // Time step for winches controller [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &owcsContTimeStep, bufferLine);    // Time step for OWCS controller [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &simulationTime, bufferLine);      // Total time of simulation [s]
    fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    rotSimpFlag = dummyBool; // Flag to use simplification for rigid body rotation dynamics [0 No, 1 Yes]
    fscanf(file_pointer, "%d %[^\n]\n", &timeIntMethod,
           bufferLine);                                             // Temporal integration alforithm [1: BDF1, 2: BDFN]
    fscanf(file_pointer, "%d %[^\n]\n", &timeIntOrder, bufferLine); // Order for temporal integration
    fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    timeIntAdaptivity = dummyBool; // Time step adaptivity [0 No, 1 Yes]
    fscanf(file_pointer, "%d %[^\n]\n", &timeIntJacNumStepsMax,
           bufferLine); // Maximum number of steps for temporal integration Jacobian re-computation.
    fscanf(file_pointer, "%lf %[^\n]\n", &timeIntAbsTol, bufferLine); // Absolute tolerance for temporal integration.
    fscanf(file_pointer, "%lf %[^\n]\n", &timeIntRelTol, bufferLine); // Relative tolerance for temporal integration.
    fscanf(file_pointer, "%d %[^\n]\n", &maxIterStep,
           bufferLine); // Maximum number of iterations for one step of temporal integration.
    fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    readEquilibrium = dummyBool; // Read dataStaticIC.dat? [0 No, 1 Yes]
    fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    writeEquilibrium = dummyBool; // Write dataStaticIC.dat? [0 No, 1 Yes]
    fscanf(file_pointer, "%d %[^\n]\n", &flagStatic,
           bufferLine); // Mooring initial condition flag [0: Catenary, 1: Newton's]

    // Output format [0: txt, 1: csv] — optional, defaults to 0 (txt)
    if (fscanf(file_pointer, "%d %[^\n]\n", &outputFormat, bufferLine) != 1)
        outputFormat = 0;

    // Close file
    fclose(file_pointer);

    // Show inputs
    // TODO: Implement a logger with different levels of verbosity
    if (false)
    {
        std::cout << "Gravity: " << gravity << std::endl;
        std::cout << "Water Density: " << waterDensity << std::endl;
        std::cout << "Water Depth: " << waterDepth << std::endl;
        std::cout << "Max Time Step: " << maxTimeStep << std::endl;
        std::cout << "Simulation Time: " << simulationTime << std::endl;
        std::cout << "Time Integration Method: " << timeIntMethod << std::endl;
        std::cout << "Time Integration Absolute Tolerace: " << timeIntAbsTol << std::endl;
        std::cout << "Time Integration Relative Tolerace: " << timeIntRelTol << std::endl;
        std::cout << "Max Iterations per Step: " << maxIterStep << std::endl;
        std::cout << "Read Equilibrium: " << readEquilibrium << std::endl;
        std::cout << "Write Equilibrium: " << writeEquilibrium << std::endl;
        std::cout << "Static Equilibrium Method: " << flagStatic << std::endl;
    }

    Logger::info("--> Simulation Properties Read");
}

void Simulation::ReadPropertiesYAML()
{
    Logger::info("--> Reading Simulation Properties");

    // Open and parse YAML file
    std::string file_path = JoinPath(inputFolderPath, "dataProblem.yaml");
    try
    {
        yamlRoot = YAML::LoadFile(file_path);
    }
    catch (const YAML::Exception& e)
    {
        std::stringstream ss;
        ss << "Error parsing YAML file: " << file_path << "\n    " << e.what();
        throw IOError(ss.str());
    }

    if (!yamlRoot["problem"])
    {
        throw ValueError("YAML file missing required 'problem' section.");
    }

    YAML::Node prob = yamlRoot["problem"];
    gravity = prob["gravity"].as<double>();
    waterDensity = prob["water_density"].as<double>();
    airAtmPresDensity = prob["air_density"].as<double>();
    airAtmPres = prob["atmospheric_pressure"].as<double>();
    airAdiabaticDilation = prob["air_adiabatic_dilation"].as<double>();
    waterDepth = prob["water_depth"].as<double>();
    writeTimeStep = prob["write_time_step"].as<double>();
    maxTimeStep = prob["max_time_step"].as<double>();
    hydroTimeStep = prob["hydro_time_step"].as<double>();
    lastHydroTime = 0.0;
    lastHydroTime_old = -hydroTimeStep;
    lastHydroTime_old2 = -2.0 * hydroTimeStep;
    fastTimeStep = prob["fast_time_step"].as<double>();
    fastControllerTimeStep = prob["fast_controller_time_step"].as<double>();
    timeIRF = prob["irf_time"].as<double>();
    sinkingTimeStep = prob["sinking_time_step"].as<double>();
    winchesContTimeStep = prob["winches_controller_time_step"].as<double>();
    owcsContTimeStep = prob["owcs_controller_time_step"].as<double>();
    simulationTime = prob["simulation_time"].as<double>();
    rotSimpFlag = prob["rotation_simplification"].as<int>();
    timeIntMethod = prob["time_integration_method"].as<int>();
    timeIntOrder = prob["time_integration_order"].as<int>();
    timeIntAdaptivity = prob["time_step_adaptivity"].as<int>();
    timeIntJacNumStepsMax = prob["jacobian_recomputation_steps"].as<int>();
    timeIntAbsTol = prob["absolute_tolerance"].as<double>();
    timeIntRelTol = prob["relative_tolerance"].as<double>();
    maxIterStep = prob["max_iterations_per_step"].as<int>();
    readEquilibrium = prob["read_equilibrium"].as<int>();
    writeEquilibrium = prob["write_equilibrium"].as<int>();
    flagStatic = prob["mooring_initial_condition"].as<int>();

    // Output format — optional, defaults to "txt"
    if (prob["output_format"])
    {
        std::string fmt = prob["output_format"].as<std::string>();
        if (fmt == "csv")
            outputFormat = 1;
        else
            outputFormat = 0;
        Logger::debug("Output format: " + fmt + " (" + std::to_string(outputFormat) + ")");
    }
    else
    {
        outputFormat = 0;
        Logger::debug("Output format: txt (default)");
    }

    Logger::info("--> Simulation Properties Read");
}

void Simulation::ReadWaves()
{
    (this->*pReadWaves)();
}

void Simulation::ReadWavesASCII()
{
    Logger::info("--> Reading Waves");

    // Declare local variables
    char bufferLine[1000];
    char wave_type[1000];
    double H, T, D, rampTime;

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataWaves.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        Logger::warning("dataWaves.dat not found. Using zero-height wave.");
        pWave = new RegularWave(this, 0.0, 1.0, 0.0, 0.0);
        return;
    }

    // Ignore header lines
    for (int ii = 0; ii < 3; ii++)
    {
        fgets(bufferLine, sizeof(bufferLine), file_pointer);
    }

    // Read wave type and main parameters
    fgets(wave_type, sizeof(wave_type), file_pointer);
    fscanf(file_pointer, "%lf %[^\n]\n", &H, bufferLine);
    fscanf(file_pointer, "%lf %[^\n]\n", &T, bufferLine);
    fscanf(file_pointer, "%lf %[^\n]\n", &D, bufferLine);
    fscanf(file_pointer, "%lf %[^\n]\n", &rampTime, bufferLine);

    // Create wave object and read additional parameters if necessary
    if (strncmp(wave_type, "REG", 3) == 0)
    {
        pWave = new RegularWave(this, H, T, D, rampTime);
    }
    else if (strncmp(wave_type, "IRR", 3) == 0)
    {
        pWave = new IrregularWave(this, H, T, D, rampTime);
        for (int ii = 0; ii < 3; ii++)
        {
            fgets(bufferLine, sizeof(bufferLine), file_pointer);
        }
        fscanf(file_pointer, "%d %[^\n]\n", &pWave->specType_flag, bufferLine);
        fscanf(file_pointer, "%d %[^\n]\n", &pWave->piecewise_flag, bufferLine);
        fgets(bufferLine, sizeof(bufferLine), file_pointer);
        fscanf(file_pointer, "%lf %[^\n]\n", &pWave->gamma, bufferLine);
        fscanf(file_pointer, "%lf %[^\n]\n", &pWave->s, bufferLine);
        fscanf(file_pointer, "%lf %[^\n]\n", &pWave->dtheta, bufferLine);
        fscanf(file_pointer, "%lf %[^\n]\n", &pWave->factor, bufferLine);
        fscanf(file_pointer, "%d %[^\n]\n", &pWave->readPhases_flag, bufferLine);
        fscanf(file_pointer, "%lf %[^\n]\n", &pWave->rel_tol, bufferLine);
        fscanf(file_pointer, "%lf %[^\n]\n", &pWave->dt, bufferLine);
        char cWavePhasesFileName[1000];
        fscanf(file_pointer, "%s %[^\n]\n", cWavePhasesFileName, bufferLine);
        pWave->wavePhasesFileName = cWavePhasesFileName;
        pWave->filePhases_path = JoinPath(inputFolderPath, pWave->wavePhasesFileName);
        // Ignore one line
        fgets(bufferLine, sizeof(bufferLine), file_pointer);
        char cWaveDatabaseName[1000];
        fscanf(file_pointer, "%s %[^\n]\n", cWaveDatabaseName, bufferLine);
        pWave->waveDatabaseName = cWaveDatabaseName;
        pWave->file_path = JoinPath(inputFolderPath, pWave->waveDatabaseName);
    }
    else
    {
        std::stringstream ss;
        ss << "Error while parsing file: dataWaves.dat; Unexpected wave type. \n";
        throw ValueError(ss.str());
    }

    // Close file
    fclose(file_pointer);

    // Show inputs
    // TODO: Implement a logger with different levels of verbosity
    if (false)
    {
        if (strncmp(wave_type, "REG", 3) == 0)
        {
            std::cout << "Wave type: Regular" << std::endl;
            std::cout << "Wave height: " << H << std::endl;
            std::cout << "Wave period: " << T << std::endl;
            std::cout << "Wave heading: " << D << std::endl;
        }
        else
        {
            std::cout << "Wave type: Irregular" << std::endl;
            std::cout << "Wave significant height: " << H << std::endl;
            std::cout << "Wave peak period: " << T << std::endl;
            std::cout << "Wave heading: " << D << std::endl;
            if (pWave->specType_flag == 1)
            {
                std::cout << "Wave peak enhacement factor: " << pWave->gamma << std::endl;
                std::cout << "Wave directional spreading: " << pWave->s << std::endl;
                std::cout << "Wave directional step: " << pWave->dtheta << std::endl;
                std::cout << "Relative tolerance for wave check: " << pWave->rel_tol << std::endl;
                std::cout << "Time step for wave check: " << pWave->dt << std::endl;
            }
            else
            {
                std::cout << "Wave base data file name: " << pWave->waveDatabaseName << std::endl;
            }
        }
    }

    Logger::info("--> Wave read");

    // Preprocess wave data if H>0
    if (H > 0.0)
    {
        Logger::info("-->  Preprocessing wave...");
        pWave->CheckBreakingWave();
        pWave->GetWaveSpectrum();
        if (strncmp(wave_type, "IRR", 3) == 0)
        {
            pWave->WriteOut(outputFolderPath);
        }
    }
    else
    {
        Logger::debug("-->  Zero-height wave: skipping preprocessing.");
        pWave->SetZeroHeight();
    }

    Logger::info("--> Wave preprocessed");
}

void Simulation::ReadWavesYAML()
{
    Logger::info("--> Reading Waves");

    if (!yamlRoot["waves"])
    {
        Logger::warning("'waves' section not found in YAML. Using zero-height wave.");
        pWave = new RegularWave(this, 0.0, 1.0, 0.0, 0.0);
        return;
    }

    YAML::Node wav = yamlRoot["waves"];
    std::string wave_type = wav["type"].as<std::string>();
    double H = wav["height"].as<double>();
    double T = wav["period"].as<double>();
    double D = wav["heading"].as<double>();
    double rampTime = wav["ramp_time"].as<double>();

    if (wave_type == "REG")
    {
        pWave = new RegularWave(this, H, T, D, rampTime);
    }
    else if (wave_type == "IRR")
    {
        pWave = new IrregularWave(this, H, T, D, rampTime);
        pWave->specType_flag = wav["spectrum_type"].as<int>();
        pWave->piecewise_flag = wav["piecewise_flag"].as<int>();
        pWave->gamma = wav["gamma"].as<double>();
        pWave->s = wav["spreading"].as<double>();
        pWave->dtheta = wav["dtheta"].as<double>();
        pWave->factor = wav["factor"].as<double>();
        pWave->readPhases_flag = wav["read_phases_flag"].as<int>();
        pWave->rel_tol = wav["relative_tolerance"].as<double>();
        pWave->dt = wav["dt"].as<double>();
        pWave->wavePhasesFileName = wav["phases_file"].as<std::string>();
        pWave->filePhases_path = JoinPath(inputFolderPath, pWave->wavePhasesFileName);
        pWave->waveDatabaseName = wav["database_file"].as<std::string>();
        pWave->file_path = JoinPath(inputFolderPath, pWave->waveDatabaseName);
    }
    else
    {
        std::stringstream ss;
        ss << "Error in YAML waves section: Unexpected wave type '" << wave_type << "'.\n";
        throw ValueError(ss.str());
    }

    Logger::info("--> Wave read");

    // Preprocess wave data if H>0
    if (H > 0.0)
    {
        Logger::info("-->  Preprocessing wave...");
        pWave->CheckBreakingWave();
        pWave->GetWaveSpectrum();
        if (wave_type == "IRR")
        {
            pWave->WriteOut(outputFolderPath);
        }
    }
    else
    {
        Logger::debug("-->  Zero-height wave: skipping preprocessing.");
        pWave->SetZeroHeight();
    }

    Logger::info("--> Wave preprocessed");
}

void Simulation::ReadWinches()
{
    (this->*pReadWinches)();
}

void Simulation::ReadWinchesASCII()
{
    Logger::info("--> Reading Winches Properties");

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataWinches.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        Logger::warning("dataWinches.dat not found. Setting numWinches = 0.");
        numWinches = 0;
        pWinches = new Winchie*[0];
        WinchesController = nullptr;
        useWinches = false;
        return;
    }

    // Read number of winches defined in the file
    fscanf(file_pointer, "%d %[^\n]\n", &numWinches, bufferLine);

    if ((numWinches == 0) && useWinches)
    {
        throw ValueError(
            "Use of winches is requested when loading BCPs but there is no winches specified in dataWinches.dat\n");
    }

    // Allocate a vector of pointers to Winch class objects
    pWinches = new Winchie*[numWinches];

    // Read Winches
    for (int ii = 0; ii < numWinches; ii++)
    {
        pWinches[ii] = new Winchie(ii);
        pWinches[ii]->ReadPropertiesASCII(file_pointer);
        pWinches[ii]->LineW = pLines[pWinches[ii]->nLine - 1];
    }

    // Close the file
    fclose(file_pointer);
    Logger::info("--> Winches Properties Read");

    // Read Winches Controller
    Logger::info("--> Reading Winches Controller Properties");
    file_path = JoinPath(inputFolderPath, "dataWinchesController.dat");
    file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: dataWinchesController.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
    }
    // Skip 3 header lines
    for (int ii = 0; ii < 3; ii++)
        fgets(bufferLine, sizeof(bufferLine), file_pointer);
    // Read controller type
    int controllerType;
    fscanf(file_pointer, "%d %[^\n]\n", &controllerType, bufferLine);
    // Create appropriate controller via factory
    WinchesController = WinchieController::Create(controllerType, numWinches, pWinches, this);
    WinchesController->ReadPropertiesASCII(file_pointer);
    fclose(file_pointer);
    // Open Winches Controller output files
    WinchesController->OpenOutputFiles(outputFolderPath);
    Logger::info("--> Winches Controller Properties Read");
}

void Simulation::ReadWinchesYAML()
{
    Logger::info("--> Reading Winches Properties");

    if (!yamlRoot["winches"])
    {
        Logger::warning("'winches' section not found in YAML. Setting numWinches = 0.");
        numWinches = 0;
        pWinches = new Winchie*[0];
        WinchesController = nullptr;
        useWinches = false;
        return;
    }

    YAML::Node winchesNode = yamlRoot["winches"];
    YAML::Node winchList = winchesNode["winches"];
    numWinches = winchList ? (int)winchList.size() : 0;

    if ((numWinches == 0) && useWinches)
    {
        throw ValueError("Use of winches is requested when loading BCPs but there are no winches specified in YAML\n");
    }

    pWinches = new Winchie*[numWinches];
    for (int ii = 0; ii < numWinches; ii++)
    {
        pWinches[ii] = new Winchie(ii);
        pWinches[ii]->ReadPropertiesYAML(winchList[ii]);
        pWinches[ii]->LineW = pLines[pWinches[ii]->nLine - 1];
    }
    Logger::info("--> Winches Properties Read");

    // Read Winches Controller
    Logger::info("--> Reading Winches Controller Properties");
    if (!winchesNode["controller"])
    {
        std::stringstream ss;
        ss << "YAML 'winches' section missing 'controller' subsection.";
        throw IOError(ss.str());
    }
    YAML::Node ctrlNode = winchesNode["controller"];
    int controllerType = ctrlNode["type"].as<int>();
    WinchesController = WinchieController::Create(controllerType, numWinches, pWinches, this);
    WinchesController->ReadPropertiesYAML(ctrlNode);
    WinchesController->OpenOutputFiles(outputFolderPath);
    Logger::info("--> Winches Controller Properties Read");
}

void Simulation::ReadSeaFloor()
{
    (this->*pReadSeaFloor)();
}

void Simulation::ReadSeaFloorASCII()
{
    Logger::info("--> Reading SeaFloor");

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataSeaFloor.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        Logger::warning("dataSeaFloor.dat not found. Using flat sea floor.");
        numBathymetry = 0;
        numInclined = 0;
        numFlat = 1;
        numFloor = numBathymetry + numInclined + numFlat;
        pSeaFloor = new SeaFloor*[numFloor];
        pBathymetry = new Bathymetry*[0];
        pInclined = new Inclined*[0];
        pFlat = new Flat*[numFlat];
        pFlat[0] = new Flat(0);
        pSeaFloor[0] = pFlat[0];
        pSeaFloor[0]->seabedDepth = this->waterDepth;
        Logger::info("--> Floor Properties Read");
        return;
    }

    // Read data
    fscanf(file_pointer, "%d %[^\n]\n", &numBathymetry, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numInclined, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numFlat, bufferLine);
    numFloor = numBathymetry + numInclined + numFlat;

    // TODO: Implement a logger with different levels of verbosity
    // printf("Number of bathymetry: %d\n", numBathymetry);
    // printf("Number of slopes: %d\n", numInclined);
    // printf("Number of planes: %d\n", numFlat);

    // Allocate a vector of pointers to SeaFloor class objects
    pSeaFloor = new SeaFloor*[numFloor];
    int floor_count = 0;

    // Read all SeaFloor objects
    pBathymetry = new Bathymetry*[numBathymetry];
    for (int ii = 0; ii < numBathymetry; ii++)
    {
        pBathymetry[ii] = new Bathymetry(floor_count);
        pBathymetry[ii]->ReadPropertiesASCII(file_pointer, inputFolderPath);
        pSeaFloor[floor_count] = pBathymetry[ii];
        floor_count++;
    }
    pInclined = new Inclined*[numInclined];
    for (int ii = 0; ii < numInclined; ii++)
    {
        pInclined[ii] = new Inclined(floor_count);
        pInclined[ii]->ReadPropertiesASCII(file_pointer);
        pSeaFloor[floor_count] = pInclined[ii];
        floor_count++;
    }
    pFlat = new Flat*[numFlat];
    for (int ii = 0; ii < numFlat; ii++)
    {
        pFlat[ii] = new Flat(floor_count);
        pFlat[ii]->ReadPropertiesASCII(file_pointer);
        pSeaFloor[floor_count] = pFlat[ii];
        floor_count++;
    }

    // Close file
    fclose(file_pointer);

    Logger::info("--> Floor Properties Read");
}

void Simulation::ReadSeaFloorYAML()
{
    Logger::info("--> Reading SeaFloor");

    if (!yamlRoot["seafloor"])
    {
        Logger::warning("'seafloor' section not found in YAML. Using flat sea floor.");
        numBathymetry = 0;
        numInclined = 0;
        numFlat = 1;
        numFloor = 1;
        pSeaFloor = new SeaFloor*[numFloor];
        pBathymetry = new Bathymetry*[0];
        pInclined = new Inclined*[0];
        pFlat = new Flat*[numFlat];
        pFlat[0] = new Flat(0);
        pSeaFloor[0] = pFlat[0];
        pSeaFloor[0]->seabedDepth = this->waterDepth;
        Logger::info("--> Floor Properties Read");
        return;
    }

    YAML::Node sfNode = yamlRoot["seafloor"];
    numBathymetry = sfNode["bathymetry"] ? (int)sfNode["bathymetry"].size() : 0;
    numInclined = sfNode["inclined"] ? (int)sfNode["inclined"].size() : 0;
    numFlat = sfNode["flat"] ? (int)sfNode["flat"].size() : 0;
    numFloor = numBathymetry + numInclined + numFlat;

    pSeaFloor = new SeaFloor*[numFloor];
    int floor_count = 0;

    pBathymetry = new Bathymetry*[numBathymetry];
    for (int ii = 0; ii < numBathymetry; ii++)
    {
        pBathymetry[ii] = new Bathymetry(floor_count);
        pBathymetry[ii]->ReadPropertiesYAML(sfNode["bathymetry"][ii], inputFolderPath);
        pSeaFloor[floor_count] = pBathymetry[ii];
        floor_count++;
    }
    pInclined = new Inclined*[numInclined];
    for (int ii = 0; ii < numInclined; ii++)
    {
        pInclined[ii] = new Inclined(floor_count);
        pInclined[ii]->ReadPropertiesYAML(sfNode["inclined"][ii]);
        pSeaFloor[floor_count] = pInclined[ii];
        floor_count++;
    }
    pFlat = new Flat*[numFlat];
    for (int ii = 0; ii < numFlat; ii++)
    {
        pFlat[ii] = new Flat(floor_count);
        pFlat[ii]->ReadPropertiesYAML(sfNode["flat"][ii]);
        pSeaFloor[floor_count] = pFlat[ii];
        floor_count++;
    }

    Logger::info("--> Floor Properties Read");
}

void Simulation::ReadWindTurbines(void)
{
    (this->*pReadWindTurbines)();
}

void Simulation::ReadWindTurbinesASCII(void)
{
    Logger::info("--> Reading Wind Turbines Properties");

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataWindTurbines.dat");
    FILE* file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        Logger::warning("dataWindTurbines.dat not found. Setting numWindTurbines = 0.");
        numWindTurbines = 0;
        pWindTurbines = new WindTurbine*[0];
        return;
    }

    // Read number of Wind Turbines defined in the file
    fscanf(file_pointer, "%d %[^\n]\n", &numWindTurbines, bufferLine);

    // Allocate a vector of pointers to WindTurbine class objects
    pWindTurbines = new WindTurbine*[numWindTurbines];
    for (int ii = 0; ii < numWindTurbines; ii++)
    {
        pWindTurbines[ii] = new WindTurbine(ii, this);
        pWindTurbines[ii]->ReadPropertiesASCII(file_pointer);
    }

    // Close the file
    fclose(file_pointer);

    Logger::info("--> Wind Turbines Properties Read");
}

void Simulation::ReadWindTurbinesYAML(void)
{
    Logger::info("--> Reading Wind Turbines Properties");

    if (!yamlRoot["wind_turbines"])
    {
        Logger::warning("'wind_turbines' section not found in YAML. Setting numWindTurbines = 0.");
        numWindTurbines = 0;
        pWindTurbines = new WindTurbine*[0];
        return;
    }

    numWindTurbines = (int)yamlRoot["wind_turbines"].size();
    pWindTurbines = new WindTurbine*[numWindTurbines];
    for (int ii = 0; ii < numWindTurbines; ii++)
    {
        pWindTurbines[ii] = new WindTurbine(ii, this);
        pWindTurbines[ii]->ReadPropertiesYAML(yamlRoot["wind_turbines"][ii]);
    }

    Logger::info("--> Wind Turbines Properties Read");
}

void Simulation::ReadOWCs(void)
{
    (this->*pReadOWCs)();
}

void Simulation::ReadOWCsASCII(void)
{
    Logger::info("--> Reading OWC Turbines Properties");
    // Declare local variables
    char bufferLineT[1000];
    // Open file
    std::string file_pathT = JoinPath(inputFolderPath, "dataOWCTurbines.dat");
    FILE* pFileT = fopen(file_pathT.c_str(), "r");
    if (pFileT == NULL)
    {
        Logger::warning("dataOWCTurbines.dat not found. Setting numOWCTurbines = 0.");
        numOWCTurbines = 0;
        pOWCTurbines = new OWCTurbineType*[0];
    }
    else
    {
        // Read total number of owcs
        fscanf(pFileT, "%d %[^\n]\n", &numOWCTurbines, bufferLineT);
        if (numOWCTurbines > 0)
        {
            // Initiallice OWC array
            pOWCTurbines = new OWCTurbineType*[numOWCTurbines];
            for (int ii = 0; ii < numOWCTurbines; ii++)
            {
                pOWCTurbines[ii] = new OWCTurbineType(ii, this);
                pOWCTurbines[ii]->Initialize(pFileT);
            }
        }
        else
        {
            pOWCTurbines = new OWCTurbineType*[0];
        }
        // Close file
        fclose(pFileT);
    }
    Logger::info("--> OWC Turbines Properties Read");

    Logger::info("--> Reading OWCs Properties");
    // Declare local variables
    char bufferLine[1000];
    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataOWCs.dat");
    FILE* pFile = fopen(file_path.c_str(), "r");
    if (pFile == NULL)
    {
        Logger::warning("dataOWCs.dat not found. Setting numOWCs = 0.");
        numOWCs = 0;
        pOWCs = new OWC*[0];
    }
    else
    {
        // Read total number of owcs
        fscanf(pFile, "%d %[^\n]\n", &numOWCs, bufferLine);
        if (numOWCs > 0)
        {
            // Initiallice OWC array
            pOWCs = new OWC*[numOWCs];
            for (int ii = 0; ii < numOWCs; ii++)
            {
                pOWCs[ii] = new OWC(ii, this);
                pOWCs[ii]->Initialize(pFile);
            }
        }
        else
        {
            pOWCs = new OWC*[0];
        }
        // Close file
        fclose(pFile);
    }
    Logger::info("--> OWCs Properties Read");
}

void Simulation::ReadOWCsYAML(void)
{
    Logger::info("--> Reading OWCs Properties");

    // Read OWC Turbine Types
    if (!yamlRoot["owcs"] || !yamlRoot["owcs"]["turbine_types"])
    {
        Logger::warning("'owcs.turbine_types' section not found in YAML. Setting numOWCTurbines = 0.");
        numOWCTurbines = 0;
        pOWCTurbines = new OWCTurbineType*[0];
    }
    else
    {
        YAML::Node turbTypes = yamlRoot["owcs"]["turbine_types"];
        numOWCTurbines = (int)turbTypes.size();
        if (numOWCTurbines > 0)
        {
            pOWCTurbines = new OWCTurbineType*[numOWCTurbines];
            for (int ii = 0; ii < numOWCTurbines; ii++)
            {
                pOWCTurbines[ii] = new OWCTurbineType(ii, this);
                pOWCTurbines[ii]->ReadPropertiesYAML(turbTypes[ii]);
                pOWCTurbines[ii]->Initialize(nullptr);
            }
        }
        else
        {
            pOWCTurbines = new OWCTurbineType*[0];
        }
    }
    Logger::info("--> OWC Turbines Properties Read");

    // Read OWC Chambers
    if (!yamlRoot["owcs"] || !yamlRoot["owcs"]["chambers"])
    {
        Logger::warning("'owcs.chambers' section not found in YAML. Setting numOWCs = 0.");
        numOWCs = 0;
        pOWCs = new OWC*[0];
    }
    else
    {
        YAML::Node chambers = yamlRoot["owcs"]["chambers"];
        numOWCs = (int)chambers.size();
        if (numOWCs > 0)
        {
            pOWCs = new OWC*[numOWCs];
            for (int ii = 0; ii < numOWCs; ii++)
            {
                pOWCs[ii] = new OWC(ii, this);
                pOWCs[ii]->ReadPropertiesYAML(chambers[ii]);
                pOWCs[ii]->Initialize(nullptr);
            }
        }
        else
        {
            pOWCs = new OWC*[0];
        }
    }
    Logger::info("--> OWCs Properties Read");
}

void Simulation::Run()
{
    Logger::info("--> Starting Simulation Run...");

    // Declare local variables
    time_t tstart, tend;
    double wallTime = 0.0;
    double wallTimeHydro = 0.0;
    double wallTimeFAST = 0.0;
    double wallTimeControllerFAST = 0.0;
    double wallTimeSinking = 0.0;
    double wallTimeControllerWinches = 0.0;
    double wallTimeControllerOWCs = 0.0;
    tstart = time(0);
    bool flag_debug_lines = false;

    { std::ostringstream _p; _p << std::fixed << std::setprecision(2) << "  t = " << wallTime << " / " << simulationTime << " s"; Logger::progress(_p.str()); }

    // Check if system has degrees of freedom
    if (numSystem > 0)
    {
        // Standard ODE solver-based time stepping for systems with DOFs
        // TODO: Implement a logger with different levels of verbosity
        //  std::cout<< "In Simulation::Run --> Starting temporal integration loop "<< std::endl;
        while (pTimeSolver->t < simulationTime - 2e-14)
        {
            // std::cout<< "In Simulation::Run --> step() "<< std::endl;
            pTimeSolver->step();

            // Print lines in first time step for debugging
            if (flag_debug_lines)
            {
                flag_debug_lines = false;
                for (int ii = 0; ii < numLines; ii = ii + 1)
                {
                    fprintf(pLines[ii]->pfile_debug, "s    x    y    z \n");
                    for (int jj = 0; jj < pLines[ii]->N; jj = jj + 1)
                        fprintf(pLines[ii]->pfile_debug, "%f    %f    %f    %f \n", pLines[ii]->s(jj, 0),
                                pLines[ii]->pos(jj, 0), pLines[ii]->pos(jj, 1), pLines[ii]->pos(jj, 2));
                    fclose(pLines[ii]->pfile_debug);
                }
            }

            if (pTimeSolver->t >= wallTime + writeTimeStep - 1e-12)
            {
                // std::cout<< "In Simulation::Run --> WriteOut() "<< std::endl;
                wallTime = writeTimeStep * round((wallTime + writeTimeStep) / writeTimeStep);
                { std::ostringstream _p; _p << std::fixed << std::setprecision(2) << "  t = " << wallTime << " / " << simulationTime << " s"; Logger::progress(_p.str()); }
                for (int ii = 0; ii < numLines; ii = ii + 1)
                    pLines[ii]->WriteOut(pTimeSolver->t);
                for (int ii = 0; ii < numBodies; ii = ii + 1)
                    pBodies[ii]->WriteOut(pTimeSolver->t);
                for (int ii = 0; ii < numOWCs; ii = ii + 1)
                    pOWCs[ii]->WriteOut(pTimeSolver->t);
                for (int ii = 0; ii < numSprings; ii = ii + 1)
                    pSprings[ii]->WriteOut(pTimeSolver->t);
            }

            for (int ii = 0; ii < numLines; ii = ii + 1)
            {
                if ((pLines[ii]->flag_visc == 1) && (pTimeSolver->t >= pLines[ii]->last_time + pLines[ii]->dt))
                {
                    pLines[ii]->update_buffer(pTimeSolver->t);
                    pLines[ii]->last_time = pTimeSolver->t;
                }
            }

            if (pTimeSolver->t >= wallTimeHydro + hydroTimeStep - 1e-12)
            {
                // std::cout<< "In Simulation::Run --> Computing hydrodynamic forces... "<< std::endl;
                wallTimeHydro = hydroTimeStep * round((wallTimeHydro + hydroTimeStep) / hydroTimeStep);
                if (numBodies > 0)
                {
                    UpdateSystem();
                }
                // Update Lagrange interpolation time stamps ONCE (outside body loop)
                lastHydroTime_old2 = lastHydroTime_old;
                lastHydroTime_old = lastHydroTime;
                lastHydroTime = pTimeSolver->t;
                for (int ii = 0; ii < numBodies; ii = ii + 1)
                {
                    pBodies[ii]->Fb_old2 = pBodies[ii]->Fb_old;
                    pBodies[ii]->Fb_old = pBodies[ii]->Fb;
                    pBodies[ii]->Fb = pBodies[ii]->pHydro->CalculateHydrodynamicForces(pTimeSolver->t) +
                                      pBodies[ii]->pHydro->CalculateHydrostaticForces(pTimeSolver->t);
                }
            }

            if (numWindTurbines > 0)
            {
                if (pTimeSolver->t >= wallTimeFAST + fastTimeStep)
                {
                    wallTimeFAST += fastTimeStep;
                    for (int ii = 0; ii < numWindTurbines; ii = ii + 1)
                    {
                        // std::cout<< "In Simulation::Run --> Computing forces for turbine " << ii << std::endl;
                        pWindTurbines[ii]->SetInputsFAST();
                        pWindTurbines[ii]->ComputeForces(wallTimeFAST);
                        pWindTurbines[ii]->WriteOut(wallTimeFAST);
                    }
                }
                if (pTimeSolver->t >= wallTimeControllerFAST + fastControllerTimeStep)
                {
                    wallTimeControllerFAST += fastControllerTimeStep;
                    for (int ii = 0; ii < numWindTurbines; ii = ii + 1)
                    {
                        // std::cout<< "In Simulation::Run --> Calling controller for turbine " << ii << std::endl;
                        pWindTurbines[ii]->SetInputsFAST();
                        pWindTurbines[ii]->ComputeControler(wallTimeControllerFAST);
                    }
                }
            }

            if (numOWCs > 0)
            {
                if (pTimeSolver->t >= wallTimeControllerOWCs + owcsContTimeStep)
                {
                    // std::cout << "In Simulation::Run --> Controlling OWCs... " << std::endl;
                    wallTimeControllerOWCs += owcsContTimeStep;
                    for (int ii = 0; ii < numOWCs; ii = ii + 1)
                    {
                        if (pOWCs[ii]->turbine_type > 0)
                        {
                            pOWCs[ii]->pOWCTurbine->ComputeGenTorque();
                        }
                    }
                }
            }

            if (numSinking > 0)
            {
                if (pTimeSolver->t >= wallTimeSinking + sinkingTimeStep)
                {
                    // std::cout << "In Simulation::Run --> Updating sinking... " << std::endl;
                    wallTimeSinking += sinkingTimeStep;
                    for (int ii = 0; ii < numSinking; ii = ii + 1)
                    {
                        pSinking[ii]->UpdateSinkingHydrodynamics(wallTime);
                    }
                    // std::cout << "In Simulation::Run --> UpdateSystemMatrix();" << std::endl;
                    UpdateSystemMatrix();
                    // std::cout << "    pSinking[ii]->WriteOut();" << std::endl;
                    for (int ii = 0; ii < numSinking; ii = ii + 1)
                        pSinking[ii]->WriteOut(wallTime);
                    // std::cout << "In Simulation::Run --> ... done updating sinking!" << std::endl;
                }
            }

            if (numWinches > 0)
            {
                if (pTimeSolver->t >= wallTimeControllerWinches + winchesContTimeStep)
                {
                    // std::cout << "In Simulation::Run --> Controlling winches... " << std::endl;
                    wallTimeControllerWinches += winchesContTimeStep;
                    WinchesController->controlWinchies(wallTimeControllerWinches);
                    WinchesController->WriteOut(wallTimeControllerWinches);
                    // std::cout << "In Simulation::Run --> ... done controlling winches!" << std::endl;
                }
            }
        }
        pTimeSolver->finalize();
    }
    else
    {
        // Simple time-stepping for systems with zero DOFs (e.g., fixed body under waves)
        Logger::info("Running static analysis (no DOFs)...");

        int num_steps = (int)std::ceil(simulationTime / writeTimeStep);
        int output_step = 0;

        for (int step = 0; step <= num_steps; step++)
        {
            double current_time = step * writeTimeStep;

            // Don't go beyond simulation time
            if (current_time > simulationTime + 1e-10)
                break;

            // Calculate forces and dynamics at current time
            CalculateSystemDynamicsStatic(current_time);

            // Write output at specified intervals
            if (current_time >= wallTime - 1e-12)
            {
                if (output_step > 0)
                {
                    std::ostringstream _p; _p << std::fixed << std::setprecision(2) << "  t = " << current_time << " / " << simulationTime << " s"; Logger::progress(_p.str());
                }
                output_step++;

                for (int ii = 0; ii < numLines; ii = ii + 1)
                    pLines[ii]->WriteOut(current_time);
                for (int ii = 0; ii < numBodies; ii = ii + 1)
                    pBodies[ii]->WriteOut(current_time);
                for (int ii = 0; ii < numOWCs; ii = ii + 1)
                    pOWCs[ii]->WriteOut(current_time);
                for (int ii = 0; ii < numSprings; ii = ii + 1)
                    pSprings[ii]->WriteOut(current_time);

                wallTime += writeTimeStep; // Increment to next output time
            }

            // Update hydrodynamic forces at specified intervals
            if (current_time >= wallTimeHydro + hydroTimeStep - 1e-12)
            {
                wallTimeHydro += hydroTimeStep;
                if (numBodies > 0)
                {
                    for (int ii = 0; ii < numBodies; ii = ii + 1)
                    {
                        pBodies[ii]->Fb = pBodies[ii]->pHydro->CalculateHydrodynamicForces(wallTimeHydro);
                    }
                }
            }

            // Update wind turbine forces at specified intervals
            if (numWindTurbines > 0)
            {
                if (current_time >= wallTimeFAST + fastTimeStep - 1e-12)
                {
                    wallTimeFAST += fastTimeStep;
                    for (int ii = 0; ii < numWindTurbines; ii = ii + 1)
                    {
                        pWindTurbines[ii]->SetInputsFAST();
                        pWindTurbines[ii]->ComputeForces(wallTimeFAST);
                        pWindTurbines[ii]->WriteOut(wallTimeFAST);
                    }
                }
                if (current_time >= wallTimeControllerFAST + fastControllerTimeStep - 1e-12)
                {
                    wallTimeControllerFAST += fastControllerTimeStep;
                    for (int ii = 0; ii < numWindTurbines; ii = ii + 1)
                    {
                        pWindTurbines[ii]->SetInputsFAST();
                        pWindTurbines[ii]->ComputeControler(wallTimeControllerFAST);
                    }
                }
            }

            // Update OWC controllers at specified intervals
            if (numOWCs > 0)
            {
                if (current_time >= wallTimeControllerOWCs + owcsContTimeStep - 1e-12)
                {
                    wallTimeControllerOWCs += owcsContTimeStep;
                    for (int ii = 0; ii < numOWCs; ii = ii + 1)
                    {
                        if (pOWCs[ii]->turbine_type > 0)
                        {
                            pOWCs[ii]->pOWCTurbine->ComputeGenTorque();
                        }
                    }
                }
            }

            // Update sinking dynamics at specified intervals
            if (numSinking > 0)
            {
                if (current_time >= wallTimeSinking + sinkingTimeStep - 1e-12)
                {
                    wallTimeSinking += sinkingTimeStep;
                    for (int ii = 0; ii < numSinking; ii = ii + 1)
                    {
                        pSinking[ii]->UpdateSinkingHydrodynamics(current_time);
                    }
                    UpdateSystemMatrix();
                    for (int ii = 0; ii < numSinking; ii = ii + 1)
                        pSinking[ii]->WriteOut(current_time);
                }
            }

            // Update winch controllers at specified intervals
            if (numWinches > 0)
            {
                if (current_time >= wallTimeControllerWinches + winchesContTimeStep - 1e-12)
                {
                    wallTimeControllerWinches += winchesContTimeStep;
                    WinchesController->controlWinchies(wallTimeControllerWinches);
                    WinchesController->WriteOut(wallTimeControllerWinches);
                }
            }

            // Note: time is incremented via the for loop counter in step
        }
    }

    tend = time(0);
    double computational_time = difftime(tend, tstart);
    {
        std::ostringstream _s;
        if (computational_time < 60)
            _s << std::fixed << std::setprecision(1) << "  Computational time : " << computational_time << " s";
        else if (computational_time < 3600)
            _s << std::fixed << std::setprecision(2) << "  Computational time : " << computational_time / 60.0 << " min";
        else
            _s << std::fixed << std::setprecision(2) << "  Computational time : " << computational_time / 3600.0 << " h";
        Logger::info(_s.str());
    }
    if (numSystem > 0)
    {
        Logger::info("  Total function calls  : " + std::to_string(numCallsSysFun));
        Logger::info("  Total Jacobian updates: " + std::to_string(pTimeSolver->iJ));
        Logger::info("  Avg Newton iterations : " + std::to_string(pTimeSolver->nNewtonIterAvg));
        if (pTimeSolver->nConvergenceFailed > 0)
            Logger::warning("  Convergence failed " + std::to_string(pTimeSolver->nConvergenceFailed) + " time(s).");
    }
    if (writeEquilibrium == 1)
    {
        Logger::info("  Writing initial condition to dataStaticIC.dat...");
        std::string filename = JoinPath(outputFolderPath, "dataStaticIC.dat");
        pTimeSolver->y.save(filename, arma::arma_ascii);
    }
}

void Simulation::SetupCase()
{
    Logger::info("--> Setting up case configuration...");

    // Count the number of BCP in each body and create pointer array
    if (numBodies > 0)
    {
        Logger::debug("Counting BCPs per body");
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        for (int jj = 0; jj < pBodies[ii]->numBcps; jj++)
        {
            if (pBodies[ii]->pIndexBcps[jj] + 1 > numBcps)
            {
                std::stringstream ss;
                ss << "BCP index: " << pBodies[ii]->pIndexBcps[jj] << " in Body: " << pBodies[ii]->GetId()
                   << " is out of range when compare with the Number of BCPs(" << numBcps << ") defined in"
                   << " dataBCPs.dat";
                throw ValueError(ss.str());
            }
            pBcps[pBodies[ii]->pIndexBcps[jj]]->numBodiesBcp++;
        }
    }
    bool* pDefined_body_bcps = new bool[numBcps];
    for (int ii = 0; ii < numBcps; ii++)
    {
        pDefined_body_bcps[ii] = 0;
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        for (int jj = 0; jj < pBodies[ii]->numBcps; jj++)
        {
            if (!pDefined_body_bcps[pBodies[ii]->pIndexBcps[jj]])
            {
                pBcps[pBodies[ii]->pIndexBcps[jj]]->pBodies =
                    new Body*[pBcps[pBodies[ii]->pIndexBcps[jj]]->numBodiesBcp];
                pDefined_body_bcps[pBodies[ii]->pIndexBcps[jj]] = true;
            }
        }
    }
    delete[] pDefined_body_bcps;
    if (numBodies > 0)
    {
        Logger::debug("BCPs per body done");
    }

    // Floor triangulation initialization
    if (numFloor > 0)
    {
        Logger::debug("Initializing floor");
    }
    for (int ii = 0; ii < numFloor; ii++)
    {
        if (pSeaFloor[ii]->GetType() == 3)
        {
            dynamic_cast<Bathymetry*>(pSeaFloor[ii])->getVertexNormals();
            dynamic_cast<Bathymetry*>(pSeaFloor[ii])->getProjectionMatrix();
        }
        else if (pSeaFloor[ii]->GetType() == 2)
        {
            dynamic_cast<Inclined*>(pSeaFloor[ii])->getPlaneEquation();
        }
        else if (pSeaFloor[ii]->GetType() == 1)
        {
            Logger::debug("Floor level: " + std::to_string(pSeaFloor[ii]->seabedDepth) + " m");
        }
    }
    if (numFloor > 0)
    {
        Logger::debug("Floor initialization done");
    }

    // Assing to each BCP the corresponding Body pointer
    if (numBodies > 0)
    {
        Logger::debug("Assigning BCPs to bodies");
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        {
            std::ostringstream oss;
            oss << "Body " << ii + 1 << " pos = " << pBodies[ii]->pos.t();
            Logger::debug(oss.str());
        }
        for (int jj = 0; jj < pBodies[ii]->numBcps; jj++)
        {
            pBodies[ii]->pBodyBcps[jj] = pBcps[pBodies[ii]->pIndexBcps[jj]];
            pBodies[ii]->pBodyBcps[jj]->pBodies[pBodies[ii]->pBodyBcps[jj]->countBody] = pBodies[ii];
            pBodies[ii]->pBodyBcps[jj]->countBody++;
        }
        pBodies[ii]->UpdateBcps();
    }
    if (numBodies > 0)
    {
        Logger::debug("BCPs assigned to bodies");
    }

    // Count the number of Wind Turbines in each body and create pointer array
    if (numBodies > 0)
    {
        Logger::debug("Checking wind turbines per body");
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        for (int jj = 0; jj < pBodies[ii]->numWindTurbs; jj++)
        {
            if (pBodies[ii]->pIndexWindTurbs[jj] + 1 > numWindTurbines)
            {
                std::stringstream ss;
                ss << "Wind Turbine index: " << pBodies[ii]->pIndexWindTurbs[jj] << " in Body: " << pBodies[ii]->GetId()
                   << " is out of range when compare with the Number of Wind Turbines (" << numWindTurbines
                   << ") defined in"
                   << " dataWindTurbines.dat";
                throw ValueError(ss.str());
            }
        }
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        for (int jj = 0; jj < numBodies; jj++)
        {
            if (ii != jj && pBodies[ii]->numWindTurbs > 0 && pBodies[jj]->numWindTurbs > 0)
            {
                for (int kii = 0; kii < pBodies[ii]->numWindTurbs; kii++)
                {
                    int indWTloc = pBodies[ii]->pIndexWindTurbs[kii];
                    for (int kjj = 0; kjj < numBodies; kjj++)
                    {
                        if (indWTloc == pBodies[jj]->pIndexWindTurbs[kjj])
                        {
                            std::stringstream ss;
                            ss << "Bodies " << ii + 1 << " and " << jj + 1 << " share wind turbine" << indWTloc << "!";
                            throw ValueError(ss.str());
                        }
                    }
                }
            }
        }
    }
    if (numBodies > 0)
    {
        Logger::debug("Wind turbine check done");
    }

    // Assing to each Body the corresponding Wind Turbine pointers
    if (numBodies > 0)
    {
        Logger::debug("Assigning wind turbines to bodies");
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        for (int jj = 0; jj < pBodies[ii]->numWindTurbs; jj++)
        {
            pBodies[ii]->pBodyWindTurbs[jj] = pWindTurbines[pBodies[ii]->pIndexWindTurbs[jj]];
        }
    }
    if (numBodies > 0)
    {
        Logger::debug("Wind turbines assigned to bodies");
    }

    // Count the number of lines in each joint BCP and elastic anchor BCP
    if (numBcps > 0)
    {
        Logger::debug("Counting lines per joint/elastic BCP");
    }
    for (int ii = 0; ii < numBcps; ii++)
    {
        if (pBcps[ii]->GetType() == 3 || pBcps[ii]->GetType() == 5)
        {
            int temp_nL = 0;
            int temp_BCP_Id = pBcps[ii]->GetId();
            for (int jj = 0; jj < numLines; jj++)
            {
                if (pLines[jj]->indexBcps[0] == temp_BCP_Id)
                {
                    temp_nL++;
                }
                if (pLines[jj]->indexBcps[1] == temp_BCP_Id)
                {
                    temp_nL++;
                }
            }
            pBcps[ii]->posLines = arma::zeros(temp_nL, 3);
            pBcps[ii]->velLines = arma::zeros(temp_nL, 3);
        }
    }
    if (numBcps > 0)
    {
        Logger::debug("Joint/elastic BCP line count done");
    }

    // Count the number of Lines in each body and create pointer array
    if (numLines > 0)
    {
        Logger::debug("Counting lines per BCP");
    }
    for (int ii = 0; ii < numLines; ii++)
    {
        for (int jj = 0; jj < pLines[ii]->numBcps; jj++)
        {
            if (pLines[ii]->indexBcps[jj] + 1 > numBcps)
            {
                std::stringstream ss;
                ss << "BCP index: " << pLines[ii]->indexBcps[jj] << " in Line: " << pLines[ii]->GetId()
                   << " is out of range when compare with the Number of BCPs(" << numBcps << ") defined in"
                   << " dataBCPs.dat";
                throw ValueError(ss.str());
            }
            pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp++;
        }
    }
    if (numLines > 0)
    {
        Logger::debug("Lines per BCP count done");
    }
    if (numLines > 0)
    {
        Logger::debug("Initializing BCP line data");
    }
    bool* pDefined_lines_bcps = new bool[numBcps];
    for (int ii = 0; ii < numBcps; ii++)
    {
        pDefined_lines_bcps[ii] = false;
    }
    for (int ii = 0; ii < numLines; ii++)
    {
        for (int jj = 0; jj < pLines[ii]->numBcps; jj++)
        {
            if (!pDefined_lines_bcps[pLines[ii]->indexBcps[jj]])
            {
                pBcps[pLines[ii]->indexBcps[jj]]->pLines = new Line*[pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp];
                pBcps[pLines[ii]->indexBcps[jj]]->pBcpLineIndex =
                    new int[pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp];
                pBcps[pLines[ii]->indexBcps[jj]]->pBcpLineNode = new int[pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp];
                pDefined_lines_bcps[pLines[ii]->indexBcps[jj]] = true;
            }
        }
    }
    delete[] pDefined_lines_bcps;
    if (numLines > 0)
    {
        Logger::debug("BCP line data init done");
    }

    // Assing to each Line the corresponding BCP pointer
    if (numLines > 0)
    {
        Logger::debug("Assigning BCPs to lines");
    }
    for (int ii = 0; ii < numLines; ii++)
    {
        for (int jj = 0; jj < pLines[ii]->numBcps; jj++)
        {
            pLines[ii]->pLineBcps[jj] = pBcps[pLines[ii]->indexBcps[jj]];
            pLines[ii]->pLineBcps[jj]->pLines[pLines[ii]->pLineBcps[jj]->countLine] = pLines[ii];
            pLines[ii]->pLineBcps[jj]->pBcpLineIndex[pLines[ii]->pLineBcps[jj]->countLine] = ii;
            pLines[ii]->pLineBcps[jj]->pBcpLineNode[pLines[ii]->pLineBcps[jj]->countLine] = jj;
            pLines[ii]->pLineBcps[jj]->countLine++;
        }

        int type0 = pLines[ii]->pLineBcps[0]->GetType();
        int type1 = pLines[ii]->pLineBcps[1]->GetType();
        bool isType0Dynamic = (type0 == 3 || type0 == 5); // Joint or ElasticAnchor
        bool isType1Dynamic = (type1 == 3 || type1 == 5); // Joint or ElasticAnchor

        if (!isType0Dynamic && !isType1Dynamic)
        {
            numDofLinesTotal += (pLines[ii]->N - 2);
            pLines[ii]->first_node = 1;
            pLines[ii]->last_node = pLines[ii]->N - 1;
        }
        else if (isType0Dynamic && isType1Dynamic)
        {
            numDofLinesTotal += pLines[ii]->N;
            pLines[ii]->first_node = 0;
            pLines[ii]->last_node = pLines[ii]->N;
        }
        else if (isType0Dynamic && !isType1Dynamic)
        {
            numDofLinesTotal += (pLines[ii]->N - 1);
            pLines[ii]->first_node = 0;
            pLines[ii]->last_node = pLines[ii]->N - 1;
        }
        else if (!isType0Dynamic && isType1Dynamic)
        {
            numDofLinesTotal += (pLines[ii]->N - 1);
            pLines[ii]->first_node = 1;
            pLines[ii]->last_node = pLines[ii]->N;
        }

        // TODO: Remove this try-catch block and manage the exceptions inside the Line class
        try
        {
            // TODO: Check how this interacts with the FEM initial position
            if (!readEquilibrium)
            {
                pLines[ii]->initLine();
            }
            pLines[ii]->SEM_getBaseFunctions();

            // TODO: Implement a logger with different levels of verbosity
            // pLines[ii]->print_out();
        }
        catch (int e)
        {
            if (e == 0)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is under the floor level.");
            if (e == 1)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " touches the seafloor and it shouldn't.");
            if (e == 2)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is not tense and laying on the seafloor. It should be pretensed.");
            if (e == 3)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is not tense and vertical. It should be pretensed.");
            if (e == 4)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " is not tense and it should.");
            if (e == 5)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " initial shape can't be computed with QS method.");
            if (e == 6)
                Logger::error("Line " + std::to_string(pLines[ii]->nLine) + " touches the seafloor although none of its ends are there.");
        }

        pLines[ii]->pLineSeaFloor = pSeaFloor[pLines[ii]->indexSeaFloor];
    }
    if (numLines > 0)
    {
        Logger::debug("BCPs assigned to lines");
    }

    // Count the number of line nodes without repetition of joint nodes
    if (numLines > 0)
    {
        Logger::debug("Counting unique line nodes");
    }
    int numUsedJointBCPs = 0;
    for (int jj = 0; jj < numLines; jj++)
    {
        numAllLinesNodes += pLines[jj]->N;
        for (int kk = 0; kk < 2; kk++)
        {
            if (pLines[jj]->pLineBcps[kk]->GetType() == 3 || pLines[jj]->pLineBcps[kk]->GetType() == 5)
            {
                if (pLines[jj]->pLineBcps[kk]->flag_counted == 0)
                {
                    pLines[jj]->pLineBcps[kk]->flag_counted = 1;
                }
                else
                {
                    numUsedJointBCPs++;
                }
            }
        }
    }
    numAllLinesNodes -= numUsedJointBCPs;
    if (numLines > 0)
    {
        Logger::debug("Unique line node count done");
    }

    // Build lines coupling sparse matrix
    if (numLines > 0)
    {
        Logger::debug("Building lines coupling sparse matrix");
    }
    int indFirstNodeAvail = 0;
    for (int jj = 0; jj < numLines; jj++)
    {
        // Store the number of nodes in the line
        int numLineNodes_tmp = pLines[jj]->N;

        // Allocate memory for the index vector
        pLines[jj]->ind4CouplingMat = arma::zeros<arma::uvec>(numLineNodes_tmp);

        // First node
        if (pLines[jj]->pLineBcps[0]->flag_assigned == 0)
        {
            pLines[jj]->pLineBcps[0]->flag_assigned = 1;
            pLines[jj]->pLineBcps[0]->couplingMatIndex = indFirstNodeAvail;
            indFirstNodeAvail++;
        }
        pLines[jj]->ind4CouplingMat(0) = pLines[jj]->pLineBcps[0]->couplingMatIndex;

        // Intermediate nodes
        for (int kk = 1; kk < numLineNodes_tmp - 1; kk++)
        {
            pLines[jj]->ind4CouplingMat(kk) = indFirstNodeAvail;
            indFirstNodeAvail++;
        }

        // Last node
        if (pLines[jj]->pLineBcps[1]->flag_assigned == 0)
        {
            pLines[jj]->pLineBcps[1]->flag_assigned = 1;
            pLines[jj]->pLineBcps[1]->couplingMatIndex = indFirstNodeAvail;
            indFirstNodeAvail++;
        }
        pLines[jj]->ind4CouplingMat(numLineNodes_tmp - 1) = pLines[jj]->pLineBcps[1]->couplingMatIndex;
    }
    if (numLines > 0)
    {
        Logger::debug("Coupling matrix index assignment done");
    }

    // Compute the lines coupling matrices
    if (numLines > 0)
    {
        Logger::debug("Computing lines coupling matrices");
    }
    pLinesCouplingMatrix = new arma::mat(numAllLinesNodes, numAllLinesNodes, arma::fill::zeros);
    pLinesCouplingMatrixInv = new arma::mat(numAllLinesNodes, numAllLinesNodes);
    pLinesCouplingMatrix_sp = new arma::sp_mat(numAllLinesNodes, numAllLinesNodes);
    ComputeLinesCouplingMatrix();
    // TODO: 100 should be a parameter
    if (numAllLinesNodes < 100 && !useWinches)
    {
        *pLinesCouplingMatrixInv = arma::inv(*pLinesCouplingMatrix);
    }
    if (numLines > 0)
    {
        Logger::debug("Lines coupling matrices done");
    }

    // Compute equilibrium with FEM for all lines at the same time
    if (!readEquilibrium && flagStatic == 1 && numLines > 0) //&& flagStatic == 0?
    {
        Logger::info("Computing FEM equilibrium for all lines...");

        // Store the lines friction model and tension flag
        arma::umat flagLineas;
        arma::umat flagTension;
        flagLineas.zeros(numLines, 1);
        flagTension.zeros(numLines, 1);
        for (int ii = 0; ii < numLines; ii++)
        {
            // Call GetValues only for non-dynamic BCPs (skip type 5 elastic anchors to preserve their initial pos)
            if (pLines[ii]->pLineBcps[0]->GetType() != 5)
                pLines[ii]->pLineBcps[0]->GetValues(0.0);
            if (pLines[ii]->pLineBcps[1]->GetType() != 5)
                pLines[ii]->pLineBcps[1]->GetValues(0.0);
            // Store the original values
            flagLineas(ii) = pLines[ii]->frictionModel;
            flagTension(ii) = pLines[ii]->flag_tension;
            // Set them to a value that allows the equilibrium computation
            pLines[ii]->frictionModel = 0;
            pLines[ii]->flag_tension = 1;
        }

        // Start the equilibrium computation with an initial guess
        arma::mat posicionInicial = ComputeLinesInitialPoint();
        // Perform the equilibrium computation
        ComputeLinesEquilibrium(posicionInicial);
        Logger::info("--> Initial static position computed");
        // Return the original values of the friction model and tension flag
        for (int i = 0; i < numLines; i++)
        {
            pLines[i]->frictionModel = flagLineas(i);
            pLines[i]->flag_tension = flagTension(i);
        }

        Logger::debug("FEM equilibrium done");
    }

    // Initialize strain history vector if viscoelasticity model is used
    for (int i = 0; i < numLines; i++)
    {
        if (pLines[i]->flag_stiffness == 1)
        {
            pLines[i]->initiallize_strain_memory();
        }
    }
    Logger::debug("Strain memory initialized");
    // Store the initial position, required for the stick-slip friction model
    if (numLines > 0)
    {
        Logger::debug("Storing initial line positions");
    }
    for (int i = 0; i < numLines; i++)
    {
        pLines[i]->frictionPos = pLines[i]->pLineSeaFloor->projectPoints(pLines[i]->pos)(0, 0);
    }
    if (numLines > 0)
    {
        Logger::debug("Initial line positions stored");
    }

    // Setup Springs
    if (numSprings > 0)
    {
        Logger::debug("Setting up springs");
    }
    for (int ii = 0; ii < numSprings; ii++)
    {
        Logger::debug("Spring " + std::to_string(ii) + ": BCP_1=" + std::to_string(pSprings[ii]->BCP_1) + " BCP_2=" + std::to_string(pSprings[ii]->BCP_2));
        pSprings[ii]->SpringBCP[0] = pBcps[pSprings[ii]->BCP_1];
        pSprings[ii]->SpringBCP[1] = pBcps[pSprings[ii]->BCP_2];
        pSprings[ii]->pSim = this;
        pSprings[ii]->OpenOutputFiles(outputFolderPath);
    }
    if (numSprings > 0)
    {
        Logger::debug("Springs setup done");
    }

    // Setup hidro data bases
    if (numBodies > 0)
    {
        Logger::debug("Setting up hydro databases");
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        pBodies[ii]->pHydro->SetUp();
    }
    if (numBodies > 0)
    {
        Logger::debug("Hydro databases setup done");
    }

    // Setup winchies controller
    if (useWinches)
    {
        Logger::debug("Setting up winch controller");
        WinchesController->SetUpWinchiesController();
        Logger::debug("Winch controller setup done");
    }

    if (numWindTurbines > 0)
    {
        // Setup wind turbines
        Logger::debug("Setting up wind turbines");
        for (int ii = 0; ii < numWindTurbines; ii++)
        {
            pWindTurbines[ii]->Initialize();
        }
        Logger::debug("Wind turbines setup done");

        // Computing bodies structural mass considering wind turbines
        Logger::debug("Computing body inertia with wind turbine contributions");
        for (int ii = 0; ii < numBodies; ii++)
        {
            if (pBodies[ii]->numWindTurbs > 0)
            {
                // TODO: Here, it would be convenient to check that for the different turbines
                // the inertia of the HDB are at least very similar.
                pBodies[ii]->inertia = pBodies[ii]->pBodyWindTurbs[0]->bodyInerMat;
                for (int jj = 0; jj < pBodies[ii]->numWindTurbs; jj++)
                {

                    pBodies[ii]->inertia += pBodies[ii]->pBodyWindTurbs[0]->towrInerMat;
                    pBodies[ii]->inertia += pBodies[ii]->pBodyWindTurbs[0]->turbInerMat;
                }
            }
        }
        Logger::debug("Body inertia with turbines done");
    }

    Logger::info("--> Case configuration done");
}

void Simulation::ComputeLinesCouplingMatrix(void)
{
    // Compute the matrices used to couple the lines which are joined together with a Joint type BCP

    // Loop over all the lines to compute the mass matrix from the lines
    for (int jj = 0; jj < numLines; jj++)
    {
        int numLineNodes_tmp = pLines[jj]->N;

        Logger::debug("Assembling local mass matrix for line " + std::to_string(jj));
        arma::mat LineMassMat_tmp = pLines[jj]->MM * pLines[jj]->dL;
        if (pLines[jj]->pLineBcps[0]->GetType() != 3 && pLines[jj]->pLineBcps[0]->GetType() != 5)
        {
            LineMassMat_tmp.row(0) = arma::zeros(1, numLineNodes_tmp);
            LineMassMat_tmp(0, 0) = 1.0;
        }
        if (pLines[jj]->pLineBcps[1]->GetType() != 3 && pLines[jj]->pLineBcps[1]->GetType() != 5)
        {
            LineMassMat_tmp.row(numLineNodes_tmp - 1) = arma::zeros(1, numLineNodes_tmp);
            LineMassMat_tmp(numLineNodes_tmp - 1, numLineNodes_tmp - 1) = 1.0;
        }

        // Add the mass matrix to the global matrix in sparse or full format depending on the size
        // TODO: 100 should be a parameter
        Logger::debug("Adding local mass matrix to global matrix");
        if (numAllLinesNodes >= 100)
        {
            for (int irow = 0; irow < numLineNodes_tmp; irow++)
            {
                for (int icol = 0; icol < numLineNodes_tmp; icol++)
                {
                    (*pLinesCouplingMatrix_sp)(pLines[jj]->ind4CouplingMat(irow), pLines[jj]->ind4CouplingMat(icol)) +=
                        LineMassMat_tmp(irow, icol);
                }
            }
        }
        else
        {
            (*pLinesCouplingMatrix).submat(pLines[jj]->ind4CouplingMat, pLines[jj]->ind4CouplingMat) += LineMassMat_tmp;
        }
    }

    // Loop over all the BCPs to compute the mass matrix from the joints and elastic anchors
    for (int jj = 0; jj < numBcps; jj++)
    {
        if ((pBcps[jj]->GetType() == 3 || pBcps[jj]->GetType() == 5) && (pBcps[jj]->flag_assigned == 1))
        {
            double node_mass = (pBcps[jj]->GetType() == 3) ? pBcps[jj]->mass_Joint
                                                           : dynamic_cast<ElasticAnchorBCP*>(pBcps[jj])->anchor_mass;
            // TODO: 100 should be a parameter
            if (numAllLinesNodes >= 100)
            {
                (*pLinesCouplingMatrix_sp)(pBcps[jj]->couplingMatIndex, pBcps[jj]->couplingMatIndex) += node_mass;
            }
            else
            {
                (*pLinesCouplingMatrix)(pBcps[jj]->couplingMatIndex, pBcps[jj]->couplingMatIndex) += node_mass;
            }
        }
    }
}

Simulation::Simulation(std::string incProjectPath, std::string incDataFormat)
{
    // Assign direct variables
    projectFolderPath = incProjectPath;

    // Process indirect variables
    if (!incDataFormat.compare("ASCII"))
    {
        dataFormat = 0;
        dataFormatStr = "ASCII";
        inputFolderPath = JoinPath(incProjectPath, "input");
        outputFolderPath = JoinPath(incProjectPath, "output");
        pReadProperties = &Simulation::ReadPropertiesASCII;
        pReadWaves = &Simulation::ReadWavesASCII;
        pReadBcps = &Simulation::ReadBcpsASCII;
        pReadBodies = &Simulation::ReadBodiesASCII;
        pReadSinking = &Simulation::ReadSinkingASCII;
        pReadLines = &Simulation::ReadLinesASCII;
        pReadSprings = &Simulation::ReadSpringsASCII;
        pReadWinches = &Simulation::ReadWinchesASCII;
        pReadSeaFloor = &Simulation::ReadSeaFloorASCII;
        pReadWindTurbines = &Simulation::ReadWindTurbinesASCII;
        pReadOWCs = &Simulation::ReadOWCsASCII;
    }
    else if (!incDataFormat.compare("YAML"))
    {
        dataFormat = 1;
        dataFormatStr = "YAML";
        inputFolderPath = JoinPath(incProjectPath, "input");
        outputFolderPath = JoinPath(incProjectPath, "output");
        pReadProperties = &Simulation::ReadPropertiesYAML;
        pReadWaves = &Simulation::ReadWavesYAML;
        pReadBcps = &Simulation::ReadBcpsYAML;
        pReadBodies = &Simulation::ReadBodiesYAML;
        pReadSinking = &Simulation::ReadSinkingYAML;
        pReadLines = &Simulation::ReadLinesYAML;
        pReadSprings = &Simulation::ReadSpringsYAML;
        pReadWinches = &Simulation::ReadWinchesYAML;
        pReadSeaFloor = &Simulation::ReadSeaFloorYAML;
        pReadWindTurbines = &Simulation::ReadWindTurbinesYAML;
        pReadOWCs = &Simulation::ReadOWCsYAML;
    }
    else
    {
        std::stringstream ss;
        ss << "Simulation data format --> " << incDataFormat
           << " is not available.\n    Available formats: ASCII | YAML.";
        throw ValueError(ss.str());
    }
}

void Simulation::UpdateSystem()
{
    // Call overloaded version with pTimeSolver->t as the time parameter
    if (pTimeSolver != nullptr)
    {
        UpdateSystem(pTimeSolver->t);
    }
}

void Simulation::UpdateSystem(double time)
{
    // Update the bodies velocity buffers for the radiation forces computation
    if (numBodies > 0)
    {
        // Update time vector if any
        bool restoreMatrix = false;
        timeBufferCount++;
        if (timeBufferCount < timeBufferSize)
        {
            timeBuffer(0, timeBufferCount) = time;
        }
        else
        {
            arma::mat timeBufferNew = arma::zeros(1, timeBufferSize);
            timeBufferNew.cols(0, pBodies[0]->pHydro->GetNumPointsIrf() - 1) =
                timeBuffer.cols(timeBufferSize - pBodies[0]->pHydro->GetNumPointsIrf(), timeBufferSize - 1);
            timeBuffer = timeBufferNew;
            timeBufferCount = pBodies[0]->pHydro->GetNumPointsIrf() - 1;
            restoreMatrix = true;
        }

        // Update bodies velocity
        int ini = 0;
        for (int ii = 0; ii < numBodies; ii++)
        {
            pBodies[ii]->StoreVelocities(restoreMatrix);
            ini = ini + 6;
        }
    }
}

void Simulation::UpdateSystemMatrix()
{
    // Update the system matrix and its inverse (used when some body is sinking)
    *pSystemMatrix = arma::zeros(6 * numBodies, 6 * numBodies);
    for (int ii = 0; ii < this->numBodies; ii++)
    {
        (*pSystemMatrix)(pBodies[ii]->sysMatSpan1, pBodies[ii]->sysMatSpan2) += pBodies[ii]->pHydro->GetTotalMass();
    }
    *pSystemMatrixInv = arma::solve(*pSystemMatrix, eye(size(*pSystemMatrix)));
}

arma::mat Simulation::ComputeLinesInitialPoint()
{
    // Matrix to store the positions of the nodes from all lines
    arma::mat positionLinesCouplingVector = arma::zeros(numAllLinesNodes, 3);
    // Vector to store the indexes of the lines that are not joint, start assuming that all are joints
    arma::mat isNotJoint = arma::zeros(numAllLinesNodes, 1);

    // Loop over all the lines to compute the initial position matrix
    for (int ii = 0; ii < numLines; ii++)
    {
        // Initialize the initial position matrix with zeros for current line
        arma::mat initial_positions = arma::zeros(pLines[ii]->N, 3);

        // Get the position of the first and last nodes
        arma::mat pos1 = arma::strans(pLines[ii]->pLineBcps[0]->pos);
        arma::mat posN = arma::strans(pLines[ii]->pLineBcps[1]->pos);

        // Fill the initial position matrix with the first and last nodes positions and
        // the rest of the nodes with the previous values
        initial_positions = pLines[ii]->pos;
        initial_positions.row(0) = pos1;
        initial_positions.row(pLines[ii]->N - 1) = posN;

        // Insert current line positions in the global vector
        positionLinesCouplingVector.rows(pLines[ii]->ind4CouplingMat) = initial_positions;

        // Store the initial position in the line object
        pLines[ii]->pos = initial_positions;

        // Look for the nodes that are not joints or elastic anchors, and store the indexes
        for (int k = 0; k < 2; k++)
        {
            if (pLines[ii]->pLineBcps[k]->GetType() != 3 && pLines[ii]->pLineBcps[k]->GetType() != 5)
            {
                isNotJoint(pLines[ii]->pLineBcps[k]->couplingMatIndex) = 1;
            }
        }
    }

    // A vector without anchors or fairleads is to be returned
    indexesFairAnchor = arma::find(isNotJoint);          // Indices of the nodes that are not joints
    indexesNoFairNoAnchor = arma::find(isNotJoint == 0); // Indices of the nodes that are joints
    // Store the removed rows
    removed_rows = positionLinesCouplingVector.rows(indexesFairAnchor);
    // Remove rows that are not joints
    positionLinesCouplingVector.shed_rows(indexesFairAnchor);
    // Convert to a column vector
    arma::mat output_position =
        arma::reshape(arma::strans(positionLinesCouplingVector), 3 * positionLinesCouplingVector.n_rows, 1);

    return output_position;
}

arma::mat Simulation::ComputeLinesForces(arma::mat posicion)
{
    // Input is a column vector with the positions of the nodes that are not anchors or fairleads
    // We need to convert it to a matrix with 3 columns

    // Initiallize the matrix with zeros
    arma::mat full_matrix = arma::zeros(numAllLinesNodes, 3);
    // Include the static nodes positions
    full_matrix.rows(indexesFairAnchor) = removed_rows;
    // Reshape the input vector to a matrix
    arma::mat posicionMatriz = arma::strans(arma::reshape(posicion, 3, numAllLinesNodes - removed_rows.n_rows));
    // Include the dynamic nodes positions
    full_matrix.rows(indexesNoFairNoAnchor) = posicionMatriz;

    // Initiallize the vector to store the forces
    arma::mat forcesLinesCouplingVector = arma::zeros(numAllLinesNodes, 3);

    // Loop over all the lines to compute the forces
    for (int ii = 0; ii < numLines; ii++)
    {
        // Get the initial positions of the current line
        pLines[ii]->pos = full_matrix.rows(pLines[ii]->ind4CouplingMat);
        // Set the line nodes velocity to zero
        pLines[ii]->vel = arma::zeros(pLines[ii]->N, 3);
        // Compute the forces
        // std::cout << "Simulation::ComputeLinesForces: Compute the forces..." << std::endl;
        pLines[ii]->SEM_computeF(0.0);
        // Store the forces in the global vector
        forcesLinesCouplingVector.rows(pLines[ii]->ind4CouplingMat) += pLines[ii]->F;
    }

    // Populate posLines/velLines for dynamic BCPs (joints and elastic anchors)
    for (int ii = 0; ii < numLines; ii++)
    {
        for (int kk = 0; kk < 2; kk++)
        {
            int bcpType = pLines[ii]->pLineBcps[kk]->GetType();
            if (bcpType == 3 || bcpType == 5)
            {
                int lineNode = (kk == 0) ? 0 : pLines[ii]->N - 1;
                pLines[ii]->pLineBcps[kk]->posLines.row(pLines[ii]->pLineBcps[kk]->iLJ) = pLines[ii]->pos.row(lineNode);
                pLines[ii]->pLineBcps[kk]->velLines.row(pLines[ii]->pLineBcps[kk]->iLJ) = pLines[ii]->vel.row(lineNode);
                pLines[ii]->pLineBcps[kk]->iLJ++;
            }
        }
    }
    // Compute and add joint/elastic anchor forces
    for (int ii = 0; ii < numBcps; ii++)
    {
        if ((pBcps[ii]->GetType() == 3 || pBcps[ii]->GetType() == 5) && (pBcps[ii]->flag_assigned == 1))
        {
            pBcps[ii]->GetValues(0.0);
            forcesLinesCouplingVector.row(pBcps[ii]->couplingMatIndex) += pBcps[ii]->JointForce;
        }
    }

    // Remove the forces of the nodes that are not joints
    forcesLinesCouplingVector.shed_rows(indexesFairAnchor);
    // Convert the matrix to a column vector
    arma::mat output_force =
        arma::reshape(arma::strans(forcesLinesCouplingVector), 3 * forcesLinesCouplingVector.n_rows, 1);
    // Return output force
    return output_force;

    // TODO: Check why using lines accelerations does not work properly to compute initial position
    // // Compute lines accelerations
    // arma::mat LinesAccelerations = arma::solve(*pLinesCouplingMatrix, forcesLinesCouplingVector);
    // // Remove the forces of the nodes that are not joints
    // LinesAccelerations.shed_rows(indexesFairAnchor);
    // // Convert the matrix to a column vector
    // arma::mat output_force = arma::reshape(arma::strans(LinesAccelerations), 3 * LinesAccelerations.n_rows, 1);
    // // Return output force
    // return output_force;
}

arma::mat Simulation::ComputeLinesJacobian(arma::mat position, arma::mat force)
{
    // Get the number of variables from the input
    int numVariables = position.n_rows;
    // Set the identity matrix to get the canonic base vectors
    arma::mat canonic_base = arma::eye(numVariables, numVariables);
    // Initialling the jacobian matrix
    arma::mat jacobian = arma::zeros(numVariables, numVariables);
    // Set the finite differences step parameter
    double h = 1e-12;
    // Loop over all the variables to compute the jacobian
    for (int j = 0; j < numVariables; j++)
    {
        jacobian.col(j) = (ComputeLinesForces(position + h * canonic_base.col(j)) - force) / h;
    }
    return jacobian;
}

void Simulation::ComputeLinesEquilibrium(arma::mat x)
{
    // Newton-Raphson method for lines equilibrium
    // TODO: Use a method implemented in a library instead or move this to MathTools
    // TODO: Fix bug: this does not work for unsteady initial conditions (falling lines)

    // Parameters
    int maxIter = 1000;
    double tol = 1.0e-9;
    double relTol = 1.0e-6;
    double minStep = 1e-14;

    // Initial stop criterion variables
    double cantidadRel = 2 * relTol;
    double cantidadAbs = 2 * tol;
    double step = 2 * minStep;
    int iter = 0;

    // Initialization of solution vectos
    arma::mat xsol = arma::zeros(x.n_rows, 1);
    arma::mat fxsol = arma::zeros(x.n_rows, 1);

    // Declare local variables
    arma::mat L;
    arma::mat U;
    arma::mat P;
    arma::mat jacobiano;
    arma::mat y;
    arma::mat dk;

    // Store initial force for the stop criterion
    // std::cout << "Simulation::ComputeLinesEquilibrium: Store initial force for the stop criterion." << std::endl;
    arma::mat fx = ComputeLinesForces(x);
    double normafxInicial = arma::norm(fx, 2);

    // Newton-Raphson loop
    // std::cout << "Simulation::ComputeLinesEquilibrium: Newton-Raphson loop." << std::endl;
    while ((iter < maxIter) && (step > minStep) && ((cantidadAbs > tol) || (cantidadRel > relTol)))
    {

        // Step 1: compute the Jacobian Matrix
        jacobiano = ComputeLinesJacobian(x, fx);

        // Step 2: Solve the linear system
        dk = arma::solve(jacobiano, -fx);

        // Step 3: Find the step size
        double rho = 1; // Initial step
        double sigma = 1e-4;
        double beta = 0.5;
        arma::mat new_step_vec = ComputeLinesForces(x + rho * dk);
        while (arma::norm(new_step_vec, "inf") > ((1.0 - sigma) * arma::norm(fx, "inf")) && rho > 0.01)
        {
            rho = beta * rho;
            new_step_vec = ComputeLinesForces(x + rho * dk);
        }

        // Step 3: Compute the new solution
        xsol = x + rho * dk;
        fxsol = ComputeLinesForces(xsol); // Here the position of the Line objects is updated

        // Step 4: Update the stop criterion variables
        // cantidadAbs = arma::norm(dk, 2);
        // cantidadRel = cantidadAbs / arma::norm(x, 2);
        step = arma::norm(rho * dk, "inf");
        cantidadAbs = arma::norm(fxsol, "inf");
        cantidadRel = arma::norm(dk / arma::as_scalar(arma::sqrt(arma::mean(arma::pow(x, 2)))), "inf");

        if ((iter % 100) == 0)
        {
            std::ostringstream oss;
            oss << "Iteration " << iter << "/" << maxIter << ". Max. force: " << cantidadAbs << " N (goal " << tol
                << ").  Max. rel. pos. step: " << cantidadRel << " (goal " << relTol << ").";
            Logger::debug(oss.str());
        }
        iter = iter + 1;

        // Update the solution for next iteration
        x = xsol;
        fx = fxsol;
        // TODO: Check that forces are not computed twice for the same point
    }

    // Check if the maximum number of iterations was exceeded
    if ((iter >= maxIter) || (step <= minStep))
    {
        Logger::warning("Convergence for initial condition was not perfectly achieved!");
    }
    {
        std::ostringstream oss;
        oss << "Final Iteration: " << iter << "/" << maxIter << ". Max. force: " << cantidadAbs << " N (goal " << tol
            << ").  Max. rel. pos. step: " << cantidadRel << " (goal " << relTol << ").";
        Logger::debug(oss.str());
    }
}