
#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <ctime>

#include "Simulation.hpp"
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
#if __cplusplus >= 201703L && __has_include(<filesystem>)
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
        if (this->pOWCs[ii]->turbine_type == 0)
        {
            pOWCs[ii]->rel_pressure = arma::as_scalar(y.row(ini));
            ini = ini + 1;
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
        if (pLines[ii]->pLineBcps[0]->GetType() != 3)
        {
            pLines[ii]->pLineBcps[0]->GetValues(time);
            pLines[ii]->pos.row(0) = pLines[ii]->pLineBcps[0]->pos.t();
            pLines[ii]->vel.row(0) = pLines[ii]->pLineBcps[0]->vel.t();
        }
        if (pLines[ii]->pLineBcps[1]->GetType() != 3)
        {
            pLines[ii]->pLineBcps[1]->GetValues(time);
            pLines[ii]->pos.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->pos.t();
            pLines[ii]->vel.row(pLines[ii]->N - 1) = pLines[ii]->pLineBcps[1]->vel.t();
        }
    }
    // Set boundary conditions on pos and vel of Lines if the BCP is a joint
    // std::cout << "Simulation::CalculateSystemDynamics - Set Boundary conditios if BCP is a Joint" << std::endl;
    for (int ii = 0; ii < numLines; ii++)
    {
        if (pLines[ii]->pLineBcps[0]->GetType() == 3)
        {
            pLines[ii]->pLineBcps[0]->posLines.row(pLines[ii]->pLineBcps[0]->iLJ) = pLines[ii]->pos.row(0);
            pLines[ii]->pLineBcps[0]->velLines.row(pLines[ii]->pLineBcps[0]->iLJ) = pLines[ii]->vel.row(0);
            pLines[ii]->pLineBcps[0]->iLJ = pLines[ii]->pLineBcps[0]->iLJ + 1;
        }
        if (pLines[ii]->pLineBcps[1]->GetType() == 3)
        {
            pLines[ii]->pLineBcps[1]->posLines.row(pLines[ii]->pLineBcps[1]->iLJ) = pLines[ii]->pos.row(pLines[ii]->N - 1);
            pLines[ii]->pLineBcps[1]->velLines.row(pLines[ii]->pLineBcps[1]->iLJ) = pLines[ii]->vel.row(pLines[ii]->N - 1);
            pLines[ii]->pLineBcps[1]->iLJ = pLines[ii]->pLineBcps[1]->iLJ + 1;
        }
    }
    for (int ii = 0; ii < numLines; ii++)
    {
        if (pLines[ii]->pLineBcps[0]->GetType() == 3)
        {
            pLines[ii]->pLineBcps[0]->GetValues(time);
            pLines[ii]->pos.row(0) = pLines[ii]->pLineBcps[0]->pos.t();
            pLines[ii]->vel.row(0) = pLines[ii]->pLineBcps[0]->vel.t();
        }
        if (pLines[ii]->pLineBcps[1]->GetType() == 3)
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
        pLines[ii]->SEM_computeF();
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

    // Compute hydrostatic and hydrodynamic forces
    // std::cout << "Simulation::CalculateSystemDynamics - Compute hydrodynamic and hydrostatic forces" << std::endl;
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) = pBodiesFree[ii]->Fb + pBodiesFree[ii]->pHydro->CalculateHydrostaticForces(time);
    }
    if (Fb.has_nan() | Fb.has_inf())
    {
        std::cout << std::endl
                  << "ERROR: NaN or Inf detected in Fb for hydrostatic or hydrodynamic forces" << std::endl;
        throw std::exception();
    }

    // Compute forces on BCPs
    // std::cout << "Simulation::CalculateSystemDynamics - Compute forces on BCPs" << std::endl;
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        pBodiesFree[ii]->ComputeBcpForces();
        Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) = Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) + pBodiesFree[ii]->bcpForces;
    }
    if (Fb.has_nan() | Fb.has_inf())
    {
        std::cout << std::endl
                  << "ERROR: NaN or Inf detected in Fb for BCP forces" << std::endl;
        throw std::exception();
    }

    // Add wind turbine forces
    // std::cout << "Simulation::CalculateSystemDynamics - Add wind turbine forces" << std::endl;
    for (int ii = 0; ii < numBodiesFree; ii++)
    {
        pBodiesFree[ii]->ComputeWindTurbForces();
        Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) = Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) + pBodiesFree[ii]->windTurbForces;
    }
    if (Fb.has_nan() | Fb.has_inf())
    {
        std::cout << std::endl
                  << "ERROR: NaN or Inf detected in Fb for wind turbine forces" << std::endl;
        throw std::exception();
    }

    // Compute OWCs dynamics
    // std::cout << "Simulation::CalculateSystemDynamics - Compute OWC dynamics" << std::endl;
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
        Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) = Fb(arma::span(6 * ii, 6 * (ii + 1) - 1), 0) + pBodiesFree[ii]->owcForces;
    }

    // Compute also everything for locked bodies so it can be displayed on the output files
    arma::mat dummy;
    for (int ii = 0; ii < numBodiesLock; ii++)
    {
        dummy = pBodiesLock[ii]->pHydro->CalculateHydrostaticForces(time);
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

                accB = (*pSystemMatrixFFInv) * (Fb.rows(sysMatIndFree) - (*pSystemMatrixFL) * accAll.rows(sysMatIndLock));
            }
            else
            { // If no simplification, update matrix for new rotation states

                arma::mat tmpMat = (*pSystemMatrix);
                arma::mat tmpVec = Fb;
                arma::mat auxM, auxV, aMat, aMat_dot, invRotMat, phi_dot;
                for (int ii = 0; ii < numBodiesFree; ii++)
                {
                    aMat = pBodiesFree[ii]->aMat;
                    aMat_dot = pBodies[ii]->aMat_dot;
                    invRotMat = pBodiesFree[ii]->invRotMat;
                    phi_dot = pBodiesFree[ii]->vel.rows(arma::span(3, 5));

                    auxM = tmpMat(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), arma::span(6 * ii + 3, 6 * (ii + 1) - 1));
                    tmpMat(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), arma::span(6 * ii + 3, 6 * (ii + 1) - 1)) = auxM * invRotMat * aMat;

                    auxV = tmpVec(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), 0);
                    tmpVec(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), 0) = auxV - auxM * invRotMat * aMat_dot * phi_dot -
                                                                          arma::cross(invRotMat * aMat * phi_dot, auxM * invRotMat * aMat * phi_dot);
                }

                arma::mat tmpMat_FF = tmpMat(sysMatIndFree, sysMatIndFree);
                arma::mat tmpMat_FL = tmpMat(sysMatIndFree, sysMatIndLock);
                arma::mat tmpVec_FF = tmpVec.rows(sysMatIndFree);

                accB = arma::solve(tmpMat, tmpVec - tmpMat_FL * accAll.rows(sysMatIndLock));
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
                    tmpMat(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), arma::span(6 * ii + 3, 6 * (ii + 1) - 1)) = auxM * invRotMat * aMat;

                    auxV = tmpVec(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), 0);
                    tmpVec(arma::span(6 * ii + 3, 6 * (ii + 1) - 1), 0) = auxV - auxM * invRotMat * aMat_dot * phi_dot -
                                                                          arma::cross(invRotMat * aMat * phi_dot, auxM * invRotMat * aMat * phi_dot);
                }

                accB = arma::solve(tmpMat, tmpVec);
            }
        }
    }
    if (accB.has_nan() | accB.has_inf())
    {
        std::cout << std::endl
                  << "ERROR: NaN or Inf detected in bodies accelerations" << std::endl;
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

    // Obtain Lines accelerations, imposing boundary conditions if the BCP is not a joint
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
        // Impose BCP accelerations on lines forces vectors
        if (pLines[ii]->pLineBcps[0]->GetType() != 3)
        {
            pLines[ii]->F.row(0) = pLines[ii]->pLineBcps[0]->acc.t();
        }
        if (pLines[ii]->pLineBcps[1]->GetType() != 3)
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
            std::cout << std::endl
                      << "ERROR: NaN or Inf detected in Force vector for Line " << ii << std::endl;
            throw std::exception();
        }
        LinesCouplingVector.rows(pLines[ii]->ind4CouplingMat) += pLines[ii]->F;
    }
    for (int ii = 0; ii < numBcps; ii++)
    {
        if ((pBcps[ii]->GetType() == 3) && (pBcps[ii]->flag_assigned == 1))
        {
            LinesCouplingVector.row(pBcps[ii]->couplingMatIndex) += pBcps[ii]->JointForce;
        }
    }
    if (LinesCouplingVector.has_nan() | LinesCouplingVector.has_inf())
    {
        std::cout << std::endl
                  << "ERROR: NaN or Inf detected in lines force vector" << std::endl;
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
        std::cout << std::endl
                  << "ERROR: NaN or Inf detected in lines accelerations" << std::endl;
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
    // std::cout << "Simulation::CalculateSystemDynamics - Recompute Lines coupling matrix for new dL values" << std::endl;
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
        if (this->pOWCs[ii]->turbine_type == 0)
        {
            yprime.row(ini) = pOWCs[ii]->rel_pressure_dot;
            ini = ini + 1;
        }
    }

    // Check if yprime has a NaN
    // std::cout << "Simulation::CalculateSystemDynamics - Check if yprime has a NaN" << std::endl;
    if (yprime.has_nan() | yprime.has_inf())
    {
        std::cout << std::endl
                  << "ERROR: NaN or Inf Detected! yprime = " << std::endl;
        std::cout << yprime << std::endl;
        throw std::exception();
    }

    // std::cout << "Simulation::CalculateSystemDynamics - End of fcn" << std::endl;
    return yprime;
}

void Simulation::CloseCase()
{
    for (int ii = 0; ii < numLines; ii++)
    {
        pLines[ii]->CloseOutputFilesASCII();
    }

    for (int ii = 0; ii < numBodies; ii++)
    {
        pBodies[ii]->CloseOutputFilesASCII();
    }

    for (int ii = 0; ii < numSinking; ii++)
    {
        pSinking[ii]->CloseOutputFilesASCII();
    }

    if (useWinches)
    {
        WinchesController.CloseOutputFilesASCII();
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
    std::cout << "Num. Bodies Free: " << this->numBodiesFree << std::endl;
    std::cout << "Num. Mooring DOFs Total: " << this->numDofLinesTotal << std::endl;
    std::cout << "Num. Winchies: " << this->numWinches << std::endl;
    std::cout << "Num. Wind Turbines: " << this->numWindTurbines << std::endl;
    std::cout << "Num. OWCs: " << this->numOWCs << std::endl;

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
        if (this->pOWCs[ii]->turbine_type == 0)
        {
            numSystem = numSystem + 1;
        }
        if (this->pOWCs[ii]->turbine_type > 0)
        {
            std::stringstream ss;
            ss << "Turbines are not implemented yet for OWCs! \n";
            throw ValueError(ss.str());
        }
    }

    std::cout << "System size: " << numSystem << std::endl;

    arma::mat y = arma::zeros(numSystem, 1);
    arma::mat yprime = arma::zeros(numSystem, 1);

    int ini = 0;
    std::string file_path;
    std::cout << "Initiallizing system vector..." << std::endl;
    if (this->readEquilibrium == 0)
    {
        for (int ii = 0; ii < this->numBodiesFree; ii = ii + 1)
        {
            std::cout << "  ... Including Body: " << ii + 1 << std::endl;
            for (int jj = 0; jj < this->pBodiesFree[ii]->numDofs; jj = jj + 1)
            {
                int itemp = this->pBodiesFree[ii]->pDofs[jj];
                y(ini) = this->pBodiesFree[ii]->pos(itemp, 0);
                ini = ini + 1;
            }
        }
        for (int ii = 0; ii < this->numLines; ii = ii + 1)
        {
            std::cout << "  ... Including Line: " << ii + 1 << std::endl;
            for (int jj = pLines[ii]->first_node; jj < pLines[ii]->last_node; jj = jj + 1)
            {
                y.rows(ini, ini + 2) = this->pLines[ii]->pos.row(jj).t();
                ini = ini + 3;
            }
        }
        for (int ii = 0; ii < this->numWinches; ii = ii + 1)
        {
            std::cout << "  ... Including Winch: " << ii + 1 << std::endl;
            ini = ini + 1;
        }
        for (int ii = 0; ii < this->numWindTurbines; ii = ii + 1)
        {
            std::cout << "  ... Including Wind Turbine: " << ii + 1 << std::endl;
            y(ini) = this->pWindTurbines[ii]->rotPos;
            y(numSystem2 + ini) = this->pWindTurbines[ii]->rotSpeed;
            ini = ini + 1;
        }
        ini = 0;
        for (int ii = 0; ii < this->numOWCs; ii = ii + 1)
        {
            std::cout << "  ... Including OWC: " << ii + 1 << std::endl;
            if (this->pOWCs[ii]->turbine_type == 0)
            {
                y(2 * numSystem2 + ini) = this->pOWCs[ii]->rel_pressure;
                ini = ini + 1;
            }
        }
    }
    else
    {
        std::cout << "Reading Equilibrium.dat..." << std::endl;
        file_path = JoinPath(inputFolderPath, "Equilibrio.dat");
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
    if (this->timeIntMethod == 1)
    {
        std::cout << "Initializing BDF2 temporal solver..." << std::endl;
        pTimeSolver = new BDF2(start_time, this->simulationTime, this->maxTimeStep, y, this);
        std::cout << "  BDF2 constructor done!" << std::endl;
        pTimeSolver->init();
        pTimeSolver->atol = this->timeIntAbsTol;
        pTimeSolver->rtol = this->timeIntRelTol;
        pTimeSolver->nIterMax = this->maxIterStep;
        std::cout << "  BDF2 initiallized!" << std::endl;
    }
    else if (this->timeIntMethod == 2)
    {
        std::cout << "Initializing BDF" << timeIntOrder << " temporal solver..." << std::endl;
        pTimeSolver = new BDFN(this->timeIntOrder, this->timeIntAdaptivity, start_time, this->simulationTime, this->maxTimeStep, y, this);
        std::cout << "  BDF" << timeIntOrder << " constructor done!" << std::endl;
        pTimeSolver->init();
        pTimeSolver->atol = this->timeIntAbsTol;
        pTimeSolver->rtol = this->timeIntRelTol;
        pTimeSolver->nIterMax = this->maxIterStep;
        std::cout << "  BDF" << timeIntOrder << " initiallized!" << std::endl;
    }

    // Write initial condition to files
    for (int ii = 0; ii < this->numLines; ii = ii + 1)
        this->pLines[ii]->WriteOut(start_time);
    for (int ii = 0; ii < this->numBodies; ii = ii + 1)
        this->pBodies[ii]->WriteOut(start_time);

    std::cout << "Updating system..." << std::endl;
    // Save first data
    this->UpdateSystem();
    std::cout << "  System updated!" << std::endl;
}

void Simulation::LoadCase()
{
    fs::remove_all(outputFolderPath);
    fs::create_directory(outputFolderPath);

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
    std::cout << "--> Reading BCPs Properties (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosBCPs.dat");
    FILE *file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: datosBCPs.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
    }

    // Read data
    fscanf(file_pointer, "%d %[^\n]\n", &numFairBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numAnchorBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numJointBcps, bufferLine);
    fscanf(file_pointer, "%d %[^\n]\n", &numBodyBcps, bufferLine);
    numBcps = numFairBcps + numAnchorBcps + numJointBcps + numBodyBcps;

    // printf("Number of fairleads: %d\n", numFairBcps);
    // printf("Number of anchor: %d\n", numAnchorBcps);
    // printf("Number of joint: %d\n", numJointBcps);
    // printf("Number of body: %d\n", numBodyBcps);

    // Allocate a vector of BCP pointers
    pBcps = new BCP *[numBcps];
    int bcp_count = 0;

    // Read each BCP
    pFairleadBcps = new FairleadBCP *[numFairBcps];
    for (int ii = 0; ii < numFairBcps; ii++)
    {
        pFairleadBcps[ii] = new FairleadBCP(bcp_count);
        pFairleadBcps[ii]->ReadPropertiesASCII(file_pointer, inputFolderPath);
        dynamic_cast<FairleadBCP *>(pFairleadBcps[ii])->Initialize(inputFolderPath);
        pBcps[bcp_count] = pFairleadBcps[ii];
        bcp_count++;
    }
    pAnchorBcps = new AnchorBCP *[numAnchorBcps];
    for (int ii = 0; ii < numAnchorBcps; ii++)
    {
        pAnchorBcps[ii] = new AnchorBCP(bcp_count);
        pAnchorBcps[ii]->ReadPropertiesASCII(file_pointer);
        pBcps[bcp_count] = pAnchorBcps[ii];
        bcp_count++;
    }
    pJointBcps = new JointBCP *[numJointBcps];
    for (int ii = 0; ii < numJointBcps; ii++)
    {
        pJointBcps[ii] = new JointBCP(bcp_count);
        pJointBcps[ii]->ReadPropertiesASCII(file_pointer);
        dynamic_cast<JointBCP *>(pJointBcps[ii])->Initialize(this->gravity, this->waterDensity, this->waterDepth);
        pBcps[bcp_count] = pJointBcps[ii];
        bcp_count++;
    }
    pBodyBcps = new BodyBCP *[numBodyBcps];
    for (int ii = 0; ii < numBodyBcps; ii++)
    {
        pBodyBcps[ii] = new BodyBCP(bcp_count);
        pBodyBcps[ii]->ReadPropertiesASCII(file_pointer);
        pBcps[bcp_count] = pBodyBcps[ii];
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
    std::cout << "--> BCPs Properties Read" << std::endl;
}

void Simulation::ReadBcpsHDF5()
{
    std::cout << "--> Reading BCPs Properties (ASCII format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadBcpsHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> BCPs Properties Read" << std::endl;
}

void Simulation::ReadBodies()
{
    (this->*pReadBodies)();
}

void Simulation::ReadBodiesASCII()
{
    std::cout << "--> Reading Bodies Properties (ASCII format)" << std::endl;

    // Declare local variables
    int body_count = 0;
    char buffer_line[1000];
    int diff_count = 0;
    int hydro_database_count = 0;
    std::string hydro_databases_name[300];
    int max_num_bodies_database = 100;
    int pos_body = 0;
    int pos_database = 0;

    // Parse file in order to guess the number of bodies
    std::cout << "Parsing file: datosBodies.dat" << std::endl;
    std::string file_path = JoinPath(inputFolderPath, "datosBodies.dat");
    this->numBodies = parse_file(file_path);
    std::cout << "Number of bodies: " << this->numBodies << std::endl;

    if (numBodies > 0)
    {
        // Open file
        std::cout << "Opening file: datosBodies.dat" << std::endl;
        FILE *pFile = fopen(file_path.c_str(), "r");
        if (pFile == NULL)
        {
            std::stringstream ss;
            ss << "Not possible to open the file: datosBodies.dat\n    ->Dir: " << inputFolderPath << std::endl;
            throw IOError(ss.str());
        }

        // Read all bodies
        pBodies = new Body *[numBodies];
        for (int ii = 0; ii < numBodies; ii++)
        {
            std::cout << "  ... Including Body: " << ii + 1 << std::endl;
            // Discard header lines and check for body type
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
                pBodies[ii]->OpenOutputFilesASCII(outputFolderPath);
            }
            else
            {
                std::stringstream ss;
                ss << "Error while parsing file: datosBodies.dat\n --> Expected body: " << ii << " type definition\n";
                throw ValueError(ss.str());
            }
        }
        // Close file
        fclose(pFile);
        std::cout << "--> ... done!" << std::endl;

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

        std::cout << "--> HDB files considered:" << std::endl;
        for (int ii = 0; ii < hydro_database_count + 1; ii++)
        {
            std::cout << "    -->" << hydro_databases_name[ii].c_str() << std::endl;
        }

        // Arrange all the bodies by database
        Body **pBodiesSort = new Body *[numBodies];
        int *pBody_found = new int[numBodies];
        for (int ii = 0; ii < numBodies; ii++)
        {
            pBody_found[ii] = 0;
        }
        for (int ii = 0; ii < hydro_database_count; ii++)
        {
            for (int jj = 0; jj < numBodies; jj++)
            {
                if ((hydro_databases_name[ii].compare(pBodies[jj]->hydroDatabaseName) == 0) && (pBody_found[jj] == 0))
                {
                    pBody_found[jj] = 1;
                    pBodiesSort[body_count] = pBodies[jj];
                    body_count++;
                }
            }
        }
        delete[] pBody_found;
        for (int ii = 0; ii < numBodies; ii++)
        {
            pBodies[ii] = pBodiesSort[ii];
        }
        delete[] pBodiesSort;

        // Checking free bodies
        std::cout << "Checking for free bodies ..." << std::endl;
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
        pBodiesLock = new Body *[numBodiesLock];
        pBodiesFree = new Body *[numBodiesFree];
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

        // Create an array in order to store the indexes of the bodies in each database
        int **check_hydro_bodies_id = new int *[hydro_database_count];
        Body ***check_hydro_bodies = new Body **[hydro_database_count];
        for (int ii = 0; ii < hydro_database_count; ii++)
        {
            check_hydro_bodies_id[ii] = new int[max_num_bodies_database + 1];
            check_hydro_bodies[ii] = new Body *[max_num_bodies_database];
        }
        for (int ii = 0; ii < hydro_database_count; ii++)
        {
            for (int jj = 0; jj < max_num_bodies_database + 1; jj++)
            {
                check_hydro_bodies_id[ii][jj] = 0;
            }
        }

        // Check if there is some repeated body definition in each database
        std::cout << "Looking for body repetition..." << std::endl;
        for (int ii = 0; ii < this->numBodies; ii++)
        {
            // Find database position inside the array of names generated previously
            pos_database = 0;
            while (true)
            {
                if (pBodies[ii]->hydroDatabaseName.compare(hydro_databases_name[pos_database]) == 0)
                {
                    break;
                }
                pos_database++;
                if (pos_database >= hydro_database_count)
                {
                    std::stringstream ss;
                    ss << "It is not possible to find the name of the hydro database: " << pBodies[ii]->hydroDatabaseName;
                    ss << " in the list of the hydrodatabase names done with bodies definition.";
                    throw ValueError(ss.str());
                }
            }

            // Check if the body id already exist
            for (int jj = 1; jj <= check_hydro_bodies_id[pos_database][0]; jj++)
            {
                if (check_hydro_bodies_id[pos_database][jj] == pBodies[ii]->hydroDatabaseIndex)
                {
                    std::stringstream ss;
                    ss << "Repeated Hydrodynamic Bodoy Index(" << pBodies[ii]->hydroDatabaseIndex << ") definition for Body: ";
                    ss << ii << " and hydrodynamic database: " << hydro_databases_name[pos_database];
                    throw ValueError(ss.str());
                }
            }
            std::cout << check_hydro_bodies_id[pos_database][0] << std::endl;
            check_hydro_bodies_id[pos_database][0]++;
            check_hydro_bodies_id[pos_database][check_hydro_bodies_id[pos_database][0]] = pBodies[ii]->hydroDatabaseIndex;
            check_hydro_bodies[pos_database][check_hydro_bodies_id[pos_database][0] - 1] = pBodies[ii];
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

            // Set database to the target Body object
            hydro_file_path = JoinPath(this->inputFolderPath, hydro_databases_name[pos_database]);
            this->pBodies[ii]->LoadHydrodynamicDatabase(check_hydro_bodies[pos_database]);
        }

        // Fill System Matrix
        std::cout << "Fill system matrix...\n";
        arma::span a1;
        arma::span a2;
        int nn2;
        int db_shift;
        int body_shift;
        this->pSystemMatrix = new arma::mat(6 * this->numBodies, 6 * this->numBodies, arma::fill::zeros);
        this->pSystemMatrixInv = new arma::mat(6 * this->numBodies, 6 * this->numBodies, arma::fill::zeros);
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

            db_shift = 0;
            for (int jj = 0; jj < pos_database; jj++)
            {
                db_shift += 6 * check_hydro_bodies_id[pos_database][0];
            }

            // Look for position of the body
            body_shift = 6 * pBodies[ii]->hydroDatabaseIndex;

            // Fill system matrix
            a1 = arma::span(db_shift + body_shift, db_shift + body_shift + 5);
            a2 = arma::span(db_shift, db_shift + 6 * check_hydro_bodies_id[pos_database][0] - 1);
            std::cout << "  ... Filling system matrix for body: " << ii + 1 << std::endl;
            (*pSystemMatrix)(a1, a2) += pBodies[ii]->pHydro->GetTotalMass();
            std::cout << "  ... done!" << std::endl;
            pBodies[ii]->sysMatSpan1 = a1;
            pBodies[ii]->sysMatSpan2 = a2;
            pBodies[ii]->sysMatInd1 = arma::regspace<arma::uvec>(db_shift + body_shift, db_shift + body_shift + 5);
        }

        // Take free dofs index vectors from the system matrix
        sysMatIndFree = arma::zeros<arma::uvec>(6 * numBodiesFree);
        std::cout << "Take free dofs index vectors from the system matrix...\n";
        for (int ii = 0; ii < this->numBodiesFree; ii++)
        {
            sysMatIndFree(arma::span(6 * ii, 6 * (ii + 1) - 1)) = pBodiesFree[ii]->sysMatInd1;
        }
        sysMatIndLock = arma::zeros<arma::uvec>(6 * numBodiesLock);
        std::cout << "Take locked dofs index from matrix...\n";
        for (int ii = 0; ii < this->numBodiesLock; ii++)
        {
            sysMatIndLock(arma::span(6 * ii, 6 * (ii + 1) - 1)) = pBodiesLock[ii]->sysMatInd1;
        }

        // Extract submatrices from system matrix
        std::cout << "Extract submatrices from system matrix...\n";
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
        // Invert system matrix
        std::cout << "Inverting system matrix...\n";
        *pSystemMatrixInv = arma::solve(*pSystemMatrix, eye(size(*pSystemMatrix)));
        if (numBodiesFree > 0)
        {
            pSystemMatrixFFInv = new arma::mat(6 * numBodiesFree, 6 * numBodiesFree, arma::fill::zeros);
            *pSystemMatrixFFInv = arma::solve(*pSystemMatrixFF, eye(size(*pSystemMatrixFF)));
        }
        std::cout << "System matrix inverted...\n";

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
        // TODO: Check memory leaks through the code
        delete[] check_hydro_bodies_id;
        delete[] check_hydro_bodies;

        std::cout << "--> Bodies Properties Read" << std::endl;
    }
}

void Simulation::ReadBodiesHDF5()
{
    std::cout << "--> Reading Bodies Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadBodiesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> Bodies Properties Read" << std::endl;
}

void Simulation::ReadLines()
{
    (this->*pReadLines)();
}

void Simulation::ReadLinesASCII()
{
    std::cout << "--> Reading Lines Properties (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosLines.dat");
    FILE *file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: datosLines.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
    }

    // Read total number of springs to read
    fscanf(file_pointer, "%d %[^\n]\n", &numLines, bufferLine);

    // Allocate a vector of Line pointers
    pLines = new Line *[numLines];

    // Read all lines
    for (int ii = 0; ii < numLines; ii++)
    {
        pLines[ii] = new Line(ii, gravity, waterDensity, waterDepth);
        // TODO: This try-catch block should be removed and the exceptions should be handled in the Line class
        try
        {
            pLines[ii]->ReadPropertiesASCII(file_pointer);
            pLines[ii]->OpenOutputFilesASCII(outputFolderPath);
        }
        catch (int e)
        {
            if (e == 0)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " is under the floor level." << std::endl
                          << std::endl;
            if (e == 1)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor and it shouldn't. " << std::endl
                          << std::endl;
            if (e == 2)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " is not tense and laying on the seafloor. It should be pretensed. " << std::endl
                          << std::endl;
            if (e == 3)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " is not tense and vertical. It should be pretensed. " << std::endl
                          << std::endl;
            if (e == 4)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " is not tense and it should. " << std::endl
                          << std::endl;
            if (e == 5)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " initial shape can't be computed with QS method. " << std::endl;
            if (e == 6)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor althoug none of its ends are there. " << std::endl
                          << std::endl;
        }
    }

    // Close the file
    fclose(file_pointer);
    std::cout << "--> Lines Properties Read" << std::endl;
}

void Simulation::ReadLinesHDF5()
{
    std::cout << "--> Reading Lines Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadLinesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> Lines Properties Read" << std::endl;
}

void Simulation::ReadSinking()
{
    (this->*pReadSinking)();
}

void Simulation::ReadSinkingASCII()
{
    std::cout << "--> Reading Sinking Properties (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine[1000];
    int sinkingBodyIndex;

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataSinking.dat");
    FILE *pFile = fopen(file_path.c_str(), "r");
    if (pFile == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: dataSinking.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
    }

    // Read total number of sinking bodies to read
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
        pSinking = new Sinking *[numSinking];
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
            pSinking[ii]->OpenOutputFilesASCII(outputFolderPath);
        }
    }

    // Close file
    fclose(pFile);

    std::cout << "--> Sinking Properties Read" << std::endl;
}

void Simulation::ReadSinkingHDF5()
{
    std::cout << "--> Reading Sinking Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadSinkingHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> Sinking Properties Read" << std::endl;
}

void Simulation::ReadSprings()
{
    (this->*pReadSprings)();
}

void Simulation::ReadSpringsASCII()
{
    std::cout << "--> Reading Spring Properties (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosSprings.dat");
    FILE *file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: datosSprings.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
    }

    // Read file contents
    fscanf(file_pointer, "%d %[^\n]", &numSprings, bufferLine);

    // Allocate a vector of Spring pointers
    pSprings = new Spring *[numSprings];

    // Read each Spring
    for (int ii = 0; ii < numSprings; ii++)
    {
        pSprings[ii] = new Spring(ii);
        pSprings[ii]->ReadPropertiesASCII(file_path);
    }

    // Close the file
    fclose(file_pointer);

    std::cout << "--> Spring Properties Read" << std::endl;
}

void Simulation::ReadSpringsHDF5()
{
    std::cout << "--> Reading Spring Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadSpringsHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> Spring Properties Read" << std::endl;
}

void Simulation::ReadProperties()
{
    (this->*pReadProperties)();
}

void Simulation::ReadPropertiesASCII()
{
    std::cout << "--> Reading Simulation Properties (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine[1000];
    int dummyBool;

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosProblema.dat");
    FILE *file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: datosProblema.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
    }

    // Read data
    fscanf(file_pointer, "%lf %[^\n]\n", &gravity, bufferLine);                // Gravity acceleration [m/s^2]
    fscanf(file_pointer, "%lf %[^\n]\n", &waterDensity, bufferLine);           // Water density [kg/m^3]
    fscanf(file_pointer, "%lf %[^\n]\n", &airAtmPresDensity, bufferLine);      // Air density at atmospheric pressure [kg/m^3]
    fscanf(file_pointer, "%lf %[^\n]\n", &airAtmPres, bufferLine);             // Atmospheric pressure [Pa]
    fscanf(file_pointer, "%lf %[^\n]\n", &airAdiabaticDilation, bufferLine);   // Air adiabatic dilation [-]
    fscanf(file_pointer, "%lf %[^\n]\n", &waterDepth, bufferLine);             // Seabed vertical coordinate [m]
    fscanf(file_pointer, "%lf %[^\n]\n", &writeTimeStep, bufferLine);          // Output time step [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &maxTimeStep, bufferLine);            // Maximum time step for time integration [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &hydroTimeStep, bufferLine);          // Time step for hydrodynamic forces computation [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &fastTimeStep, bufferLine);           // Time step for FAST wind turbines forces computation [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &fastControllerTimeStep, bufferLine); // Time step for FAST wind turbines controller update [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &timeIRF, bufferLine);                // IRF time [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &sinkingTimeStep, bufferLine);        // Time step for synking hydrodinamic data bases update [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &controllerTimeStep, bufferLine);     // Time step for winches controller [s]
    fscanf(file_pointer, "%lf %[^\n]\n", &simulationTime, bufferLine);         // Total time of simulation [s]
    fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    rotSimpFlag = dummyBool;                                         // Flag to use simplification for rigid body rotation dynamics [0 No, 1 Yes]
    fscanf(file_pointer, "%d %[^\n]\n", &timeIntMethod, bufferLine); // Temporal integration alforithm [1: BDF1, 2: BDFN]
    fscanf(file_pointer, "%d %[^\n]\n", &timeIntOrder, bufferLine);  // Order for temporal integration
    fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    timeIntAdaptivity = dummyBool;                                    // Time step adaptivity [0 No, 1 Yes]
    fscanf(file_pointer, "%lf %[^\n]\n", &timeIntAbsTol, bufferLine); // Absolute tolerance for temporal integration.
    fscanf(file_pointer, "%lf %[^\n]\n", &timeIntRelTol, bufferLine); // Relative tolerance for temporal integration.
    fscanf(file_pointer, "%d %[^\n]\n", &maxIterStep, bufferLine);    // Maximum number of iterations for one step of temporal integration.
    fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    readEquilibrium = dummyBool; // Read Equilibrio.dat? [0 No, 1 Yes]
    fscanf(file_pointer, "%d %[^\n]\n", &dummyBool, bufferLine);
    writeEquilibrium = dummyBool;                                 // Write Equilibrio.dat? [0 No, 1 Yes]
    fscanf(file_pointer, "%d %[^\n]\n", &flagStatic, bufferLine); // Mooring initial condition flag [0: Catenary, 1: Newton's]

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

    std::cout << "--> Simulation Properties Read" << std::endl;
}

void Simulation::ReadPropertiesHDF5()
{
    std::cout << "--> Reading Simulation Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadPropertiesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> Simulation Properties Read" << std::endl;
}

void Simulation::ReadWaves()
{
    (this->*pReadWaves)();
}

void Simulation::ReadWavesASCII()
{
    std::cout << "--> Reading Waves (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine[1000];
    char wave_type[1000];
    double H, T, D;

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataWaves.dat");
    FILE *file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: dataWaves.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
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

    // Create wave object and read additional parameters if necessary
    if (strncmp(wave_type, "REG", 3) == 0)
    {
        pWave = new RegularWave(this, H, T, D);
    }
    else if (strncmp(wave_type, "IRR", 3) == 0)
    {
        pWave = new IrregularWave(this, H, T, D);
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

    std::cout << "--> ... wave read!" << std::endl;

    // Preprocess wave data if H>0
    if (H > 0.0)
    {
        std::cout << "-->  Preprocessing wave..." << std::endl;
        pWave->CheckBreakingWave();
        pWave->GetWaveSpectrum();
        if (strncmp(wave_type, "IRR", 3) == 0)
        {
            pWave->WriteOut(outputFolderPath);
        }
    }
    else
    {
        std::cout << "-->  Preprocessing dummy wave..." << std::endl;
        pWave->SetZeroHeight();
    }

    std::cout << "--> ... wave preprocessed!" << std::endl;
}

void Simulation::ReadWavesHDF5()
{
    std::cout << "--> Reading Waves (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWavesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> Simulation Properties Read" << std::endl;
}

void Simulation::ReadWinches()
{
    (this->*pReadWinches)();
}

void Simulation::ReadWinchesASCII()
{
    std::cout << "--> Reading Winches Properties (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosWinchies.dat");
    FILE *file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: datosWinchies.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
    }

    // Read number of winches defined in the file
    fscanf(file_pointer, "%d %[^\n]\n", &numWinches, bufferLine);
    printf("NumWinches: %d - UseWinches: %d\n", numWinches, useWinches);
    if ((numWinches == 0) && useWinches)
    {
        throw ValueError("Use of winches is requested when loading BCPs but there is no winches specified in datosWinches.dat\n");
    }

    // Allocate a vector of pointers to Winch class objects
    pWinches = new Winchie *[numWinches];

    // Read Winches
    for (int ii = 0; ii < numWinches; ii++)
    {
        pWinches[ii] = new Winchie(ii);
        pWinches[ii]->ReadPropertiesASCII(file_pointer);
        pWinches[ii]->LineW = pLines[pWinches[ii]->nLine - 1];
    }

    // Close the file
    fclose(file_pointer);
    std::cout << "--> Winches Properties Read" << std::endl;

    // Read Winches Controller
    std::cout << "--> Reading Winches Controller Properties (ASCII format)" << std::endl;
    file_path = JoinPath(inputFolderPath, "datosWinchiesController.dat");
    file_pointer = fopen(file_path.c_str(), "r");
    WinchesController = WinchieController(numWinches, pWinches, this);
    WinchesController.ReadPropertiesASCII(file_pointer);
    fclose(file_pointer);
    // Open Winches Controller output files
    WinchesController.OpenOutputFilesASCII(outputFolderPath);
    std::cout << "--> Winches Controller Properties Read" << std::endl;
}

void Simulation::ReadWinchesHDF5()
{
    std::cout << "--> Reading Winches Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWinchesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> Winches Properties Read" << std::endl;
}

void Simulation::ReadSeaFloor()
{
    (this->*pReadSeaFloor)();
}

void Simulation::ReadSeaFloorASCII()
{
    std::cout << "--> Reading SeaFloor (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataSeaFloor.dat");
    FILE *file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: dataSeaFloor.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
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
    pSeaFloor = new SeaFloor *[numFloor];
    int floor_count = 0;

    // Read all SeaFloor objects
    pBathymetry = new Bathymetry *[numBathymetry];
    for (int ii = 0; ii < numBathymetry; ii++)
    {
        pBathymetry[ii] = new Bathymetry(floor_count);
        pBathymetry[ii]->ReadPropertiesASCII(file_pointer, inputFolderPath);
        pSeaFloor[floor_count] = pBathymetry[ii];
        floor_count++;
    }
    pInclined = new Inclined *[numInclined];
    for (int ii = 0; ii < numInclined; ii++)
    {
        pInclined[ii] = new Inclined(floor_count);
        pInclined[ii]->ReadPropertiesASCII(file_pointer);
        pSeaFloor[floor_count] = pInclined[ii];
        floor_count++;
    }
    pFlat = new Flat *[numFlat];
    for (int ii = 0; ii < numFlat; ii++)
    {
        pFlat[ii] = new Flat(floor_count);
        pFlat[ii]->ReadPropertiesASCII(file_pointer);
        pSeaFloor[floor_count] = pFlat[ii];
        floor_count++;
    }

    // Close file
    fclose(file_pointer);

    std::cout << "--> Floor Properties Read" << std::endl;
}

void Simulation::ReadSeaFloorHDF5()
{
    std::cout << "--> Reading SeaFloor Properties (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadSeaFloorHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> SeaFloor Properties Read" << std::endl;
}

void Simulation::ReadWindTurbines(void)
{
    (this->*pReadWindTurbines)();
}

void Simulation::ReadWindTurbinesASCII(void)
{
    std::cout << "--> Reading Wind Turbines Properties (ASCII format)" << std::endl;

    // Declare local variables
    char bufferLine[1000];

    // Open file
    std::string file_path = JoinPath(inputFolderPath, "datosWindTurbines.dat");
    FILE *file_pointer = fopen(file_path.c_str(), "r");
    if (file_pointer == NULL)
    {
        std::stringstream ss;
        ss << "Not possible to open the file: datosWindTurbines.dat\n    ->Dir: " << inputFolderPath << std::endl;
        throw IOError(ss.str());
    }

    // Read number of Wind Turbines defined in the file
    fscanf(file_pointer, "%d %[^\n]\n", &numWindTurbines, bufferLine);

    // Allocate a vector of pointers to WindTurbine class objects
    pWindTurbines = new WindTurbine *[numWindTurbines];
    for (int ii = 0; ii < numWindTurbines; ii++)
    {
        pWindTurbines[ii] = new WindTurbine(ii, this);
        pWindTurbines[ii]->ReadPropertiesASCII(file_pointer);
    }

    // Close the file
    fclose(file_pointer);

    std::cout << "--> Wind Turbines Properties Read" << std::endl;
}

void Simulation::ReadWindTurbinesHDF5(void)
{
    std::cout << "--> Reading Wind Turbines (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadWindTurbinesHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> Wind Turbines Properties Read" << std::endl;
}

void Simulation::ReadOWCs(void)
{
    (this->*pReadOWCs)();
}

void Simulation::ReadOWCsASCII(void)
{
    std::cout << "--> Reading OWCs Properties (ASCII format)" << std::endl;
    // Declare local variables
    char bufferLine[1000];
    // Open file
    std::string file_path = JoinPath(inputFolderPath, "dataOWCs.dat");
    FILE *pFile = fopen(file_path.c_str(), "r");
    if (pFile == NULL)
    {
        std::cout << "    --> WARNING: dataOWCs.dat was not found! Setting numOWCs = 0!" << std::endl;
        numOWCs = 0;
    }
    else
    {
        // Read total number of owcs
        fscanf(pFile, "%d %[^\n]\n", &numOWCs, bufferLine);
        if (numOWCs > 0)
        {
            // Initiallice OWC array
            pOWCs = new OWC *[numOWCs];
            for (int ii = 0; ii < numOWCs; ii++)
            {
                pOWCs[ii] = new OWC(ii, this);
                pOWCs[ii]->Initialize(pFile);
            }
        }
        // Close file
        fclose(pFile);
    }
    std::cout << "--> OWCs Properties Read" << std::endl;
}

void Simulation::ReadOWCsHDF5(void)
{
    std::cout << "--> Reading OWCs (HDF5 format)" << std::endl;
    std::stringstream ss;
    ss << "Method ReadOWCsHDF5 in class Simulation not implemented yet.";
    throw NotImplementedError(ss.str());
    std::cout << "--> OWCs Properties Read" << std::endl;
}

void Simulation::Run()
{
    std::cout << "--> Starting Simulation Run..." << std::endl;

    // Declare local variables
    time_t tstart, tend;
    double wallTime = 0.0;
    double wallTimeHydro = 0.0;
    double wallTimeFAST = 0.0;
    double wallTimeControllerFAST = 0.0;
    double wallTimeSinking = 0.0;
    double wallTimeController = 0.0;
    tstart = time(0);
    bool flag_debug_lines = true;

    std::cout << "    t = " << wallTime << " s" << std::endl;

    // TODO: Implement a logger with different levels of verbosity
    //  std::cout<< "In Simulation::Run --> Starting temporal integration loop "<< std::endl;
    do
    {
        // std::cout<< "In Simulation::Run --> step() "<< std::endl;
        pTimeSolver->step();

        // Print lines in first time step for debugging
        if (flag_debug_lines)
        {
            flag_debug_lines = false;
            for (int ii = 0; ii < numLines; ii = ii + 1)
            {
                fprintf(pLines[ii]->pfile_line_debug, "s    x    y    z \n");
                for (int jj = 0; jj < pLines[ii]->N; jj = jj + 1)
                    fprintf(pLines[ii]->pfile_line_debug, "%f    %f    %f    %f \n",
                            pLines[ii]->s(jj, 0),
                            pLines[ii]->pos(jj, 0),
                            pLines[ii]->pos(jj, 1),
                            pLines[ii]->pos(jj, 2));
                fclose(pLines[ii]->pfile_line_debug);
            }
        }

        if (pTimeSolver->t >= wallTime + writeTimeStep)
        {
            // std::cout<< "In Simulation::Run --> WriteOut() "<< std::endl;
            wallTime = wallTime + writeTimeStep;
            std::cout << "    t = " << wallTime << " s" << std::endl;
            for (int ii = 0; ii < numLines; ii = ii + 1)
                pLines[ii]->WriteOut(wallTime);
            for (int ii = 0; ii < numBodies; ii = ii + 1)
                pBodies[ii]->WriteOut(wallTime);
            for (int ii = 0; ii < numOWCs; ii = ii + 1)
                pOWCs[ii]->WriteOut(wallTime);
        }

        if (pTimeSolver->t >= wallTimeHydro + hydroTimeStep)
        {
            // std::cout<< "In Simulation::Run --> Computing hydrodynamic forces... "<< std::endl;
            wallTimeHydro += hydroTimeStep;
            if (numBodies > 0)
            {
                UpdateSystem();
            }
            for (int ii = 0; ii < numBodies; ii = ii + 1)
            {
                pBodies[ii]->Fb = pBodies[ii]->pHydro->CalculateHydrodynamicForces(wallTimeHydro);
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
            if (pTimeSolver->t >= wallTimeController + controllerTimeStep)
            {
                // std::cout << "In Simulation::Run --> Controlling winches... " << std::endl;
                wallTimeController += controllerTimeStep;
                WinchesController.controlWinchies(wallTimeController);
                WinchesController.WriteOut(wallTimeController);
                // std::cout << "In Simulation::Run --> ... done controlling winches!" << std::endl;
            }
        }

    } while (pTimeSolver->t <= simulationTime);

    tend = time(0);
    double computational_time = difftime(tend, tstart);
    if (computational_time < 60)
    {
        std::cout << "    Computational time  : " << computational_time << " seconds" << std::endl;
    }
    else if (computational_time < 3600)
    {
        std::cout << "    Computational time  : " << computational_time / 60 << " minutes" << std::endl;
    }
    else
    {
        std::cout << "    Computational time  : " << computational_time / 3600 << " hours" << std::endl;
    }
    std::cout << "    Total function calls: " << numCallsSysFun << std::endl;
    std::cout << "    Total jac calls: " << pTimeSolver->iJ << std::endl
              << std::endl;

    if (writeEquilibrium == 1)
    {
        std::cout << "  Writting data to Equilibrio.dat ..." << std::endl
                  << std::endl;
        std::string filename = JoinPath(outputFolderPath, "Equilibrio.dat");
        pTimeSolver->y.save(filename, arma::arma_ascii);
    }
}

void Simulation::SetupCase()
{
    std::cout << "--> Setting up the case configuration ..." << std::endl;

    // Count the number of BCP in each body and create pointer array
    if (numBodies > 0)
    {
        std::cout << "  --> Counting the number of BCPs in each body ..." << std::endl;
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
                   << " datosBCPs.dat";
                throw ValueError(ss.str());
            }
            pBcps[pBodies[ii]->pIndexBcps[jj]]->numBodiesBcp++;
        }
    }
    bool *pDefined_body_bcps = new bool[numBcps];
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
                pBcps[pBodies[ii]->pIndexBcps[jj]]->pBodies = new Body *[pBcps[pBodies[ii]->pIndexBcps[jj]]->numBodiesBcp];
                pDefined_body_bcps[pBodies[ii]->pIndexBcps[jj]] = true;
            }
        }
    }
    delete[] pDefined_body_bcps;
    if (numBodies > 0)
    {
        std::cout << "  --> ... done!" << std::endl;
    }

    // Floor triangulation initialization
    if (numFloor > 0)
    {
        std::cout << "  --> Initializing floor..." << std::endl;
    }
    for (int ii = 0; ii < numFloor; ii++)
    {
        if (pSeaFloor[ii]->GetType() == 3)
        {
            dynamic_cast<Bathymetry *>(pSeaFloor[ii])->getVertexNormals();
            dynamic_cast<Bathymetry *>(pSeaFloor[ii])->getProjectionMatrix();
        }
        else if (pSeaFloor[ii]->GetType() == 2)
        {
            dynamic_cast<Inclined *>(pSeaFloor[ii])->getPlaneEquation();
        }
        else if (pSeaFloor[ii]->GetType() == 1)
        {
            std::cout << "    --> Floor level: " << std::endl
                      << pSeaFloor[ii]->fondo << " m" << std::endl;
        }
    }
    if (numFloor > 0)
    {
        std::cout << "  --> ... done!" << std::endl;
    }

    // Assing to each BCP the corresponding Body pointer
    if (numBodies > 0)
    {
        std::cout << "  --> Assigning to each BCP the corresponding body  ..." << std::endl;
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        std::cout << "    --> Body " << ii + 1 << " pos = " << pBodies[ii]->pos.t() << std::endl;
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
        std::cout << "  --> ... done!" << std::endl;
    }

    // Count the number of Wind Turbines in each body and create pointer array
    if (numBodies > 0)
    {
        std::cout << "  --> Checking the number of Wind Turbines in each body ..." << std::endl;
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        for (int jj = 0; jj < pBodies[ii]->numWindTurbs; jj++)
        {
            if (pBodies[ii]->pIndexWindTurbs[jj] + 1 > numWindTurbines)
            {
                std::stringstream ss;
                ss << "Wind Turbine index: " << pBodies[ii]->pIndexWindTurbs[jj] << " in Body: " << pBodies[ii]->GetId()
                   << " is out of range when compare with the Number of Wind Turbines (" << numWindTurbines << ") defined in"
                   << " datosWindTurbines.dat";
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
        std::cout << "  --> ... done!" << std::endl;
    }

    // Assing to each Body the corresponding Wind Turbine pointers
    if (numBodies > 0)
    {
        std::cout << "       Assingning to each Body the corresponding Wind Turbine pointers  ..." << std::endl;
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
        std::cout << "  --> ... done!" << std::endl;
    }

    // Count the number of lines in each joint BCP
    if (numBcps > 0)
    {
        std::cout << "  --> Counting the number of lines in each joint BCP ..." << std::endl;
    }
    for (int ii = 0; ii < numBcps; ii++)
    {
        if (pBcps[ii]->GetType() == 3)
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
        std::cout << "  --> ... done!" << std::endl;
    }

    // Count the number of Lines in each body and create pointer array
    if (numLines > 0)
    {
        std::cout << "        Counting the number of Lines in each body ..." << std::endl;
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
                   << " datosBCPs.dat";
                throw ValueError(ss.str());
            }
            pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp++;
        }
    }
    bool *pDefined_lines_bcps = new bool[numBcps];
    for (int ii = 0; ii < numBcps; ii++)
    {
        pDefined_lines_bcps[ii] = 0;
    }
    for (int ii = 0; ii < numLines; ii++)
    {
        for (int jj = 0; jj < pLines[ii]->numBcps; jj++)
        {
            if (!pDefined_lines_bcps[pLines[ii]->indexBcps[jj]])
            {
                pBcps[pLines[ii]->indexBcps[jj]]->pLines = new Line *[pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp];
                pBcps[pLines[ii]->indexBcps[jj]]->pBcpLineIndex = new int[pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp];
                pBcps[pLines[ii]->indexBcps[jj]]->pBcpLineNode = new int[pBcps[pLines[ii]->indexBcps[jj]]->numLinesBcp];
                pDefined_lines_bcps[pLines[ii]->indexBcps[jj]] = true;
            }
        }
    }
    delete[] pDefined_lines_bcps;
    if (numLines > 0)
    {
        std::cout << "  --> ... done!" << std::endl;
    }

    // Assing to each Line the corresponding BCP pointer
    if (numLines > 0)
    {
        std::cout << "        Assigning to each Line the corresponding BCP ..." << std::endl;
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

        if (pLines[ii]->pLineBcps[0]->GetType() != 3 && pLines[ii]->pLineBcps[1]->GetType() != 3)
        {
            numDofLinesTotal += (pLines[ii]->N - 2);
            pLines[ii]->first_node = 1;
            pLines[ii]->last_node = pLines[ii]->N - 1;
        }
        else if (pLines[ii]->pLineBcps[0]->GetType() == 3 && pLines[ii]->pLineBcps[1]->GetType() == 3)
        {
            numDofLinesTotal += pLines[ii]->N;
            pLines[ii]->first_node = 0;
            pLines[ii]->last_node = pLines[ii]->N;
        }
        else if (pLines[ii]->pLineBcps[0]->GetType() == 3 && pLines[ii]->pLineBcps[1]->GetType() != 3)
        {
            numDofLinesTotal += (pLines[ii]->N - 1);
            pLines[ii]->first_node = 0;
            pLines[ii]->last_node = pLines[ii]->N - 1;
        }
        else if (pLines[ii]->pLineBcps[0]->GetType() != 3 && pLines[ii]->pLineBcps[1]->GetType() == 3)
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
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " is under the floor level." << std::endl
                          << std::endl;
            if (e == 1)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor and it shouldn't. " << std::endl
                          << std::endl;
            if (e == 2)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " is not tense and laying on the seafloor. It should be pretensed. " << std::endl
                          << std::endl;
            if (e == 3)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " is not tense and vertical. It should be pretensed. " << std::endl
                          << std::endl;
            if (e == 4)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " is not tense and it should. " << std::endl
                          << std::endl;
            if (e == 5)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " initial shape can't be computed with QS method. " << std::endl;
            if (e == 6)
                std::cout << "ERROR: Line " << pLines[ii]->nLine << " touches the seafloor althoug none of its ends are there. " << std::endl
                          << std::endl;
        }

        pLines[ii]->pLineSeaFloor = pSeaFloor[pLines[ii]->indexSeaFloor];
    }
    if (numLines > 0)
    {
        std::cout << "  --> ... done!" << std::endl;
    }

    // Count the number of line nodes without repetition of joint nodes
    if (numLines > 0)
    {
        std::cout << "        Counting the number of line nodes without repetition of joint nodes ..." << std::endl;
    }
    int numUsedJointBCPs = 0;
    for (int jj = 0; jj < numLines; jj++)
    {
        numAllLinesNodes += pLines[jj]->N;
        for (int kk = 0; kk < 2; kk++)
        {
            if (pLines[jj]->pLineBcps[kk]->GetType() == 3)
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
        std::cout << "  --> ... done!" << std::endl;
    }

    // Build the lines coupling sparse matrix and store the lines index vectors
    if (numLines > 0)
    {
        std::cout << "  --> Building the lines coupling sparse matrix and storing the lines index vectors ..." << std::endl;
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
        std::cout << "  --> ... done!" << std::endl;
    }

    // Compute the lines coupling matrices
    if (numLines > 0)
    {
        std::cout << "  --> Computing the lines coupling matrices ..." << std::endl;
    }
    pLinesCouplingMatrix = new arma::mat(numAllLinesNodes, numAllLinesNodes, arma::fill::zeros);
    pLinesCouplingMatrixInv = new arma::mat(numAllLinesNodes, numAllLinesNodes);
    pLinesCouplingMatrix_sp = new arma::sp_mat(numAllLinesNodes, numAllLinesNodes);
    ComputeLinesCouplingMatrix();
    if (numLines > 0)
    {
        std::cout << "  --> ... done!" << std::endl;
    }
    // TODO: 100 should be a parameter
    if (numAllLinesNodes < 100)
    {
        *pLinesCouplingMatrixInv = arma::solve(*pLinesCouplingMatrix, eye(size(*pLinesCouplingMatrix)));
        std::string filename = JoinPath(outputFolderPath, "LinesCouplingMatrix.dat");
        (*pLinesCouplingMatrix).save(filename, arma::arma_ascii);
    }
    if (numLines > 0)
    {
        std::cout << "  --> ... done!" << std::endl;
    }

    // Compute equilibrium with FEM for all lines at the same time
    if (!readEquilibrium && flagStatic == 1)
    {
        std::cout << "  --> Computing the equilibrium with FEM for all lines at the same time ..." << std::endl;

        // Store the lines friction model and tension flag
        arma::umat flagLineas;
        arma::umat flagTension;
        flagLineas.zeros(numLines, 1);
        flagTension.zeros(numLines, 1);
        for (int i = 0; i < numLines; i++)
        {
            // Store the original values
            flagLineas(i) = pLines[i]->frictionModel;
            flagTension(i) = pLines[i]->flag_tension;
            // Set them to a value that allows the equilibrium computation
            pLines[i]->frictionModel = 0;
            pLines[i]->flag_tension = 1;
        }

        // Start the equilibrium computation with an initial guess
        arma::mat posicionInicial = ComputeLinesInitialPoint();

        // Perform the equilibrium computation
        ComputeLinesEquilibrium(posicionInicial);

        // Return the original values of the friction model and tension flag
        for (int i = 0; i < numLines; i++)
        {
            pLines[i]->frictionModel = flagLineas(i);
            pLines[i]->flag_tension = flagTension(i);
        }

        std::cout << "  --> ... done!" << std::endl;
    }

    // Store the initial position, required for the stick-slip friction model
    if (numLines > 0)
    {
        std::cout << "  --> Storing the initial position of the lines ..." << std::endl;
    }
    for (int i = 0; i < numLines; i++)
    {
        pLines[i]->posFriccion = pLines[i]->pLineSeaFloor->projectPoints(pLines[i]->pos)(0, 0);
    }
    if (numLines > 0)
    {
        std::cout << "  --> ... done!" << std::endl;
    }

    // Setup Springs
    if (numSprings > 0)
    {
        std::cout << "  --> Setting up springs ..." << std::endl;
    }
    for (int ii = 0; ii < numSprings; ii++)
    {
        std::cout << "        Spring: " << ii << "\n";
        std::cout << "            ->BCP_1 " << pSprings[ii]->BCP_1 << "\n";
        std::cout << "            ->BCP_2 " << pSprings[ii]->BCP_2 << "\n";
        pSprings[ii]->SpringBCP[0] = pBcps[pSprings[ii]->BCP_1];
        pSprings[ii]->SpringBCP[1] = pBcps[pSprings[ii]->BCP_2];
    }
    if (numSprings > 0)
    {
        std::cout << "  --> ... done!" << std::endl;
    }

    // Setup hidro data bases
    if (numBodies > 0)
    {
        std::cout << "  --> Setting up hydro data bases ..." << std::endl;
    }
    for (int ii = 0; ii < numBodies; ii++)
    {
        pBodies[ii]->pHydro->SetUp();
    }
    if (numBodies > 0)
    {
        std::cout << "  --> ... done!" << std::endl;
    }

    // Setup winchies controller
    if (useWinches)
    {
        std::cout << "  --> Setting up winchies controller ..." << std::endl;
        WinchesController.SetUpWinchiesController();
        std::cout << "  --> ... done!" << std::endl;
    }

    if (numWindTurbines > 0)
    {
        // Setup wind turbines
        std::cout << "  --> Setting up wind turbines ..." << std::endl;
        for (int ii = 0; ii < numWindTurbines; ii++)
        {
            pWindTurbines[ii]->Initialize();
        }
        std::cout << "  --> ... done!" << std::endl;

        // Computing bodies structural mass considering wind turbines
        std::cout << "  --> Computing bodies structural mass considering wind turbines ..." << std::endl;
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
        std::cout << "  --> ... done!" << std::endl;
    }

    std::cout << "--> ... case configuration done!" << std::endl;
}

void Simulation::ComputeLinesCouplingMatrix(void)
{
    // Compute the matrices used to couple the lines which are joined together with a Joint type BCP

    // Loop over all the lines to compute the mass matrix from the lines
    for (int jj = 0; jj < numLines; jj++)
    {
        int numLineNodes_tmp = pLines[jj]->N;

        std::cout << "    Assembling local mass matrix for line " << jj << std::endl;
        arma::mat LineMassMat_tmp = pLines[jj]->MM * pLines[jj]->dL;
        if (pLines[jj]->pLineBcps[0]->GetType() != 3)
        {
            LineMassMat_tmp.row(0) = arma::zeros(1, numLineNodes_tmp);
            LineMassMat_tmp(0, 0) = 1.0;
        }
        if (pLines[jj]->pLineBcps[1]->GetType() != 3)
        {
            LineMassMat_tmp.row(numLineNodes_tmp - 1) = arma::zeros(1, numLineNodes_tmp);
            LineMassMat_tmp(numLineNodes_tmp - 1, numLineNodes_tmp - 1) = 1.0;
        }

        // Add the mass matrix to the global matrix in sparse or full format depending on the size
        // TODO: 100 should be a parameter
        std::cout << "    Adding local mass matrix to the global matrix" << std::endl;
        if (numAllLinesNodes >= 100)
        {
            for (int irow = 0; irow < numLineNodes_tmp; irow++)
            {
                for (int icol = 0; icol < numLineNodes_tmp; icol++)
                {
                    (*pLinesCouplingMatrix_sp)(pLines[jj]->ind4CouplingMat(irow), pLines[jj]->ind4CouplingMat(icol)) += LineMassMat_tmp(irow, icol);
                }
            }
        }
        else
        {
            (*pLinesCouplingMatrix).submat(pLines[jj]->ind4CouplingMat, pLines[jj]->ind4CouplingMat) += LineMassMat_tmp;
        }
    }

    // Loop over all the BCPs to compute the mass matrix from the joints
    for (int jj = 0; jj < numBcps; jj++)
    {
        if ((pBcps[jj]->GetType() == 3) && (pBcps[jj]->flag_assigned == 1))
        {
            // TODO: 100 should be a parameter
            if (numAllLinesNodes >= 100)
            {
                (*pLinesCouplingMatrix_sp)(pBcps[jj]->couplingMatIndex, pBcps[jj]->couplingMatIndex) += pBcps[jj]->mass_Joint;
            }
            else
            {
                (*pLinesCouplingMatrix)(pBcps[jj]->couplingMatIndex, pBcps[jj]->couplingMatIndex) += pBcps[jj]->mass_Joint;
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
    else if (!incDataFormat.compare("HDF5"))
    {
        dataFormat = 1;
        dataFormatStr = "HDF5";
        inputFolderPath = incProjectPath;
        outputFolderPath = incProjectPath;
        pReadProperties = &Simulation::ReadPropertiesHDF5;
        pReadWaves = &Simulation::ReadWavesHDF5;
        pReadBcps = &Simulation::ReadBcpsHDF5;
        pReadBodies = &Simulation::ReadBodiesHDF5;
        pReadSinking = &Simulation::ReadSinkingHDF5;
        pReadLines = &Simulation::ReadLinesHDF5;
        pReadSprings = &Simulation::ReadSpringsHDF5;
        pReadWinches = &Simulation::ReadWinchesHDF5;
        pReadSeaFloor = &Simulation::ReadSeaFloorHDF5;
        pReadWindTurbines = &Simulation::ReadWindTurbinesHDF5;
        pReadOWCs = &Simulation::ReadOWCsHDF5;
    }
    else
    {
        std::stringstream ss;
        ss << "Simulation data format --> " << incDataFormat << " is not available.\n    Available formats: ASCII | HDF5.";
        throw ValueError(ss.str());
    }
}

void Simulation::UpdateSystem()
{
    // Update the bodies velocity buffers for the radiation forces computation
    if (numBodies > 0)
    {
        // Update time vector if any
        bool restoreMatrix = false;
        timeBufferCount++;
        if (timeBufferCount < timeBufferSize)
        {
            timeBuffer(0, timeBufferCount) = pTimeSolver->t;
        }
        else
        {
            arma::mat timeBufferNew = arma::zeros(1, timeBufferSize);
            timeBufferNew.cols(0, pBodies[0]->pHydro->GetNumPointsIrf() - 1) = timeBuffer.cols(timeBufferSize - pBodies[0]->pHydro->GetNumPointsIrf(), timeBufferSize - 1);
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

        // Look for the nodes that are not joints, and store the indexes
        for (int k = 0; k < 2; k++)
        {
            if (pLines[ii]->pLineBcps[k]->GetType() != 3)
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
    arma::mat output_position = arma::reshape(arma::strans(positionLinesCouplingVector), 3 * positionLinesCouplingVector.n_rows, 1);

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
    for (int i = 0; i < numLines; i++)
    {
        // Get the initial positions of the current line
        arma::mat initial_positions = full_matrix.rows(pLines[i]->ind4CouplingMat);
        pLines[i]->pos = initial_positions;
        // Set the line nodes velocity to zero
        pLines[i]->vel = arma::zeros(pLines[i]->N, 3);
        // Compute the forces
        pLines[i]->SEM_computeF();
        // Store the forces in the global vector
        forcesLinesCouplingVector.rows(pLines[i]->ind4CouplingMat) += pLines[i]->F;
    }

    // Remove the forces of the nodes that are not joints
    forcesLinesCouplingVector.shed_rows(indexesFairAnchor);
    // Convert the matrix to a column vector
    arma::mat output_force = arma::reshape(arma::strans(forcesLinesCouplingVector), 3 * forcesLinesCouplingVector.n_rows, 1);

    return output_force;
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
    double tol = timeIntAbsTol;
    double tolrelativa = timeIntRelTol;

    // Initial stop criterion variables
    double cantidadRel = 2 * tolrelativa;
    double cantidadAbs = 2 * tol;
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
    arma::mat fx = ComputeLinesForces(x);
    double normafxInicial = arma::norm(fx, 2);

    // Newton-Raphson loop
    while (iter < maxIter && ((cantidadAbs > tol) || (cantidadRel > tolrelativa)))
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
        while (arma::norm(new_step_vec, 2) > ((1 - sigma * rho) * arma::norm(fx)) && rho > 0.01)
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
        cantidadAbs = arma::norm(dk, "inf");
        cantidadRel = arma::norm(dk / x, "inf");
        iter = iter + 1;

        // Update the solution for next iteration
        x = xsol;
        fx = fxsol;
        // TODO: Check that forces are not computed twice for the same point
    }

    // Check if the maximum number of iterations was exceeded
    if (iter >= maxIter)
    {
        std::stringstream ss;
        ss << "The maximum number of iterations was exceeded. Convergence was not achieved \n";
        throw ValueError(ss.str());
    }
}
