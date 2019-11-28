
#include <armadillo>
#include "CummingsIrregular.hpp"


arma::mat CumminsIrregular::ComputeRadiationForces()
{	
	// Allocate radiation force solution vector
	arma::mat radiation_force = arma::zeros(6, 1);

	// Declare local auxiliary variables
	arma::mat time_local;
	arma::mat time_irf;
	arma::mat irf_local;
	arma::mat irf_local_interp;
	arma::mat vel_local;
	arma::mat vel_local_interp;
	arma::mat vel_local_filter;
	int idx_begin, idx_end;
	double dt = IRFTime(0, 1) - IRFTime(0, 0);
	int aux;
	// Compute radiation forces for the 6DOFs
	for(int ib=0; ib<numBodies; ib++)
	{
		for(int i=0; i<6; i++)
		{
			for(int j=0; j<6; j++)
			{
				if (pSim->timeBuffer(0, pSim->timeBufferCount) > IRFTime(0, (*pIRFPoints[ib])(i, j)))
				{
					//std::cout << "HydroDatabase::ComputeRadiationForces - TimeBuffer exceed" << std::endl;
					//std::cout << "HydroDatabase::ComputeRadiationForces - timeBuffer size: " << arma::size(pSim->timeBuffer) << std::endl;
					//std::cout << "HydroDatabase::ComputeRadiationForces - IRFTime: " << IRFTime(0, (*pIRFPoints[ib])(i, j)) << std::endl;
					
					// Get IRF function from the storage
					irf_local = (*pIRF[ib]).subcube(0, i, j, numPointsIRF-1, i, j);
					
					// Calculate begin and end indexes for matrix slicing
					for (int k=0; k<pSim->timeBufferCount; k++)
					{
						//std::cout << "here" << std::endl;
						if (pSim->timeBuffer(0, k) > (pSim->timeBuffer(0, pSim->timeBufferCount)-IRFTime(0, (*pIRFPoints[ib])(i, j))))
						{
							//std::cout << "here" << std::endl;
							(*pTimeStartPos)(ib, i, j) = k-1;
							//std::cout << "here" << std::endl;
							break;
						}
					}
					idx_begin = (*pTimeStartPos)(ib, i, j);
					idx_end = pSim->timeBufferCount;
					//std::cout << "HydroDatabase::ComputeRadiationForces - idx_begin: " << idx_begin << std::endl;
					//std::cout << "HydroDatabase::ComputeRadiationForces - idx_end: " << idx_end << std::endl;

					// Get velocity chunck from the storage
					vel_local = (*pBodies[ib]).velBuffer.submat(j, idx_begin, j, idx_end);

					// Interpolate local velocity vector in order to fit with IRF resolution
					//pSim->timeBuffer.cols(idx_begin, idx_end).print();
					//IRFTime.cols(0, (*pIRFPoints[ib])(i, j)).print();
					time_local = pSim->timeBuffer.cols(idx_begin, idx_end)-pSim->timeBuffer(0, idx_begin);
					time_irf = IRFTime.cols(0, (*pIRFPoints[ib])(i, j));
					vel_local_interp = interp1(time_local, vel_local, time_irf);

					// Calulate Duhamel integral
					vel_local_filter = arma::flipud(irf_local)%(vel_local_interp.t());
					vel_local_filter = vel_local_filter.t();
					radiation_force(i) += trapz(vel_local_filter, dt);
					/**
					if ((i==2) && (j==2))
					{
						std::cout << "dt: " << dt << std::endl;
						std::cout << "Radidation force (2,2): " << trapz(vel_local_filter, dt) << std::endl;
						std::cout << "Radidation force Unit time (2,2): " << trapz(vel_local_filter, 1.0) << std::endl;
						std::cout << "Radiation force: " << radiation_force(i) << std::endl;
						time_local.save(arma::hdf5_name("time_local_x.h5","irf"));
						time_irf.save(arma::hdf5_name("time_irf.h5","time_irf"));
						irf_local.save(arma::hdf5_name("irf_x.h5","irf"));
						vel_local.save(arma::hdf5_name("vel_x.h5","vel"));
						vel_local_interp.save(arma::hdf5_name("vel_local_interp_x.h5", "vel_local_interp"));
						vel_local_filter.save(arma::hdf5_name("vel_local_filter_x.h5","vel"));
						//time_local = (*SD.timeBuffer.cols(0, idx_end);
						//time_local.save(arma::hdf5_name("time_x.h5","time"));
					}
					**/
				}
				else if (pSim->timeBufferCount > 0)
				{
					// Calculate begin and end indexes for matrix slicing
					idx_end = pSim->timeBufferCount;

					// Get IRF function from the storage
					irf_local = (*pIRF[ib]).subcube(0, i, j, numPointsIRF-1, i, j);

					// Get velocity chunk from the storage
					vel_local = (*pBodies[ib]).velBuffer.submat(j, 0, j, idx_end);

					// Interpolate local velocity vector in order to fit with IRF resolution
					time_local = IRFTime.cols(0, numPointsIRF-1) - IRFTime(0, numPointsIRF-1) + pSim->timeBuffer(0, idx_end);
					irf_local_interp = interp1(time_local, arma::flipud(irf_local), pSim->timeBuffer.cols(0, idx_end));

					// Calulate Duhamel integral
					vel_local_filter = irf_local_interp%vel_local;
					radiation_force(i) += trapzi(pSim->timeBuffer.cols(0, idx_end), vel_local_filter);

					/**
					if ((idx_end == 500) && (i==2) && (j==2))
					{
						time_local.save(arma::hdf5_name("time_local_2500.h5","irf"));
						irf_local.save(arma::hdf5_name("irf_2500.h5","irf"));
						irf_local_interp.save(arma::hdf5_name("irf_interp_2500.h5","irf_interp"));
						vel_local.save(arma::hdf5_name("vel_2500.h5","vel"));
						time_local = (*SD.timeBuffer).cols(0, idx_end);
						time_local.save(arma::hdf5_name("time_2500.h5","time"));
					}
					**/
				}
			}
		}
	}
	
	return radiation_force;
}


arma::mat CumminsIrregular::ComputeFirstWaveExcForce()
{
	// Create local variables
	double time = pSim->pTimeSolver->t;
	// Get First Order Wave exciting data from storage
	arma::mat wave_exc_mag = pWaveExcitingMag->subcube(0, numPeriodExc, numHeadingExc, 5, numPeriodExc, numHeadingExc);
	arma::mat wave_exc_pha = pWaveExcitingPha->subcube(0, numPeriodExc, numHeadingExc, 5, numPeriodExc, numHeadingExc);
	// Calculate wave force
	double time_slope = 1/(*pFrequencies)(0, numPeriodExc);
	arma::mat wave_force = arma::zeros(6, 1);
	double angular_freq = 2*M_PI*(*pFrequencies)(0, numPeriodExc);
	for (int i=0; i<6; i++)
	{
		wave_force(i, 0) = wave_exc_mag(i, 0)*waveAmplitude*cos(angular_freq*time + wave_exc_pha(i, 0));

		if (time < time_slope)
		{
			wave_force(i, 0) = wave_force(i, 0)*time/time_slope;
		}
	}
	//std::cout << "Wave Amplitude: " << waveAmplitude << std::endl;

	return wave_force;
}
