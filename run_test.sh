module purge
ml Armadillo/10.5.3-foss-2020b
ml FASTurbine_wrapper/0.0-foss-2020b

cd cases/test_case
/home/projects/energia/developments/oasis_wind_turbine/oasis/bin/oasis ./
rm core*
cd ..
cd ..