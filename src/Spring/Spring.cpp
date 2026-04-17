
#include <armadillo>
#include <string>
#include <cstdio>
#include "Spring.hpp"
#include "../MathTools.hpp"
#include "../Exceptions/Exception.hpp"

// Lee inputs de los muelles
void Spring::ReadPropertiesASCII(std::string file_path)
{

	int ii, jj, kk, ll, temp_N;
	std::string Dummy;
	const int nInored = 97; // numero de lineas que se leen para cada nueva linea
	arma::mat temp_vec;

	SpringVectors.set_size(3, 2);
	data_StressStrain.set_size(6, 2);

	// Abro el fichero
	std::ifstream datosSprings(file_path);

	// Ignoro la primera linea del fichero, que contiene el numero de muelles a estudiar
	datosSprings >> Dummy;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n'); // El ignore sirve para ignorar el texto de la linea

	// Ignoro las lineas que ya se han leido
	for (ii = 0; ii < nSpring; ii = ii + 1)
	{
		for (jj = 1; jj <= nInored; jj = jj + 1)
		{
			datosSprings >> Dummy;
			datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		}
	}

	// Ignoro las tres primeras lineas, donde pone "New spring"
	for (ii = 1; ii <= 3; ii = ii + 1)
	{
		datosSprings >> Dummy;
		datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	}

	// Leo todo
	datosSprings >> stressModelFlag;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> dampingFlag;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> frictionFlag;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> frameFlag;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> BCP_1;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> BCP_2;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> BCP_1_type;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> BCP_2_type;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	BCP_1 -= 1;
	BCP_2 -= 1;

	if (BCP_1_type > 0 && BCP_2_type > 0)
	{
		std::stringstream ss;
		ss << "One of the BCPs must be a fixed point"
		   << ".\n";
		throw ValueError(ss.str());
	}

	for (jj = 0; jj < 2; jj = jj + 1)
	{
		for (ii = 0; ii < 3; ii = ii + 1)
		{
			temp_vec = arma::zeros(3, 1);
			datosSprings >> temp_vec(0, 0);
			datosSprings >> temp_vec(1, 0);
			datosSprings >> temp_vec(2, 0);
			datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
			SpringVectors(ii, jj) = temp_vec;
		}
	}

	// Check and post process the spring base
	std::cout << "        Post-processing Spring " << nSpring + 1 << " vectors..." << std::endl;
	arma::field<arma::mat> SpringVectors_new;
	SpringVectors_new.set_size(3, 2);
	double temp_norm;
	for (jj = 0; jj < 2; jj = jj + 1)
	{
		for (ii = 0; ii < 3; ii = ii + 1)
		{
			temp_vec = SpringVectors(ii, jj);
			temp_norm = arma::norm(temp_vec);
			if (temp_norm < 0.99 || temp_norm > 1.01)
			{
				std::stringstream ss;
				ss << "The norm of the spring vectors must be unitary"
				   << ".\n";
				throw ValueError(ss.str());
			}
			temp_vec = temp_vec / temp_norm;
			SpringVectors_new(ii, jj) = temp_vec;
		}
		SpringVectors_new(1, jj) = arma::cross(SpringVectors_new(2, jj), SpringVectors_new(0, jj));
		SpringVectors_new(1, jj) = SpringVectors_new(1, jj) / arma::norm(SpringVectors_new(1, jj));
		SpringVectors_new(2, jj) = arma::cross(SpringVectors_new(0, jj), SpringVectors_new(1, jj));
		for (ii = 0; ii < 3; ii = ii + 1)
		{
			temp_norm = arma::norm(SpringVectors_new(ii, jj) - SpringVectors(ii, jj));
			if (temp_norm > 0.02)
			{
				std::stringstream ss;
				ss << "The spring vectors must form an orthonormal basis"
				   << ".\n";
				throw ValueError(ss.str());
			}
		}
	}
	for (jj = 0; jj < 2; jj = jj + 1)
	{
		std::cout << "          Spring " << jj + 1 << " ..." << std::endl;
		for (ii = 0; ii < 3; ii = ii + 1)
		{
			temp_norm = arma::norm(SpringVectors_new(ii, jj) - SpringVectors(ii, jj));
			if (temp_norm > 1e-6)
			{
				std::cout << "            Vector " << ii + 1 << " changed from: " << std::endl
						  << "            " << SpringVectors(ii, jj).t() << "            to" << std::endl
						  << "            " << SpringVectors_new(ii, jj).t();
			}
		}
		std::cout << "          ... checked! " << std::endl;
	}
	SpringVectors = SpringVectors_new;

	datosSprings >> Dummy;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	for (jj = 0; jj < 6; jj = jj + 1)
	{
		for (kk = 0; kk < 6; kk = kk + 1)
		{
			datosSprings >> SpringMatrix_K(jj, kk);
		}
		datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	}

	datosSprings >> Dummy;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> mu_d;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> mu_s;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> vt;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	datosSprings >> Dt;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');

	arma::mat tmpA, tmpB = arma::zeros(4, 1);

	tmpA = {{pow(vt / 2, 3), pow(vt / 2, 2), vt / 2, 1},
			{3 * pow(vt / 2, 2), 2 * (vt / 2), 1, 0},
			{pow(vt, 3), pow(vt, 2), vt, 1},
			{3 * pow(vt, 2), 2 * vt, 1, 0}};
	tmpB(0, 0) = mu_s / 2;
	tmpB(1, 0) = mu_s / vt;
	tmpB(2, 0) = mu_s;
	tmpB(3, 0) = 0;
	a_1 = arma::solve(tmpA, tmpB);

	tmpA = {{pow(vt, 3), pow(vt, 2), vt, 1},
			{3 * pow(vt, 2), 2 * vt, 1, 0},
			{pow(vt * 2, 3), pow(vt * 2, 2), vt * 2, 1},
			{3 * pow(vt * 2, 2), 2 * (vt * 2), 1, 0}};
	tmpB(0, 0) = mu_s;
	tmpB(1, 0) = 0;
	tmpB(2, 0) = mu_d;
	tmpB(3, 0) = 0;
	a_2 = arma::solve(tmpA, tmpB);

	datosSprings >> Dummy;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	for (jj = 0; jj < 6; jj = jj + 1)
	{
		for (kk = 0; kk < 6; kk = kk + 1)
		{
			datosSprings >> SpringMatrix_M(jj, kk);
		}
		datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	}

	datosSprings >> Dummy;
	datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	for (jj = 0; jj < 6; jj = jj + 1)
	{
		datosSprings >> SpringMatrix_D(jj, 0);
		datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
	}

	for (jj = 0; jj < 6; jj = jj + 1)
	{
		datosSprings >> Dummy;
		datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		datosSprings >> temp_N;
		datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		n_StressStrain[jj] = temp_N;
		arma::mat temp_vec2 = arma::zeros(temp_N, 1);
		arma::mat temp_mat = arma::zeros(temp_N, 6);
		for (ll = 0; ll < temp_N; ll = ll + 1)
		{
			datosSprings >> temp_vec2(ll, 0);
		}
		datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		data_StressStrain(jj, 0) = temp_vec2;
		temp_mat = arma::zeros(temp_N, 6);
		for (kk = 0; kk < 6; kk = kk + 1)
		{
			for (ll = 0; ll < temp_N; ll = ll + 1)
			{
				datosSprings >> temp_mat(ll, kk);
			}
			datosSprings.ignore(std::numeric_limits<int>::max(), '\n');
		}
		data_StressStrain(jj, 1) = temp_mat;
		// std::cout << "temp_vec2  " << temp_vec2 << std::endl;
		// std::cout << "temp_mat  " << temp_mat << std::endl;
	}

	// Cierro el fichero
	datosSprings.close();
}

void Spring::ReadPropertiesYAML(YAML::Node node)
{
	arma::mat temp_vec;

	SpringVectors.set_size(3, 2);
	data_StressStrain.set_size(6, 2);

	stressModelFlag = node["stress_model_flag"].as<int>();
	dampingFlag = node["damping_flag"].as<int>();
	frictionFlag = node["friction_flag"].as<int>();
	frameFlag = node["frame_flag"].as<int>();
	BCP_1 = node["BCP_1"].as<int>() - 1;
	BCP_2 = node["BCP_2"].as<int>() - 1;
	BCP_1_type = node["BCP_1_type"].as<int>();
	BCP_2_type = node["BCP_2_type"].as<int>();

	if (BCP_1_type > 0 && BCP_2_type > 0)
	{
		std::stringstream ss;
		ss << "One of the BCPs must be a fixed point.\n";
		throw ValueError(ss.str());
	}

	YAML::Node vecNode = node["vectors"];
	for (int jj = 0; jj < 2; jj++)
	{
		for (int ii = 0; ii < 3; ii++)
		{
			temp_vec = arma::zeros(3, 1);
			YAML::Node v = vecNode[jj * 3 + ii];
			temp_vec(0, 0) = v[0].as<double>();
			temp_vec(1, 0) = v[1].as<double>();
			temp_vec(2, 0) = v[2].as<double>();
			SpringVectors(ii, jj) = temp_vec;
		}
	}

	// Check and post process the spring base
	std::cout << "        Post-processing Spring " << nSpring + 1 << " vectors..." << std::endl;
	arma::field<arma::mat> SpringVectors_new;
	SpringVectors_new.set_size(3, 2);
	double temp_norm;
	for (int jj = 0; jj < 2; jj++)
	{
		for (int ii = 0; ii < 3; ii++)
		{
			temp_vec = SpringVectors(ii, jj);
			temp_norm = arma::norm(temp_vec);
			if (temp_norm < 0.99 || temp_norm > 1.01)
			{
				std::stringstream ss;
				ss << "The norm of the spring vectors must be unitary.\n";
				throw ValueError(ss.str());
			}
			temp_vec = temp_vec / temp_norm;
			SpringVectors_new(ii, jj) = temp_vec;
		}
		SpringVectors_new(1, jj) = arma::cross(SpringVectors_new(2, jj), SpringVectors_new(0, jj));
		SpringVectors_new(1, jj) = SpringVectors_new(1, jj) / arma::norm(SpringVectors_new(1, jj));
		SpringVectors_new(2, jj) = arma::cross(SpringVectors_new(0, jj), SpringVectors_new(1, jj));
		for (int ii = 0; ii < 3; ii++)
		{
			temp_norm = arma::norm(SpringVectors_new(ii, jj) - SpringVectors(ii, jj));
			if (temp_norm > 0.02)
			{
				std::stringstream ss;
				ss << "The spring vectors must form an orthonormal basis.\n";
				throw ValueError(ss.str());
			}
		}
	}
	for (int jj = 0; jj < 2; jj++)
	{
		std::cout << "          Spring " << jj + 1 << " ..." << std::endl;
		for (int ii = 0; ii < 3; ii++)
		{
			temp_norm = arma::norm(SpringVectors_new(ii, jj) - SpringVectors(ii, jj));
			if (temp_norm > 1e-6)
			{
				std::cout << "            Vector " << ii + 1 << " changed from: " << std::endl
						  << "            " << SpringVectors(ii, jj).t() << "            to" << std::endl
						  << "            " << SpringVectors_new(ii, jj).t();
			}
		}
		std::cout << "          ... checked! " << std::endl;
	}
	SpringVectors = SpringVectors_new;

	// Read stiffness matrix
	YAML::Node kNode = node["stiffness_matrix"];
	for (int jj = 0; jj < 6; jj++)
		for (int kk = 0; kk < 6; kk++)
			SpringMatrix_K(jj, kk) = kNode[jj][kk].as<double>();

	// Read friction parameters
	mu_d = node["mu_d"].as<double>();
	mu_s = node["mu_s"].as<double>();
	vt = node["vt"].as<double>();
	Dt = node["Dt"].as<double>();

	arma::mat tmpA, tmpB = arma::zeros(4, 1);
	tmpA = {{pow(vt / 2, 3), pow(vt / 2, 2), vt / 2, 1},
			{3 * pow(vt / 2, 2), 2 * (vt / 2), 1, 0},
			{pow(vt, 3), pow(vt, 2), vt, 1},
			{3 * pow(vt, 2), 2 * vt, 1, 0}};
	tmpB(0, 0) = mu_s / 2;
	tmpB(1, 0) = mu_s / vt;
	tmpB(2, 0) = mu_s;
	tmpB(3, 0) = 0;
	a_1 = arma::solve(tmpA, tmpB);

	tmpA = {{pow(vt, 3), pow(vt, 2), vt, 1},
			{3 * pow(vt, 2), 2 * vt, 1, 0},
			{pow(vt * 2, 3), pow(vt * 2, 2), vt * 2, 1},
			{3 * pow(vt * 2, 2), 2 * (vt * 2), 1, 0}};
	tmpB(0, 0) = mu_s;
	tmpB(1, 0) = 0;
	tmpB(2, 0) = mu_d;
	tmpB(3, 0) = 0;
	a_2 = arma::solve(tmpA, tmpB);

	// Read mass matrix
	YAML::Node mNode = node["mass_matrix"];
	for (int jj = 0; jj < 6; jj++)
		for (int kk = 0; kk < 6; kk++)
			SpringMatrix_M(jj, kk) = mNode[jj][kk].as<double>();

	// Read damping vector
	YAML::Node dNode = node["damping_vector"];
	for (int jj = 0; jj < 6; jj++)
		SpringMatrix_D(jj, 0) = dNode[jj].as<double>();

	// Read stress-strain data for each DOF
	YAML::Node ssNode = node["stress_strain"];
	for (int jj = 0; jj < 6; jj++)
	{
		YAML::Node dofNode = ssNode[jj];
		int temp_N = (int)dofNode["displacements"].size();
		n_StressStrain[jj] = temp_N;
		arma::mat temp_vec2 = arma::zeros(temp_N, 1);
		for (int ll = 0; ll < temp_N; ll++)
			temp_vec2(ll, 0) = dofNode["displacements"][ll].as<double>();
		data_StressStrain(jj, 0) = temp_vec2;

		arma::mat temp_mat = arma::zeros(temp_N, 6);
		YAML::Node forcesNode = dofNode["forces"];
		for (int kk = 0; kk < 6; kk++)
			for (int ll = 0; ll < temp_N; ll++)
				temp_mat(ll, kk) = forcesNode[kk][ll].as<double>();
		data_StressStrain(jj, 1) = temp_mat;
	}
}

// Calcula las fuerzas que aplica el muelle en los BCPs y las guarda en estos
void Spring::computeSpringForces(void)
{
	// Calculo los vectores unitarios del muelle en global para cada cuerpo
	arma::field<arma::mat> SpringVectorsG;
	arma::field<arma::mat> SpringVectorsG_dot;
	SpringVectorsG.set_size(3, 2);
	SpringVectorsG_dot.set_size(3, 2);
	for (int jj = 0; jj < 2; jj = jj + 1)
	{ // Bucle sobre los dos BCPs
		for (int ii = 0; ii < 3; ii = ii + 1)
		{ // Bucle sobre los tres vectores del muelle
			SpringVectorsG(ii, jj) = SpringBCP[jj]->rotMat * SpringVectors(ii, jj);
			SpringVectorsG_dot(ii, jj) = SpringBCP[jj]->rotMat_dot * SpringVectors(ii, jj);
		}
	}

	// Actualizo las posiciones del BCP que no sea un punto fijo, si se da el caso.
	arma::mat temp_p, temp_p0, temp_v;
	if (BCP_1_type > 0 && flagStickSlip == 0)
	{
		temp_p0 = SpringBCP[0]->posG_BCP; // Point of BCP1 previous position
		temp_p = SpringBCP[1]->posG_BCP;  // BCP2 position to project over BCP1
		if (BCP_1_type == 1)
		{
			temp_v = SpringVectorsG(2, 0);																  // BCP1 z vector taken as line vector
			SpringBCP[0]->posG_BCP = temp_p0 + arma::as_scalar(temp_v.t() * (temp_p - temp_p0)) * temp_v; // Update BCP1 position with projection on line
		}
		else if (BCP_1_type == 2)
		{
			temp_v = SpringVectorsG(0, 0);																 // BCP1 x vector taken as plane normal vector
			SpringBCP[0]->posG_BCP = temp_p - arma::as_scalar(temp_v.t() * (temp_p - temp_p0)) * temp_v; // Update BCP1 position with projection on plane
		}
		SpringBCP[0]->posWrtCdgGlobal = SpringBCP[0]->posWrtCdgGlobal + SpringBCP[0]->posG_BCP - temp_p0;
		SpringBCP[0]->posWrtCdgLocal = SpringBCP[0]->rotMat.t() * SpringBCP[0]->posWrtCdgGlobal;
	}
	if (BCP_2_type > 0 && flagStickSlip == 0)
	{
		temp_p0 = SpringBCP[1]->posG_BCP; // Point of BCP2 previous position
		temp_p = SpringBCP[0]->posG_BCP;  // BCP1 position to project over BCP1
		if (BCP_2_type == 1)
		{
			temp_v = SpringVectorsG(2, 1);																  // BCP2 z vector taken as line vector
			SpringBCP[1]->posG_BCP = temp_p0 + arma::as_scalar(temp_v.t() * (temp_p - temp_p0)) * temp_v; // Update BCP2 position with projection on line
		}
		else if (BCP_2_type == 2)
		{
			temp_v = SpringVectorsG(0, 1);																 // BCP2 x vector taken as plane normal vector
			SpringBCP[1]->posG_BCP = temp_p - arma::as_scalar(temp_v.t() * (temp_p - temp_p0)) * temp_v; // Update BCP2 position with projection on plane
		}
		SpringBCP[1]->posWrtCdgGlobal = SpringBCP[1]->posWrtCdgGlobal + SpringBCP[1]->posG_BCP - temp_p0;
		SpringBCP[1]->posWrtCdgLocal = SpringBCP[1]->rotMat.t() * SpringBCP[1]->posWrtCdgGlobal;
	}

	// Calculo el vector en global que une los BCPs
	arma::mat L12 = SpringBCP[1]->posG_BCP - SpringBCP[0]->posG_BCP;
	arma::mat L12_dot = SpringBCP[1]->velG_BCP - SpringBCP[0]->velG_BCP;

	arma::mat x1 = SpringVectorsG(0, 0);
	arma::mat y1 = SpringVectorsG(1, 0);
	arma::mat z1 = SpringVectorsG(2, 0);
	arma::mat x2 = SpringVectorsG(0, 1);
	arma::mat y2 = SpringVectorsG(1, 1);
	arma::mat z2 = SpringVectorsG(2, 1);
	arma::mat x1_dot = SpringVectorsG_dot(0, 0);
	arma::mat y1_dot = SpringVectorsG_dot(1, 0);
	arma::mat z1_dot = SpringVectorsG_dot(2, 0);
	arma::mat x2_dot = SpringVectorsG_dot(0, 1);
	arma::mat y2_dot = SpringVectorsG_dot(1, 1);
	arma::mat z2_dot = SpringVectorsG_dot(2, 1);

	arma::mat mat_spring2global, mat_spring2global_dot;
	if (frameFlag == 0)
	{
		mat_spring2global = 0.5 * (arma::join_horiz(x1, y1, z1) + arma::join_horiz(x2, y2, z2));
		mat_spring2global_dot = 0.5 * (arma::join_horiz(x1_dot, y1_dot, z1_dot) + arma::join_horiz(x2_dot, y2_dot, z2_dot));
	}
	else if (frameFlag == 1)
	{
		mat_spring2global = arma::join_horiz(x1, y1, z1);
		mat_spring2global_dot = arma::join_horiz(x1_dot, y1_dot, z1_dot);
	}
	else if (frameFlag == 2)
	{
		mat_spring2global = arma::join_horiz(x2, y2, z2);
		mat_spring2global_dot = arma::join_horiz(x2_dot, y2_dot, z2_dot);
	}
	else
	{
		std::stringstream ss;
		ss << "Spring frame option not available"
		   << ".\n";
		throw ValueError(ss.str());
	}
	arma::mat mat_global2spring = mat_spring2global.t();
	arma::mat mat_global2spring_dot = mat_spring2global_dot.t();

	arma::mat SpringStrains = arma::zeros(6, 1);	 // Deformacion en el muelle
	arma::mat SpringStrains_dot = arma::zeros(6, 1); // Derivada temporal de la deformacion en el muelle

	SpringStrains.rows(0, 2) = mat_global2spring * L12;
	SpringStrains(3, 0) = 0.5 * (atan2(arma::dot(y2, z1), arma::dot(z2, z1)) - atan2(arma::dot(y1, z2), arma::dot(z2, z1)));
	SpringStrains(4, 0) = 0.5 * (asin(arma::dot(-x2, z1)) - asin(arma::dot(-x1, z2)));
	SpringStrains(5, 0) = 0.5 * (atan2(arma::dot(x2, y1), arma::dot(x2, x1)) - atan2(arma::dot(x1, y2), arma::dot(x2, x1)));

	// Lo mismo para la derivada temporal
	if (dampingFlag > 0 || frictionFlag > 0)
	{
		SpringStrains_dot.rows(0, 2) = mat_global2spring_dot * L12 + mat_global2spring * L12_dot;
		SpringStrains_dot(3, 0) = 0.5 * (((arma::dot(y2_dot, z1) + arma::dot(y2, z1_dot)) * arma::dot(z2, z1) + arma::dot(y2, z1) * (arma::dot(z2_dot, z1) + arma::dot(z2, z1_dot))) / (pow(arma::dot(z2, z1), 2) + pow(arma::dot(y2, z1), 2)) - ((arma::dot(y1_dot, z2) + arma::dot(y1, z2_dot)) * arma::dot(z1, z2) + arma::dot(y1, z2) * (arma::dot(z1_dot, z2) + arma::dot(z1, z2_dot))) / (pow(arma::dot(z1, z2), 2) + pow(arma::dot(y1, z2), 2)));
		SpringStrains_dot(4, 0) = 0.5 * ((arma::dot(x2_dot, z1) + arma::dot(x2, z1_dot)) / sqrt(1.0 - pow(arma::dot(x2, z1), 2)) - (arma::dot(x1_dot, z2) + arma::dot(x1, z2_dot)) / sqrt(1.0 - pow(arma::dot(x1, z2), 2)));
		SpringStrains_dot(5, 0) = 0.5 * (((arma::dot(x2_dot, y1) + arma::dot(x2, y1_dot)) * arma::dot(x2, x1) + arma::dot(x2, y1) * (arma::dot(x2_dot, x1) + arma::dot(x2, x1_dot))) / (pow(arma::dot(x2, y1), 2) + pow(arma::dot(x2, x1), 2)) - ((arma::dot(x1_dot, y2) + arma::dot(x1, y2_dot)) * arma::dot(x1, x2) + arma::dot(x1, y2) * (arma::dot(x1_dot, x2) + arma::dot(x1, x2_dot))) / (pow(arma::dot(x1, y2), 2) + pow(arma::dot(x1, x2), 2)));
	}

	arma::mat tempF_L = arma::zeros(6, 1); // Fuerza registrada en el muelle en coordenadas del cuerpo
	arma::mat tempF_G = arma::zeros(6, 1); // Fuerza registrada en el muelle en coordenadas globales

	// Variables temporales necesarias mas adelante
	arma::mat temp_strainData, temp_stressData, tempS;
	arma::mat tempF_K = arma::zeros(6, 1);

	if (stressModelFlag == 1)
	{ // Modelo de muelle lineal y simetrico
		// Fuerza registrada en el muelle en coordenadas del cuerpo, multiplicando matriz y strain
		tempF_L = SpringMatrix_K * SpringStrains;
		tempF_K = tempF_L;
	}
	else if (stressModelFlag == 2)
	{ // Modelo de muelle no lineal, input es curva de deformacion-fuerza.
		// Fuerza registrada en el muelle en coordenadas del cuerpo, interpolando las curvas de esfuerzo-deformacion
		for (int jj = 0; jj < 6; jj = jj + 1)
		{

			tempS = SpringStrains(jj, 0);
			temp_strainData = data_StressStrain(jj, 0);
			temp_stressData = data_StressStrain(jj, 1);

			tempF_L = tempF_L + (interp1(temp_strainData, temp_stressData, tempS)).t();
		}
		tempF_K = tempF_L;
	}

	// Fuerza de damping
	arma::mat tempF_D = arma::zeros(6, 1);
	if (dampingFlag == 1)
	{
		tempF_D = (SpringMatrix_D % arma::abs(arma::sign(tempF_K))) % SpringStrains_dot;
		tempF_L = tempF_L + tempF_D;
	}

	// Fuerza de friccion
	arma::mat tempF_F = arma::zeros(6, 1);
	if (frictionFlag == 1)
	{

		arma::mat SpringStrains_dot_abs = arma::abs(SpringStrains_dot);

		arma::mat mu_vec = arma::zeros(6, 1);
		for (int jj = 0; jj < 6; jj = jj + 1)
		{
			double vv = arma::as_scalar(SpringStrains_dot_abs(jj, 0));
			if (vv < vt / 2)
			{
				mu_vec(jj, 0) = mu_s / vt * vv;
			}
			else if (vv >= vt / 2 && vv < vt)
			{
				mu_vec(jj, 0) = arma::as_scalar(a_1(0, 0) * pow(vv, 3) + a_1(1, 0) * pow(vv, 2) + a_1(2, 0) * vv + a_1(3, 0));
			}
			else if (vv >= vt && vv < vt * 2)
			{
				mu_vec(jj, 0) = arma::as_scalar(a_2(0, 0) * pow(vv, 3) + a_2(1, 0) * pow(vv, 2) + a_2(2, 0) * vv + a_2(3, 0));
			}
			else
			{
				mu_vec(jj, 0) = mu_d;
			}
		}
		mu_vec = mu_vec % arma::sign(SpringStrains_dot);

		tempF_F = (SpringMatrix_M * arma::abs(tempF_K)) % mu_vec;
	}
	else if (frictionFlag == 2)
	{

		arma::mat fn = SpringMatrix_M * arma::abs(tempF_K);
		double fn_norm = fn.max();
		arma::mat VelSlip = arma::abs(arma::sign(fn)) % SpringStrains_dot;
		double v_norm = arma::norm(VelSlip);
		if (flagStickSlip == 0 && fn_norm > 1e-1 && v_norm < vt)
		{
			flagStickSlip = 1;
			SpringStrains_Stick = SpringStrains;
		}
		if (fn_norm <= 1e-1 || v_norm >= vt)
		{
			flagStickSlip = 0;
		}

		double bv = step(v_norm, -vt, -1.0, vt, 1.0);

		if (v_norm > 1e-5)
		{
			tempF_F = tempF_F + (VelSlip % fn) * bv * mu_d / v_norm;
		}

		if (flagStickSlip == 1)
		{
			arma::mat DeltaStick = arma::abs(arma::sign(fn)) % (SpringStrains - SpringStrains_Stick);
			double D = arma::norm(DeltaStick);
			double bD = step(D, -Dt, -1.0, Dt, 1.0);
			if (D > 1e-5)
			{
				tempF_F = tempF_F + (DeltaStick % fn) * (1 - bv) * bD * mu_s / D;
			}
		}

		// if (fn_norm>1e-1){
		// 	std::cout << "arma::sign(SpringStrains_dot) = " << std::endl << arma::sign(SpringStrains_dot);
		// 	std::cout << "fn = " << std::endl << fn;
		// 	std::cout << "v_norm = "  << v_norm << std::endl;
		// 	std::cout << "bv = "  << bv << std::endl;
		// 	std::cout << "step(0.1,-0.2,-1.0,0.2,1.0) = "  << step(0.1,-0.2,-1.0,0.2,1.0) << std::endl;
		// 	std::cout << "mu_d = " << mu_d << std::endl;
		// 	std::cout << "flagStickSlip = " << flagStickSlip << std::endl;
		// 	std::cout << "tempF_F = " << std::endl << tempF_F << std::endl;
		// }
	}

	tempF_L = tempF_L + tempF_F;

	// Paso de coordenadas locales a coordenadas globales la fuerza del muelle
	tempF_G = arma::zeros(6, 1);
	tempF_G.rows(0, 2) = mat_spring2global * (tempF_L.rows(0, 2));
	tempF_G.rows(3, 5) = mat_spring2global * (tempF_L.rows(3, 5));

	SpringBCP[0]->forceBcp = SpringBCP[0]->forceBcp + tempF_G; // Acumulo la fuerza obtenida en el BCP
	SpringBCP[1]->forceBcp = SpringBCP[1]->forceBcp - tempF_G;

	SpringBCP[0]->temp = SpringBCP[0]->temp + tempF_L;
	SpringBCP[1]->temp = SpringBCP[1]->temp - tempF_L;

	/*
	std::cout << "L12 = " << std::endl << L12;
	std::cout << "SpringStrains = " << std::endl << SpringStrains;
	std::cout << "L12_dot = " << std::endl << L12_dot;
	std::cout << "SpringStrains_dot = " << std::endl << SpringStrains_dot;
	std::cout << "tempF_K = " << std::endl << tempF_K;
	std::cout << "tempF_D = " << std::endl << tempF_D;
	std::cout << "tempF_F = " << std::endl << tempF_F;
	std::cout << "tempF_L = " << std::endl << tempF_L;
	std::cout << "tempF_G = " << std::endl << tempF_G << std::endl;
	std::stringstream ss;
	ss << "STOP" << ".\n";
	throw ValueError(ss.str());

	*/
}