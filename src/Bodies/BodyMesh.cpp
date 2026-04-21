// SPDX-License-Identifier: GPL-3.0-or-later
/*
Library for body meshes
*/

#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <armadillo>
#include "BodyMesh.hpp"
#include "../Bodies/Bodies.hpp"
#include "../Exceptions/Exception.hpp"
#include "../Simulations/Simulation.hpp"
#include "../Logger.hpp"
#include "../os_tools.hpp"
#include "../MathTools.hpp"
#ifdef OASIS_USE_STL_READER
#include <stl_reader.h>
#endif

////////////////////////////////////////////////////////////////////////////
/////////////////////////// BodyMesh CLASS DEFINITION //////////////////////
////////////////////////////////////////////////////////////////////////////

BodyMesh::BodyMesh(int incId, std::string incMeshFileName, Body* incpBody)
{
    id = incId;
    meshFileName = incMeshFileName;
    pBody = incpBody;

    pSim = pBody->pSim;
}

int BodyMesh::GetId(void)
{
    return this->id;
}

int BodyMesh::GetType(void)
{
    return this->typeMesh;
}

void BodyMesh::ReadPropertiesASCII(void)
{
    Logger::info("--> Reading Body 2D Mesh (ASCII format)");

    int numElemNodes = 3; // HARDCODED to read STL

#ifndef OASIS_USE_STL_READER
    throw std::runtime_error("BodyMesh::ReadPropertiesASCII requires stl_reader. "
                             "Rebuild with OASIS_USE_STL_READER=ON.");
#else
    // Read body mesh
    stl_reader::StlMesh<float, unsigned int> mesh(meshFileName);

    iniNumElems = mesh.num_tris();
    arma::mat tmp_nodes = arma::zeros(iniNumElems * numElemNodes, 3);
    arma::umat tmp_elems = arma::zeros<arma::umat>(iniNumElems, numElemNodes);

    for (int ielem = 0; ielem < iniNumElems; ++ielem)
    {
        for (int inode = 0; inode < numElemNodes; ++inode)
        {
            const float* c = mesh.vrt_coords(mesh.tri_corner_ind(ielem, inode));
            tmp_nodes.row(3 * ielem + inode) = {c[0], c[1], c[2]};
            tmp_elems(ielem, inode) = 3 * ielem + inode;
        }
    }

    // Duplicated vertex nodes removal
    arma::uvec ind;
    std::tie(iniNodes, ind) = unique_rows(tmp_nodes);
    iniNumNodes = iniNodes.n_rows;

    // If iniNodes contains values larger than 1000, convert from millimeters to meters and display a warning
    if (arma::abs(iniNodes).max() > 1000.0)
    {
        Logger::warning("It seems that the mesh is defined in millimeters. Converting to meters...");
        iniNodes = iniNodes / 1000.0;
    }

    // Rearrange elems index
    iniElems = indMat(ind, tmp_elems);

    Logger::debug("  --> num_nodes = " + std::to_string(iniNumNodes));
    Logger::debug("  --> num_elems = " + std::to_string(iniNumElems));
#endif // OASIS_USE_STL_READER
}

void BodyMesh::TransformMesh(void)
{
    // std::cout << "--> Transforming Body 2D Mesh" << std::endl;

    // Rotation matrix
    arma::mat rotMat = pBody->rotMat;
    arma::vec pos = pBody->pos.rows(0, 2);

    // Nodes transformation
    nodes = (pos * arma::ones(1, iniNumNodes) + rotMat * iniNodes.t()).t();

    // Normals transformation
    normals = (rotMat * iniNormals.t()).t();
}

void BodyMesh::CutMesh(double time)
{
    // std::cout << "--> Cutting Body 2D Mesh" << std::endl;

    double LAMBDA = 0.4; // RATIO FOR SOME GAUSSIAN NODES
    double TOL = 1.e-6;  // TOLERANCE FOR NUMERIC BOUND

    // Temporary variables
    arma::mat transNodes = nodes, transNormals = normals;

    // Sorting elements
    arma::uvec verticesUW;
    arma::vec eta;
    if (pBody->flag_hydrostatics == 1)
    {
        verticesUW = (transNodes.submat(0, 2, numVertices - 1, 2) <= 0.0);
    }
    else if (pBody->flag_hydrostatics == 2)
    {
        eta = pSim->pWave->GetFreeSurface(time, transNodes.col(0), transNodes.col(1));
        verticesUW = (transNodes.submat(0, 2, numVertices - 1, 2)) <= eta.subvec(0, numVertices - 1);
    }

    // Element status and UW index
    arma::uvec status = arma::sum(indMat(verticesUW, iniElems.cols(0, 2)), 1);
    arma::uvec ind1 = arma::find(status == 1);
    int numElem1 = ind1.n_rows;
    arma::uvec ind2 = arma::find(status == 2);
    int numElem2 = ind2.n_rows;
    arma::uvec ind3 = arma::find(status == 3);
    int numElem3 = ind3.n_rows;

    // Completely submerged elements and their nodes
    arma::uvec nodesUW, indMap;
    arma::umat auxNodes = iniElems.rows(ind3);
    if (pBody->flag_hydrostatics == 1)
    {
        nodesUW = arma::find(transNodes.col(2) <= 0.0);
    }
    else if (pBody->flag_hydrostatics == 2)
    {
        nodesUW = unique(arma::vectorise(auxNodes.t()));
    }
    int numNodesUW = nodesUW.n_rows;
    indMap = arma::zeros<arma::uvec>(iniNumNodes);
    indMap.rows(nodesUW) = arma::linspace<arma::uvec>(0, numNodesUW - 1, numNodesUW);

    // Including elements with 3 vertices submerged
    nodes = transNodes.rows(nodesUW);
    normals = transNormals.rows(ind3);
    elems = indMat(indMap, auxNodes);
    jacobians = iniJacobians(ind3);
    numNodes = numNodesUW;
    numElems = numElem3;

    // Loop variables initialization
    arma::rowvec p_1, p_2, p_3, p_4, p_12, p_23, p_31, p_34, p_41, p_123, p_134, p_1234;
    arma::rowvec v_1, v_2, v_3, v_4, w_1, w_2, tmpNormal;
    arma::rowvec pi_1, pi_2, pi_3, v_p, v_pi;
    arma::rowvec q_2, q_3, q_4;
    arma::urowvec auxVec, auxVec1, auxVec2, tmpElems;
    arma::uvec sortInd, icol;
    arma::vec jac, incEta;
    arma::mat incNodes, outNodes;
    double mu_1, mu_2, jac1, jac2, lambda;
    q_2 = q_3 = q_4 = arma::zeros<arma::rowvec>(3);

    // Including elements with 2 vertices submerged
    for (int ielem : ind2)
    {

        if (pBody->flag_hydrostatics == 1)
        {
            // Sorting vertices by height
            incNodes = sort_rows(transNodes.rows(iniElems(ielem, arma::span(0, 2))), 2);
        }
        else if (pBody->flag_hydrostatics == 2)
        {
            // Sorting vertices by height wrt free surface
            tmpElems = iniElems(ielem, arma::span(0, 2));
            icol = {2};
            sortInd = arma::sort_index(transNodes.submat(tmpElems, icol) - eta(tmpElems));
            incNodes = transNodes.rows(tmpElems(sortInd));
            incEta = eta.rows(tmpElems(sortInd));
        }

        // --------------------- Cutting element with two vertices underwater ---------------------

        // Initial vertex nodes
        p_1 = incNodes.row(0);
        p_2 = incNodes.row(1);
        p_3 = incNodes.row(2);

        // Triangle cut
        if (pBody->flag_hydrostatics == 1)
        {
            w_1 = p_3 - p_1;
            w_2 = p_3 - p_2;
            mu_1 = arma::as_scalar(-p_1.col(2) / w_1.col(2));
            mu_2 = arma::as_scalar(-p_2.col(2) / w_2.col(2));
            q_3.subvec(0, 1) = p_2.subvec(0, 1) + mu_2 * w_2.subvec(0, 1);
            q_3(2) = 0;
            q_4.subvec(0, 1) = p_1.subvec(0, 1) + mu_1 * w_1.subvec(0, 1);
            q_4(2) = 0;
        }
        else if (pBody->flag_hydrostatics == 2)
        {
            // Vertices projections on Free Surface
            pi_1 = {p_1(0), p_1(1), incEta(0)};
            pi_2 = {p_2(0), p_2(1), incEta(1)};
            pi_3 = {p_3(0), p_3(1), incEta(2)};

            v_p = p_3 - p_1;
            v_pi = pi_3 - pi_1;
            if (arma::norm(v_p.subvec(0, 1)) < TOL)
            {
                q_4 = pi_3;
            }
            else
            {
                lambda = (p_1(2) - pi_1(2)) / (v_pi(2) - v_p(2));
                q_4 = p_1 + lambda * v_p;
            }

            v_p = p_3 - p_2;
            v_pi = pi_3 - pi_2;
            if (arma::norm(v_p.subvec(0, 1)) < TOL)
            {
                q_3 = pi_3;
            }
            else
            {
                lambda = (p_2(2) - pi_2(2)) / (v_pi(2) - v_p(2));
                q_3 = p_2 + lambda * v_p;
            }
        }

        // New vertex nodes
        p_3 = q_3;
        p_4 = q_4;

        // Other nodes
        v_1 = p_2 - p_1;
        v_2 = p_2 - p_3;
        v_3 = p_4 - p_3;
        v_4 = p_4 - p_1;
        p_12 = p_1 + LAMBDA * v_1;
        p_23 = p_3 + LAMBDA * v_2;
        p_34 = p_3 + LAMBDA * v_3;
        p_41 = p_1 + LAMBDA * v_4;
        p_123 = (p_12 + p_23) / 2;
        p_134 = (p_34 + p_41) / 2;
        p_1234 = (p_1 + p_3) / 2;

        // Nodes output matrix
        outNodes = arma::zeros(11, 3);
        outNodes.row(0) = p_1;
        outNodes.row(1) = p_2;
        outNodes.row(2) = p_3;
        outNodes.row(3) = p_4;
        outNodes.row(4) = p_12;
        outNodes.row(5) = p_23;
        outNodes.row(6) = p_34;
        outNodes.row(7) = p_41;
        outNodes.row(8) = p_123;
        outNodes.row(9) = p_134;
        outNodes.row(10) = p_1234;

        // Jacobian computation
        jac1 = arma::norm(arma::cross(v_1, v_2)) / 4;
        jac2 = arma::norm(arma::cross(v_3, v_4)) / 4;
        jac = {jac1, jac2};

        // --------------------- Cutting element with two vertices underewater ---------------------

        // Extracting the normal of the original element
        tmpNormal = transNormals.row(ielem);

        // Copy elem info
        nodes = arma::join_vert(nodes, outNodes);
        normals = arma::join_vert(normals, tmpNormal, tmpNormal);
        auxVec1 = {2, 0, 1, 10, 4, 5, 8};
        auxVec2 = {0, 2, 3, 10, 6, 7, 9};
        elems = arma::join_vert(elems, numNodes + auxVec1, numNodes + auxVec2);
        jacobians = arma::join_vert(jacobians, jac);
        numNodes += 11;
        numElems += 2;
    }

    // Including elements with 1 vertex submerged
    for (int ielem : ind1)
    {

        if (pBody->flag_hydrostatics == 1)
        {
            // Sorting vertices by height
            incNodes = sort_rows(transNodes.rows(iniElems(ielem, arma::span(0, 2))), 2);
        }
        else if (pBody->flag_hydrostatics == 2)
        {
            // Sorting vertices by height wrt free surface
            tmpElems = iniElems(ielem, arma::span(0, 2));
            icol = {2};
            sortInd = arma::sort_index(transNodes.submat(tmpElems, icol) - eta(tmpElems));
            incNodes = transNodes.rows(tmpElems(sortInd));
            incEta = eta.rows(tmpElems(sortInd));
        }

        // --------------------- Cutting element with one vertex underewater ---------------------

        // Initial vertex nodes
        p_1 = incNodes.row(0);
        p_2 = incNodes.row(1);
        p_3 = incNodes.row(2);

        // Triangle cut
        if (pBody->flag_hydrostatics == 1)
        {
            w_1 = p_2 - p_1;
            w_2 = p_3 - p_1;
            mu_1 = arma::as_scalar(-p_1(2) / w_1(2));
            mu_2 = arma::as_scalar(-p_1(2) / w_2(2));
            q_2.subvec(0, 1) = p_1.subvec(0, 1) + mu_1 * w_1.subvec(0, 1);
            q_2(2) = 0;
            q_3.subvec(0, 1) = p_1.subvec(0, 1) + mu_2 * w_2.subvec(0, 1);
            q_3(2) = 0;
        }
        else if (pBody->flag_hydrostatics == 2)
        {
            // Vectices proyections on Free Surface
            pi_1 = {p_1(0), p_1(1), incEta(0)};
            pi_2 = {p_2(0), p_2(1), incEta(1)};
            pi_3 = {p_3(0), p_3(1), incEta(2)};

            v_p = p_2 - p_1;
            v_pi = pi_2 - pi_1;
            if (arma::norm(v_p.subvec(0, 1)) < TOL)
            {
                q_2 = pi_2;
            }
            else
            {
                lambda = (p_1(2) - pi_1(2)) / (v_pi(2) - v_p(2));
                q_2 = p_1 + lambda * v_p;
            }

            v_p = p_3 - p_1;
            v_pi = pi_3 - pi_1;
            if (arma::norm(v_p.subvec(0, 1)) < TOL)
            {
                q_3 = pi_3;
            }
            else
            {
                lambda = (p_1(2) - pi_1(2)) / (v_pi(2) - v_p(2));
                q_3 = p_1 + lambda * v_p;
            }
        }

        // New vertex nodes
        p_2 = q_2;
        p_3 = q_3;

        // Other nodes
        v_1 = p_3 - p_1;
        v_2 = p_3 - p_2;
        p_12 = (p_1 + p_2) / 2;
        p_23 = p_2 + LAMBDA * v_2;
        p_31 = p_1 + LAMBDA * v_1;
        p_123 = (p_23 + p_31) / 2;

        // Nodes output matrix
        outNodes = arma::zeros(7, 3);
        outNodes.row(0) = p_1;
        outNodes.row(1) = p_2;
        outNodes.row(2) = p_3;
        outNodes.row(3) = p_12;
        outNodes.row(4) = p_23;
        outNodes.row(5) = p_31;
        outNodes.row(6) = p_123;

        // Jacobian computation
        jac = {arma::norm(arma::cross(v_1, v_2)) / 4};

        // --------------------- Cutting element with one vertex underewater ---------------------

        // Extracting the normal of the original element
        tmpNormal = transNormals.row(ielem);

        // Copy elem info
        nodes = arma::join_vert(nodes, outNodes);
        normals = arma::join_vert(normals, tmpNormal);
        auxVec = {0, 1, 2, 3, 4, 5, 6};
        elems = arma::join_vert(elems, numNodes + auxVec);
        jacobians = arma::join_vert(jacobians, jac);
        numNodes += 7;
        numElems += 1;
    }
}

void BodyMesh::IntegrateMesh(void)
{
    // std::cout << "--> Integrating Body 2D Mesh" << std::endl;

    weightsJacNormal = arma::zeros(numNodes, 3);
    arma::vec weightsVector = {1.0 / 12.0, 1.0 / 12.0, 1.0 / 9.0, 1.0 / 3.0, 25.0 / 108.0, 25.0 / 108.0, 25.0 / 27.0};
    arma::uvec tmpInd;

    // Assembling the integration vectors with the weights, jacobians and normals
    for (int ielem = 0; ielem < numElems; ielem++)
    {
        tmpInd = elems.row(ielem).t();
        weightsJacNormal.rows(tmpInd) += weightsVector * jacobians(ielem) * normals.row(ielem);
    }
}

////////////////////////////////////////////////////////////////////////////
///////////////////////// BodyTri2DMesh CLASS DEFINITION ///////////////////
////////////////////////////////////////////////////////////////////////////

int BodyTri2DMesh::GetType(void)
{
    return this->typeMesh;
}

void BodyTri2DMesh::Preprocess(void)
{
    Logger::info("--> Preprocessing Body 2D Mesh...");

    double LAMBDA = 0.4; // RATIO FOR SOME GAUSSIAN NODES

    // Extract data from non-preprocessed mesh
    arma::mat vertices = iniNodes;
    arma::umat elems_vertices = iniElems;
    numVertices = iniNumNodes;

    // Variables initialization
    arma::mat tmpEdges = arma::zeros(iniNumElems * 3, 3);
    arma::mat faces = arma::zeros(iniNumElems, 3);
    arma::umat tmp_elems_edges = arma::zeros<arma::umat>(iniNumElems, 3);
    arma::uvec elems_faces = arma::zeros<arma::uvec>(iniNumElems);
    arma::vec edges_length = arma::zeros(iniNumElems * 3);
    iniNormals = arma::zeros(iniNumElems, 3);
    iniJacobians = arma::zeros(iniNumElems);

    // Elements gaussian nodes computation
    arma::rowvec p_1, p_2, p_3, p_12, p_23, p_31, p_123, v_1, v_2, tmpNormal;
    double norm_fN;

    for (unsigned int iface = 0; iface < iniNumElems; iface++)
    {

        p_1 = vertices.row(elems_vertices(iface, 0));
        p_2 = vertices.row(elems_vertices(iface, 1));
        p_3 = vertices.row(elems_vertices(iface, 2));

        v_1 = p_3 - p_1;
        v_2 = p_3 - p_2;
        p_12 = (p_1 + p_2) / 2;
        p_23 = p_2 + LAMBDA * v_2;
        p_31 = p_1 + LAMBDA * v_1;
        p_123 = (p_23 + p_31) / 2;

        tmpEdges.row(3 * iface) = p_12;
        tmpEdges.row(3 * iface + 1) = p_23;
        tmpEdges.row(3 * iface + 2) = p_31;

        faces.row(iface) = p_123;

        tmp_elems_edges.row(iface) = {3 * iface, 3 * iface + 1, 3 * iface + 2};
        elems_faces(iface) = iface;

        tmpNormal = arma::cross(v_1, v_2);
        norm_fN = arma::norm(tmpNormal);
        iniNormals.row(iface) = tmpNormal / norm_fN;
        iniJacobians(iface) = norm_fN / 4;

        edges_length(3 * iface) = arma::norm(p_1 - p_2);
        edges_length(3 * iface + 1) = arma::norm(p_2 - p_3);
        edges_length(3 * iface + 2) = arma::norm(p_3 - p_1);
    }

    // Duplicated edges nodes removal
    arma::mat edges;
    arma::uvec ind;
    std::tie(edges, ind) = unique_rows(tmpEdges);
    int numEdges = edges.n_rows;

    // Rearrange elems_edges index
    arma::umat elems_edges = indMat(ind, tmp_elems_edges);

    // Reorder elements
    elems_edges += numVertices;
    elems_faces += numVertices + numEdges;

    // Putting together variables
    iniNodes = arma::join_vert(vertices, edges, faces);
    iniNumNodes = numVertices + numEdges + iniNumElems;
    iniElems = arma::join_horiz(elems_vertices, elems_edges, elems_faces);

    // Maximum edges length
    maxEdgesLength = arma::max(edges_length);
    if (pBody->flag_hydrostatics == 2 && pSim->pWave->lambda_peak < 8 * maxEdgesLength)
    {
        Logger::warning("The mesh is too coarse for the selected waves.");
    }
}
