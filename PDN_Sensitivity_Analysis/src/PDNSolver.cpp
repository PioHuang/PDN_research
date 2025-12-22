#include "PDNSolver.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <sstream>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using namespace Eigen;

PDNSolver::SolverResult PDNSolver::solveSteadyState(PDNNetwork &network,
                                                    double vddVoltage,
                                                    double gndVoltage,
                                                    double padResistance)
{
    SolverResult result;
    result.success = false;

    try
    {
        // Build node coordinate map (nodeId -> (layer, row, col, pgNet))
        std::map<int, std::tuple<int, int, int, std::string>> nodeIdToCoords;
        buildNodeCoordinateMap(network, nodeIdToCoords);

        // Build pad location map (for voltage sources)
        // Note: VDD pads have voltage != 0.0, but GND pads have voltage == 0.0
        // So we need to check both VDD and GND nodes separately
        std::map<std::tuple<int, int, int>, bool> padLocations;
        const auto nodes = network.getNodes();

        // First, identify VDD pads (voltage != 0.0)
        for (const auto &node : nodes)
        {
            if (node.pgNetName == "VDD" && node.voltage != 0.0)
            {
                auto it = nodeIdToCoords.find(node.id);
                if (it != nodeIdToCoords.end())
                {
                    int layer = std::get<0>(it->second);
                    int row = std::get<1>(it->second);
                    int col = std::get<2>(it->second);
                    padLocations[{layer, row, col}] = true;
                }
            }
        }

        // Then, identify GND pads (voltage == 0.0 but explicitly set as pad)
        // Heuristic: if there's a VDD pad at the same location, assume GND is also a pad
        // This works because VDD and GND pads are typically co-located
        for (const auto &node : nodes)
        {
            if (node.pgNetName == "GND")
            {
                auto it = nodeIdToCoords.find(node.id);
                if (it != nodeIdToCoords.end())
                {
                    int layer = std::get<0>(it->second);
                    int row = std::get<1>(it->second);
                    int col = std::get<2>(it->second);
                    // If there's already a pad at this location (from VDD), mark it
                    // OR if it's on layer 0 and we want to be more inclusive,
                    // we could mark all GND nodes on layer 0 as pads, but that's too broad
                    // Better: check if there's a VDD pad at the same location
                    if (padLocations.find({layer, row, col}) != padLocations.end())
                    {
                        // Already marked by VDD pad, keep it
                    }
                    else if (layer == 0)
                    {
                        // Check if there's a VDD pad at the same location
                        auto *vddNode = network.getNode(layer, row, col, "VDD");
                        if (vddNode && vddNode->voltage != 0.0)
                        {
                            padLocations[{layer, row, col}] = true;
                        }
                    }
                }
            }
        }

        // Calculate matrix size (VoltSpot style: 2 * nl * nr * nc)
        const int nl = network.getLayerCount();
        const int nr = network.getLayerRows(0);
        const int nc = network.getLayerCols(0);
        const int numNodes = 2 * nl * nr * nc;

        // Build conductance matrix
        std::vector<int> rowIndices, colIndices;
        std::vector<double> values;
        buildConductanceMatrix(network, nodeIdToCoords, padLocations, padResistance,
                               rowIndices, colIndices, values);

        // Build RHS vector
        std::vector<double> rhs(numNodes, 0.0);
        buildRHSVector(network, nodeIdToCoords, padLocations, vddVoltage, gndVoltage, padResistance, rhs);

        // Build Eigen sparse matrix
        SparseMatrix<double> G(numNodes, numNodes);
        G.reserve(values.size());

        for (size_t i = 0; i < values.size(); i++)
        {
            G.insert(rowIndices[i], colIndices[i]) = values[i];
        }
        G.makeCompressed();

        // Build RHS vector
        VectorXd I = Map<VectorXd>(rhs.data(), numNodes);

        // Solve using SparseLU
        SparseLU<SparseMatrix<double>> solver;
        solver.compute(G);

        if (solver.info() != Success)
        {
            std::ostringstream oss;
            oss << "Matrix factorization failed (info=" << solver.info() << ")";
            oss << ". Matrix size: " << numNodes << "x" << numNodes;
            oss << ", Non-zeros: " << G.nonZeros();
            result.errorMessage = oss.str();
            return result;
        }

        VectorXd V = solver.solve(I);

        if (solver.info() != Success)
        {
            std::ostringstream oss;
            oss << "Linear system solve failed (info=" << solver.info() << ")";
            result.errorMessage = oss.str();
            return result;
        }

        // Copy results
        result.nodeVoltages.resize(numNodes);
        for (int i = 0; i < numNodes; i++)
        {
            result.nodeVoltages[i] = V[i];
        }

        // Calculate IR drop statistics (VoltSpot style: based on VDD-GND voltage difference)
        double maxDrop = 0.0;
        double sumDrop = 0.0;
        int dropCount = 0;

        for (int l = 0; l < nl; l++)
        {
            for (int r = 0; r < nr; r++)
            {
                for (int c = 0; c < nc; c++)
                {
                    int vddIdx = l * nr * nc + r * nc + c;                // VDD node index
                    int gndIdx = nl * nr * nc + l * nr * nc + r * nc + c; // GND node index

                    if (vddIdx < static_cast<int>(result.nodeVoltages.size()) &&
                        gndIdx < static_cast<int>(result.nodeVoltages.size()))
                    {
                        double vddVolt = result.nodeVoltages[vddIdx];
                        double gndVolt = result.nodeVoltages[gndIdx];
                        double voltageDiff = vddVolt - gndVolt;     // Actual VDD-GND voltage difference
                        double idealDiff = vddVoltage - gndVoltage; // Ideal voltage difference
                        double drop = 100.0 * (1.0 - voltageDiff / idealDiff);
                        if (drop > maxDrop)
                            maxDrop = drop;
                        sumDrop += drop;
                        dropCount++;
                    }
                }
            }
        }

        result.maxIRDrop = maxDrop;
        result.avgIRDrop = dropCount > 0 ? sumDrop / dropCount : 0.0;
        result.success = true;
    }
    catch (const std::exception &e)
    {
        result.errorMessage = std::string("Exception: ") + e.what();
    }

    return result;
}

void PDNSolver::buildNodeCoordinateMap(const PDNNetwork &network,
                                       std::map<int, std::tuple<int, int, int, std::string>> &nodeIdToCoords)
{
    nodeIdToCoords.clear();

    PDNNetwork &nonConstNetwork = const_cast<PDNNetwork &>(network);

    for (int l = 0; l < network.getLayerCount(); l++)
    {
        for (int r = 0; r < network.getLayerRows(l); r++)
        {
            for (int c = 0; c < network.getLayerCols(l); c++)
            {
                for (const std::string &pgNet : {"VDD", "GND"})
                {
                    auto *node = nonConstNetwork.getNode(l, r, c, pgNet);
                    if (node)
                    {
                        nodeIdToCoords[node->id] = std::make_tuple(l, r, c, pgNet);
                    }
                }
            }
        }
    }
}

int PDNSolver::nodeIdToVoltSpotIndex(int nodeId,
                                     const PDNNetwork &network,
                                     const std::map<int, std::tuple<int, int, int, std::string>> &nodeIdToCoords)
{
    auto it = nodeIdToCoords.find(nodeId);
    if (it == nodeIdToCoords.end())
    {
        return -1;
    }

    int layer = std::get<0>(it->second);
    int row = std::get<1>(it->second);
    int col = std::get<2>(it->second);
    std::string pgNet = std::get<3>(it->second);

    const int nl = network.getLayerCount();
    const int nr = network.getLayerRows(0);
    const int nc = network.getLayerCols(0);

    if (pgNet == "VDD")
    {
        return layer * nr * nc + row * nc + col;
    }
    else
    { // GND
        return nl * nr * nc + layer * nr * nc + row * nc + col;
    }
}

void PDNSolver::buildConductanceMatrix(const PDNNetwork &network,
                                       const std::map<int, std::tuple<int, int, int, std::string>> &nodeIdToCoords,
                                       const std::map<std::tuple<int, int, int>, bool> &padLocations,
                                       double padResistance,
                                       std::vector<int> &rowIndices,
                                       std::vector<int> &colIndices,
                                       std::vector<double> &values)
{
    rowIndices.clear();
    colIndices.clear();
    values.clear();

    const int nl = network.getLayerCount();
    const int nr = network.getLayerRows(0);
    const int nc = network.getLayerCols(0);
    const auto branches = network.getBranches();

    // Build matrix entries: map from (row, col) to conductance value
    std::map<std::pair<int, int>, double> matrixEntries;

    // Initialize all diagonal elements to zero (to ensure all nodes have entries)
    for (int i = 0; i < 2 * nl * nr * nc; i++)
    {
        matrixEntries[{i, i}] = 0.0;
    }

    // Process each branch to build conductance matrix
    // Note: PDNNetwork only builds VDD branches, but VoltSpot has separate GND grid
    // We'll process VDD branches and create corresponding GND branches with same resistances
    for (const auto &branch : branches)
    {
        int vsIdx1 = nodeIdToVoltSpotIndex(branch.node1Id, network, nodeIdToCoords);
        int vsIdx2 = nodeIdToVoltSpotIndex(branch.node2Id, network, nodeIdToCoords);

        if (vsIdx1 < 0 || vsIdx2 < 0)
            continue;
        if (branch.resistance <= 0.0)
            continue;

        double conductance = 1.0 / branch.resistance;

        // Off-diagonal: G[i][j] = -1/R
        matrixEntries[{vsIdx1, vsIdx2}] -= conductance;
        matrixEntries[{vsIdx2, vsIdx1}] -= conductance;

        // Diagonal: G[i][i] += 1/R (accumulate)
        matrixEntries[{vsIdx1, vsIdx1}] += conductance;
        matrixEntries[{vsIdx2, vsIdx2}] += conductance;
    }

    // Create GND branches (VoltSpot has separate GND grid with same topology)
    // Use same resistances as VDD (from PixelModel)
    for (int l = 0; l < nl; l++)
    {
        for (int r = 0; r < nr; r++)
        {
            for (int c = 0; c < nc; c++)
            {
                PixelModel *model = network.getPixelModel(l, r, c);
                if (!model)
                    continue;

                // GND node indices (offset by nl*nr*nc)
                int gndBaseIdx = nl * nr * nc;
                int gndIdx1 = gndBaseIdx + l * nr * nc + r * nc + c;

                // Horizontal (x-direction) GND connections
                if (c < nc - 1)
                {
                    int gndIdx2 = gndBaseIdx + l * nr * nc + r * nc + (c + 1);
                    double conductance = 1.0 / model->getRx();
                    matrixEntries[{gndIdx1, gndIdx2}] -= conductance;
                    matrixEntries[{gndIdx2, gndIdx1}] -= conductance;
                    matrixEntries[{gndIdx1, gndIdx1}] += conductance;
                    matrixEntries[{gndIdx2, gndIdx2}] += conductance;
                }

                // Vertical (y-direction) GND connections
                if (r < nr - 1)
                {
                    int gndIdx2 = gndBaseIdx + l * nr * nc + (r + 1) * nc + c;
                    double conductance = 1.0 / model->getRy();
                    matrixEntries[{gndIdx1, gndIdx2}] -= conductance;
                    matrixEntries[{gndIdx2, gndIdx1}] -= conductance;
                    matrixEntries[{gndIdx1, gndIdx1}] += conductance;
                    matrixEntries[{gndIdx2, gndIdx2}] += conductance;
                }

                // Vertical between layers (z-direction) GND connections
                if (l < nl - 1)
                {
                    int gndIdx2 = gndBaseIdx + (l + 1) * nr * nc + r * nc + c;
                    double conductance = 1.0 / model->getRz();
                    matrixEntries[{gndIdx1, gndIdx2}] -= conductance;
                    matrixEntries[{gndIdx2, gndIdx1}] -= conductance;
                    matrixEntries[{gndIdx1, gndIdx1}] += conductance;
                    matrixEntries[{gndIdx2, gndIdx2}] += conductance;
                }
            }
        }
    }

    // Add pad resistance contributions (VoltSpot style: add 1/Rp to diagonal for pad nodes)
    for (const auto &padEntry : padLocations)
    {
        int layer = std::get<0>(padEntry.first);
        int row = std::get<1>(padEntry.first);
        int col = std::get<2>(padEntry.first);

        // Only add pad resistance to VDD nodes on layer 0 (VoltSpot convention)
        if (layer == 0 && padResistance > 0.0)
        {
            int vsIdx = layer * nr * nc + row * nc + col; // VDD node
            matrixEntries[{vsIdx, vsIdx}] += 1.0 / padResistance;
        }
    }

    // Convert to COO format (only non-zero entries)
    for (const auto &entry : matrixEntries)
    {
        // Only add non-zero entries (or very small values that might be numerical noise)
        if (std::abs(entry.second) > 1e-15)
        {
            rowIndices.push_back(entry.first.first);
            colIndices.push_back(entry.first.second);
            values.push_back(entry.second);
        }
    }
}

void PDNSolver::buildRHSVector(const PDNNetwork &network,
                               const std::map<int, std::tuple<int, int, int, std::string>> &nodeIdToCoords,
                               const std::map<std::tuple<int, int, int>, bool> &padLocations,
                               double vddVoltage,
                               double gndVoltage,
                               double padResistance,
                               std::vector<double> &rhs)
{
    const int nl = network.getLayerCount();
    const int nr = network.getLayerRows(0);
    const int nc = network.getLayerCols(0);
    const int numNodes = 2 * nl * nr * nc;
    rhs.resize(numNodes, 0.0);

    const auto nodes = network.getNodes();

    // Build map from VoltSpot index to node
    std::map<int, const PDNNetwork::Node *> vsIndexToNode;
    for (const auto &node : nodes)
    {
        int vsIdx = nodeIdToVoltSpotIndex(node.id, network, nodeIdToCoords);
        if (vsIdx >= 0)
        {
            vsIndexToNode[vsIdx] = &node;
        }
    }

    // Build RHS following VoltSpot's approach
    for (int l = 0; l < nl; l++)
    {
        for (int r = 0; r < nr; r++)
        {
            for (int c = 0; c < nc; c++)
            {
                // VDD nodes
                int vsIdxVDD = l * nr * nc + r * nc + c;
                auto itVDD = vsIndexToNode.find(vsIdxVDD);
                if (itVDD != vsIndexToNode.end())
                {
                    const auto *node = itVDD->second;
                    bool hasPad = (l == 0 && padLocations.find({l, r, c}) != padLocations.end());

                    if (hasPad && padResistance > 0.0)
                    {
                        // VoltSpot: rhs = vdd/Rp - power/(vdd-gnd)
                        rhs[vsIdxVDD] = vddVoltage / padResistance - node->currentLoad;
                    }
                    else
                    {
                        // Regular node: rhs = -power/(vdd-gnd) = -currentLoad
                        rhs[vsIdxVDD] = -node->currentLoad;
                    }
                }

                // GND nodes
                int vsIdxGND = nl * nr * nc + l * nr * nc + r * nc + c;
                auto itGND = vsIndexToNode.find(vsIdxGND);
                if (itGND != vsIndexToNode.end())
                {
                    const auto *node = itGND->second;
                    bool hasPad = (l == 0 && padLocations.find({l, r, c}) != padLocations.end());

                    if (hasPad && padResistance > 0.0)
                    {
                        // VoltSpot: rhs = gnd/Rp + power/(vdd-gnd)
                        rhs[vsIdxGND] = gndVoltage / padResistance + node->currentLoad;
                    }
                    else
                    {
                        // Regular node: rhs = power/(vdd-gnd) = currentLoad
                        rhs[vsIdxGND] = node->currentLoad;
                    }
                }
            }
        }
    }
}

void PDNSolver::updateNetworkVoltages(PDNNetwork &network, const SolverResult &result)
{
    if (!result.success)
        return;

    std::map<int, std::tuple<int, int, int, std::string>> nodeIdToCoords;
    buildNodeCoordinateMap(network, nodeIdToCoords);

    const int nl = network.getLayerCount();
    const int nr = network.getLayerRows(0);
    const int nc = network.getLayerCols(0);

    // Update VDD node voltages
    for (int l = 0; l < nl; l++)
    {
        for (int r = 0; r < nr; r++)
        {
            for (int c = 0; c < nc; c++)
            {
                int vsIdx = l * nr * nc + r * nc + c;
                if (vsIdx < static_cast<int>(result.nodeVoltages.size()))
                {
                    auto *node = network.getNode(l, r, c, "VDD");
                    if (node && node->voltage == 0.0)
                    { // Only update non-voltage-source nodes
                        network.setVoltageSource(l, r, c, result.nodeVoltages[vsIdx], "VDD");
                    }
                }
            }
        }
    }

    // Update GND node voltages
    for (int l = 0; l < nl; l++)
    {
        for (int r = 0; r < nr; r++)
        {
            for (int c = 0; c < nc; c++)
            {
                int vsIdx = nl * nr * nc + l * nr * nc + r * nc + c;
                if (vsIdx < static_cast<int>(result.nodeVoltages.size()))
                {
                    auto *node = network.getNode(l, r, c, "GND");
                    if (node && node->voltage == 0.0)
                    {
                        network.setVoltageSource(l, r, c, result.nodeVoltages[vsIdx], "GND");
                    }
                }
            }
        }
    }
}

void PDNSolver::exportIRDropResults(const PDNNetwork &network,
                                    const SolverResult &result,
                                    const std::string &filename,
                                    double vddVoltage,
                                    double gndVoltage)
{
    if (!result.success)
    {
        std::cerr << "Cannot export: solver did not succeed" << std::endl;
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    const int nl = network.getLayerCount();
    const int nr = network.getLayerRows(0);
    const int nc = network.getLayerCols(0);

    // VoltSpot format: #Layer:N followed by col row drop_ratio%
    for (int l = 0; l < nl; l++)
    {
        file << "#Layer:" << l << "\n";

        // Export VDD nodes (VoltSpot format: col row drop_ratio%)
        // VoltSpot uses: for(i=0; i<nr; i++) fprintf(fp, "%d\t%d\t%lf\n", j, nr-i-1, drop_ratio);
        // This means: col=j, row=nr-i-1 (inverted Y axis)
        // VoltSpot calculates IR drop based on VDD-GND voltage difference
        for (int i = 0; i < nr; i++)
        { // VoltSpot's loop order
            for (int j = 0; j < nc; j++)
            {
                int r = i; // Our row index
                int c = j; // Our col index
                int vddIdx = l * nr * nc + r * nc + c;
                int gndIdx = nl * nr * nc + l * nr * nc + r * nc + c;

                if (vddIdx < static_cast<int>(result.nodeVoltages.size()) &&
                    gndIdx < static_cast<int>(result.nodeVoltages.size()))
                {
                    double vddVolt = result.nodeVoltages[vddIdx];
                    double gndVolt = result.nodeVoltages[gndIdx];
                    double voltageDiff = vddVolt - gndVolt;     // Actual VDD-GND voltage difference
                    double idealDiff = vddVoltage - gndVoltage; // Ideal voltage difference
                    double dropRatio = 100.0 * (1.0 - voltageDiff / idealDiff);

                    // VoltSpot format: col row drop_ratio% (row = nr-i-1 for inverted Y)
                    file << c << "\t" << (nr - i - 1) << "\t"
                         << std::fixed << std::setprecision(6) << dropRatio << "\n";
                }
            }
        }
        file << "\n";
    }

    file.close();
    std::cout << "IR drop results exported to " << filename << std::endl;
}
