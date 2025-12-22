#ifndef PDN_SOLVER_H
#define PDN_SOLVER_H

#include "PDNNetwork.h"
#include <vector>
#include <string>
#include <map>
#include <tuple>

/**
 * PDN Solver for IR drop analysis using Eigen
 * Follows VoltSpot's matrix building approach
 */
class PDNSolver {
public:
    struct SolverResult {
        bool success;
        std::string errorMessage;
        std::vector<double> nodeVoltages;  // Voltage for each node (VoltSpot index order)
        double maxIRDrop;                   // Maximum IR drop percentage
        double avgIRDrop;                   // Average IR drop percentage
    };

    /**
     * Solve steady-state IR drop (following VoltSpot's approach)
     * @param network The PDN network to solve
     * @param vddVoltage VDD supply voltage
     * @param gndVoltage GND voltage (usually 0.0)
     * @param padResistance Pad resistance (Rp) for voltage sources
     * @return SolverResult with node voltages and IR drop statistics
     */
    static SolverResult solveSteadyState(PDNNetwork& network, 
                                         double vddVoltage = 1.0, 
                                         double gndVoltage = 0.0,
                                         double padResistance = 0.001);

    /**
     * Update network node voltages after solving
     * @param network The PDN network to update
     * @param result The solver result containing node voltages
     */
    static void updateNetworkVoltages(PDNNetwork& network, const SolverResult& result);

    /**
     * Export IR drop results to file (VoltSpot format compatible)
     * @param network The solved network
     * @param result The solver result
     * @param filename Output filename
     * @param vddVoltage VDD supply voltage
     * @param gndVoltage GND voltage
     */
    static void exportIRDropResults(const PDNNetwork& network, 
                                    const SolverResult& result,
                                    const std::string& filename,
                                    double vddVoltage = 1.0,
                                    double gndVoltage = 0.0);

private:
    /**
     * Convert PDNNetwork nodeId to VoltSpot-style index
     * VoltSpot index: VDD: l*nr*nc + r*nc + c, GND: nl*nr*nc + l*nr*nc + r*nc + c
     */
    static int nodeIdToVoltSpotIndex(int nodeId, 
                                     const PDNNetwork& network,
                                     const std::map<int, std::tuple<int, int, int, std::string>>& nodeIdToCoords);

    /**
     * Build node coordinate map (reverse lookup from nodeId to coordinates)
     */
    static void buildNodeCoordinateMap(const PDNNetwork& network,
                                      std::map<int, std::tuple<int, int, int, std::string>>& nodeIdToCoords);

    /**
     * Build conductance matrix following VoltSpot's approach
     * Matrix layout: [L0v, L1v, ..., Ln-1v, L0g, L1g, ..., Ln-1g]
     */
    static void buildConductanceMatrix(const PDNNetwork& network,
                                      const std::map<int, std::tuple<int, int, int, std::string>>& nodeIdToCoords,
                                      const std::map<std::tuple<int, int, int>, bool>& padLocations,
                                      double padResistance,
                                      std::vector<int>& rowIndices,
                                      std::vector<int>& colIndices,
                                      std::vector<double>& values);

    /**
     * Build RHS vector following VoltSpot's approach
     * RHS = current loads + voltage source contributions
     */
    static void buildRHSVector(const PDNNetwork& network,
                              const std::map<int, std::tuple<int, int, int, std::string>>& nodeIdToCoords,
                              const std::map<std::tuple<int, int, int>, bool>& padLocations,
                              double vddVoltage,
                              double gndVoltage,
                              double padResistance,
                              std::vector<double>& rhs);
};

#endif // PDN_SOLVER_H

