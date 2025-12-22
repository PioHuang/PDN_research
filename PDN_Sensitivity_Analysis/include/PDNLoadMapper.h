#ifndef PDN_LOAD_MAPPER_H
#define PDN_LOAD_MAPPER_H

#include "PDNNetwork.h"
#include "VoltSpotVirtualGrid.h"
#include <string>

/**
 * Utility class for mapping voltage sources and current loads
 * from VoltSpot input files to PDNNetwork
 */
class PDNLoadMapper {
public:
    /**
     * Load voltage sources from padloc file
     * @param network Target PDN network
     * @param padlocPath Path to padloc file
     * @param gridResult Virtual grid configuration
     * @param vddVoltage VDD voltage value
     * @param gndVoltage GND voltage value (usually 0.0)
     * @param padlocFormat 0 or 1: virtual grid coordinates, 2 or 3: pad grid coordinates
     */
    static void loadVoltageSources(PDNNetwork& network, const std::string& padlocPath,
                                   const VoltSpotVirtualGrid::Result& gridResult,
                                   double vddVoltage = 1.0, double gndVoltage = 0.0,
                                   int padlocFormat = 0);

    /**
     * Load current loads from FLP and PTrace files
     * @param network Target PDN network
     * @param flpPath Path to FLP (floorplan) file
     * @param ptracePath Path to PTrace (power trace) file
     * @param gridResult Virtual grid configuration
     * @param supplyVoltage Supply voltage for current calculation (P = V * I)
     */
    static void loadCurrentLoads(PDNNetwork& network, const std::string& flpPath,
                                const std::string& ptracePath,
                                const VoltSpotVirtualGrid::Result& gridResult,
                                double supplyVoltage = 1.0);

    /**
     * Allocate P/G pads to all available pad seats (VoltSpot behavior when PDN_padconfig != 0).
     * Uses a simple interleaved checkerboard pattern on the pad grid.
     */
    static void allocatePadsAllSeats(PDNNetwork& network,
                                     const VoltSpotVirtualGrid::Result& gridResult,
                                     double vddVoltage = 1.0,
                                     double gndVoltage = 0.0);
};

#endif // PDN_LOAD_MAPPER_H

