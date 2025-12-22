#ifndef VOLTSPOT_VIRTUAL_GRID_H
#define VOLTSPOT_VIRTUAL_GRID_H

#include <string>
#include <unordered_map>
#include <vector>

class VoltSpotVirtualGrid {
public:
    enum class Direction { X = 0, Y = 1 };

    struct MetalLayerSpec {
        int layerId = -1;
        double pitch_m = 0.0;
        double width_m = 0.0;
        double thickness_m = 0.0;
        double rho_ohm_m = 0.0;
        Direction direction = Direction::X;
    };

    struct ChipSize {
        double width_m = 0.0;
        double height_m = 0.0;
    };

    struct GridSize {
        int rows = 0;
        int cols = 0;
        double dx_m = 0.0;
        double dy_m = 0.0;
    };

    struct ResistanceSummary {
        double rx_edge_ohm = 0.0;
        double ry_edge_ohm = 0.0;
        std::vector<double> perMetal_edgeR_ohm;
    };

    struct Result {
        ChipSize chip;
        GridSize grid;
        double padPitch_m = 0.0;
        int gridIntv = 0;
        std::vector<MetalLayerSpec> metal;
        ResistanceSummary r;
    };

    struct Config {
        double padPitch_m = 0.0;
        double padR_ohm = 10e-3; // VoltSpot default PDN_padR
        int padConfig = 1;       // VoltSpot default PDN_padconfig
        int gridIntv = 0;
        std::string mlcfPath;
    };

    static ChipSize parseVoltSpotFlpChipSize(const std::string& flpPath);
    static std::vector<MetalLayerSpec> parseMlcf(const std::string& mlcfPath);
    static std::vector<MetalLayerSpec> defaultMetalStack();
    static std::unordered_map<std::string, std::string> parseVoltSpotConfigFile(const std::string& configPath);
    static Config configFromMap(const std::unordered_map<std::string, std::string>& kv);

    static Result build(const ChipSize& chip, const Config& cfg);
};

#endif // VOLTSPOT_VIRTUAL_GRID_H

