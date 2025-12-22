#include "PDNNetwork.h"
#include "PixelModel.h"
#include "VoltSpotVirtualGrid.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static void usage() {
    std::cout << "Usage:\n"
              << "  voltspot_grid_example --flp <path> --config <pdn.config> [--mlcf <path>] [--out <csv>]\n"
              << "\n"
              << "Notes:\n"
              << "  - This reproduces VoltSpot's virtual grid sizing and metal-to-virtual-grid resistance model\n"
              << "    (up to computing rx/ry edge resistances).\n"
              << "  - Units follow VoltSpot: meters and ohms.\n";
}

int main(int argc, char** argv) {
    std::string flpPath;
    std::string cfgPath;
    std::string mlcfPath;
    std::string outCsv = "voltspot_virtual_grid_branches.csv";

    for (int i = 1; i < argc; i++) {
        const std::string arg = argv[i];
        auto needValue = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for " + name);
            }
            return argv[++i];
        };

        if (arg == "--flp") {
            flpPath = needValue(arg);
        } else if (arg == "--config") {
            cfgPath = needValue(arg);
        } else if (arg == "--mlcf") {
            mlcfPath = needValue(arg);
        } else if (arg == "--out") {
            outCsv = needValue(arg);
        } else if (arg == "-h" || arg == "--help") {
            usage();
            return 0;
        } else {
            std::cerr << "Unknown arg: " << arg << "\n";
            usage();
            return 2;
        }
    }

    if (flpPath.empty() || cfgPath.empty()) {
        usage();
        return 2;
    }

    try {
        const auto chip = VoltSpotVirtualGrid::parseVoltSpotFlpChipSize(flpPath);
        const auto kv = VoltSpotVirtualGrid::parseVoltSpotConfigFile(cfgPath);
        auto cfg = VoltSpotVirtualGrid::configFromMap(kv);
        if (!mlcfPath.empty()) {
            cfg.mlcfPath = mlcfPath;
        }

        auto result = VoltSpotVirtualGrid::build(chip, cfg);

        std::cout << "VoltSpot virtual grid summary\n";
        std::cout << "  chip: " << result.chip.width_m << " m x " << result.chip.height_m << " m\n";
        std::cout << "  padPitch: " << result.padPitch_m << " m\n";
        std::cout << "  gridIntv: " << result.gridIntv << "\n";
        std::cout << "  grid: rows=" << result.grid.rows << " cols=" << result.grid.cols << "\n";
        std::cout << "  dx=" << result.grid.dx_m << " m  dy=" << result.grid.dy_m << " m\n";
        std::cout << "  rx_edge=" << result.r.rx_edge_ohm << " ohm  ry_edge=" << result.r.ry_edge_ohm << " ohm\n";

        // Build a PDNNetwork with uniform per-node pixel models (edge resistance model).
        PDNNetwork network(result.grid.rows, result.grid.cols, 1);
        std::vector<std::vector<std::unique_ptr<PixelModel>>> pixels;
        pixels.resize(static_cast<size_t>(result.grid.rows));
        for (int r = 0; r < result.grid.rows; r++) {
            pixels[static_cast<size_t>(r)].resize(static_cast<size_t>(result.grid.cols));
            for (int c = 0; c < result.grid.cols; c++) {
                const std::string name = "vs_r" + std::to_string(r) + "_c" + std::to_string(c);
                pixels[static_cast<size_t>(r)][static_cast<size_t>(c)] =
                    std::make_unique<PixelModel>(name, static_cast<float>(result.r.rx_edge_ohm),
                                                 static_cast<float>(result.r.ry_edge_ohm), 0.0f, 0.0f);
                network.setPixelModel(0, r, c, pixels[static_cast<size_t>(r)][static_cast<size_t>(c)].get());
            }
        }

        network.buildNetwork();
        network.printNetworkInfo();

        const std::filesystem::path outPath(outCsv);
        if (outPath.has_parent_path()) {
            std::filesystem::create_directories(outPath.parent_path());
        }
        network.exportBranchDataForVisualization(outCsv);
        std::cout << "Done.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
