#include "PDNNetwork.h"
#include "PDNLoadMapper.h"
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
              << "  pdn_load_example --flp <path> --config <path> [options]\n"
              << "\n"
              << "Required:\n"
              << "  --flp <path>          Floorplan file (FLP format)\n"
              << "                        Example: ../../VoltSpot/example.flp\n"
              << "  --config <path>       VoltSpot config file\n"
              << "                        Example: ../../VoltSpot/pdn.config\n"
              << "\n"
              << "Optional:\n"
              << "  --padloc <path>       Pad location file (for voltage sources)\n"
              << "                        Example: ../../VoltSpot/example.vgrid.padloc\n"
              << "  --ptrace <path>       Power trace file (for current loads)\n"
              << "                        Example: ../../VoltSpot/example.ptrace\n"
              << "  --mlcf <path>         Metal layer config file\n"
              << "                        Example: ../../VoltSpot/example.mlcf\n"
              << "  --vdd <voltage>       VDD voltage (default: 1.0)\n"
              << "  --gnd <voltage>       GND voltage (default: 0.0)\n"
              << "  --padloc_format <n>   Pad location format: 0/1=virtual grid, 2/3=pad grid (default: 0)\n"
              << "  --out <csv>           Output CSV file for branch data (default: pdn_with_loads_branches.csv)\n"
              << "  -h, --help            Show this help message\n"
              << "\n"
              << "Example (from build directory):\n"
              << "  ./pdn_load_example \\\n"
              << "    --flp ../../VoltSpot/example.flp \\\n"
              << "    --config ../../VoltSpot/pdn.config \\\n"
              << "    --padloc ../../VoltSpot/example.vgrid.padloc \\\n"
              << "    --ptrace ../../VoltSpot/example.ptrace\n";
}

int main(int argc, char** argv) {
    std::string flpPath;
    std::string cfgPath;
    std::string padlocPath;
    std::string ptracePath;
    std::string mlcfPath;
    std::string outCsv = "../out_voltspot/pdn_with_loads_branches.csv";
    double vddVoltage = 1.0;
    double gndVoltage = 0.0;
    int padlocFormat = 0;

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
        } else if (arg == "--padloc") {
            padlocPath = needValue(arg);
        } else if (arg == "--ptrace") {
            ptracePath = needValue(arg);
        } else if (arg == "--mlcf") {
            mlcfPath = needValue(arg);
        } else if (arg == "--vdd") {
            vddVoltage = std::stod(needValue(arg));
        } else if (arg == "--gnd") {
            gndVoltage = std::stod(needValue(arg));
        } else if (arg == "--padloc_format") {
            padlocFormat = std::stoi(needValue(arg));
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
        std::cerr << "Error: --flp and --config are required\n";
        usage();
        return 2;
    }

    try {
        // Step 1: Parse chip size from FLP file
        std::cout << "Step 1: Parsing chip size from FLP file..." << std::endl;
        const auto chip = VoltSpotVirtualGrid::parseVoltSpotFlpChipSize(flpPath);
        std::cout << "  Chip size: " << chip.width_m << " m x " << chip.height_m << " m\n";

        // Step 2: Parse config file
        std::cout << "\nStep 2: Parsing config file..." << std::endl;
        const auto kv = VoltSpotVirtualGrid::parseVoltSpotConfigFile(cfgPath);
        auto cfg = VoltSpotVirtualGrid::configFromMap(kv);
        if (!mlcfPath.empty()) {
            cfg.mlcfPath = mlcfPath;
        }

        // Step 3: Build virtual grid
        std::cout << "\nStep 3: Building virtual grid..." << std::endl;
        auto gridResult = VoltSpotVirtualGrid::build(chip, cfg);
        std::cout << "  Grid: " << gridResult.grid.rows << " rows x " << gridResult.grid.cols << " cols\n";
        std::cout << "  Grid spacing: dx=" << gridResult.grid.dx_m << " m, dy=" << gridResult.grid.dy_m << " m\n";
        std::cout << "  Edge resistances: rx=" << gridResult.r.rx_edge_ohm << " ohm, ry=" << gridResult.r.ry_edge_ohm << " ohm\n";

        // Step 4: Create PDN network with pixel models
        std::cout << "\nStep 4: Creating PDN network..." << std::endl;
        PDNNetwork network(gridResult.grid.rows, gridResult.grid.cols, 1);
        std::vector<std::vector<std::unique_ptr<PixelModel>>> pixels;
        pixels.resize(static_cast<size_t>(gridResult.grid.rows));
        for (int r = 0; r < gridResult.grid.rows; r++) {
            pixels[static_cast<size_t>(r)].resize(static_cast<size_t>(gridResult.grid.cols));
            for (int c = 0; c < gridResult.grid.cols; c++) {
                const std::string name = "pixel_r" + std::to_string(r) + "_c" + std::to_string(c);
                pixels[static_cast<size_t>(r)][static_cast<size_t>(c)] =
                    std::make_unique<PixelModel>(name, static_cast<float>(gridResult.r.rx_edge_ohm),
                                                 static_cast<float>(gridResult.r.ry_edge_ohm), 0.0f, 0.0f);
                network.setPixelModel(0, r, c, pixels[static_cast<size_t>(r)][static_cast<size_t>(c)].get());
            }
        }

        // Step 5: Build network structure
        std::cout << "\nStep 5: Building network structure..." << std::endl;
        network.buildNetwork();
        network.printNetworkInfo();

        // Step 6: Load voltage sources from padloc file
        if (!padlocPath.empty()) {
            std::cout << "\nStep 6: Loading voltage sources from padloc file..." << std::endl;
            if (!std::filesystem::exists(padlocPath)) {
                std::cerr << "Warning: Padloc file not found: " << padlocPath << std::endl;
            } else {
                PDNLoadMapper::loadVoltageSources(network, padlocPath, gridResult, 
                                                  vddVoltage, gndVoltage, padlocFormat);
            }
        } else {
            std::cout << "\nStep 6: Skipping voltage sources (--padloc not specified)" << std::endl;
        }

        // Step 7: Load current loads from FLP and PTrace files
        if (!ptracePath.empty()) {
            std::cout << "\nStep 7: Loading current loads from FLP and PTrace files..." << std::endl;
            if (!std::filesystem::exists(ptracePath)) {
                std::cerr << "Warning: PTrace file not found: " << ptracePath << std::endl;
            } else {
                PDNLoadMapper::loadCurrentLoads(network, flpPath, ptracePath, gridResult, vddVoltage);
            }
        } else {
            std::cout << "\nStep 7: Skipping current loads (--ptrace not specified)" << std::endl;
        }

        // Step 8: Export results
        std::cout << "\nStep 8: Exporting branch and node data..." << std::endl;
        const std::filesystem::path outPath(outCsv);
        if (outPath.has_parent_path()) {
            std::filesystem::create_directories(outPath.parent_path());
        }
        network.exportBranchDataForVisualization(outCsv);
        
        // Export node data (for visualization of voltage sources and current loads)
        std::filesystem::path nodeCsvPath = outPath;
        nodeCsvPath.replace_filename(outPath.stem().string() + "_nodes.csv");
        network.exportNodeDataForVisualization(nodeCsvPath.string());

        // Summary
        std::cout << "\n========================================" << std::endl;
        std::cout << "Summary" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Grid size: " << gridResult.grid.rows << " x " << gridResult.grid.cols << std::endl;
        std::cout << "Total nodes: " << network.getNodeCount() << std::endl;
        std::cout << "Total branches: " << network.getBranchCount() << std::endl;
        
        // Count voltage sources and current loads
        int voltageSourceCount = 0;
        int currentLoadCount = 0;
        auto nodes = network.getNodes();
        for (const auto& node : nodes) {
            if (node.voltage != 0.0) {
                voltageSourceCount++;
            }
            if (node.currentLoad != 0.0) {
                currentLoadCount++;
            }
        }
        std::cout << "Nodes with voltage source: " << voltageSourceCount << std::endl;
        std::cout << "Nodes with current load: " << currentLoadCount << std::endl;
        std::cout << "Branch CSV: " << std::filesystem::absolute(outCsv) << std::endl;
        std::cout << "Node CSV: " << std::filesystem::absolute(nodeCsvPath) << std::endl;
        std::cout << "========================================" << std::endl;

        std::cout << "\nDone!" << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

