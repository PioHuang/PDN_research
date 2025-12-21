#include "PDNNetwork.h"
#include "PixelModel.h"
#include <iostream>
#include <memory>
#include <vector>

/**
 * Simple example demonstrating how to build a PDN (Power Delivery Network) model
 */
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "PDN Model Building Example" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Define grid dimensions
    int rows = 5;
    int cols = 5;
    int layers = 2;
    
    std::cout << "\nCreating PDN network: " << rows << "x" << cols 
              << " grid with " << layers << " layer(s)" << std::endl;
    
    // Create PDN network
    PDNNetwork network(rows, cols, layers);
    
    // Create pixel models for each grid point
    // Each pixel model has resistances: Rx (horizontal), Ry (vertical), Rz (between layers)
    std::vector<std::vector<std::vector<std::unique_ptr<PixelModel>>>> pixelModels(layers);
    
    // Initialize pixel models with resistance values
    // Example values: Rx = 0.1 ohm, Ry = 0.1 ohm, Rz = 0.5 ohm
    float baseRx = 0.1f;
    float baseRy = 0.1f;
    float baseRz = 0.5f;
    
    for (int l = 0; l < layers; l++) {
        pixelModels[l].resize(rows);
        for (int r = 0; r < rows; r++) {
            pixelModels[l][r].resize(cols);
            for (int c = 0; c < cols; c++) {
                std::string name = "pixel_L" + std::to_string(l) + "_R" + std::to_string(r) + "_C" + std::to_string(c);
                pixelModels[l][r][c] = std::make_unique<PixelModel>(name, baseRx, baseRy, baseRz, 1.0f);
                network.setPixelModel(l, r, c, pixelModels[l][r][c].get());
            }
        }
    }
    
    std::cout << "\nPixel models created and assigned to network" << std::endl;
    std::cout << "  Rx (horizontal): " << baseRx << " ohm" << std::endl;
    std::cout << "  Ry (vertical): " << baseRy << " ohm" << std::endl;
    std::cout << "  Rz (between layers): " << baseRz << " ohm" << std::endl;
    
    // Build the network (creates nodes and branches from pixel models)
    std::cout << "\nBuilding network..." << std::endl;
    network.buildNetwork();
    
    // Print network information
    std::cout << "\n";
    network.printNetworkInfo();
    
    // Example: Set a voltage source at one corner
    std::cout << "\nSetting voltage source at (0, 0, 0) = 1.0V" << std::endl;
    network.setVoltageSource(0, 0, 0, 1.0, "VDD");
    
    // Example: Set a current load at another location
    std::cout << "Setting current load at (" << (rows-1) << ", " << (cols-1) << ", 0) = 0.1A" << std::endl;
    network.setCurrentLoad(layers-1, rows-1, cols-1, 0.1, "VDD");
    
    // Get network statistics
    std::cout << "\nNetwork Statistics:" << std::endl;
    std::cout << "  Total nodes: " << network.getNodeCount() << std::endl;
    std::cout << "  Total branches: " << network.getBranchCount() << std::endl;
    
    // Example: Access branch information
    auto branches = network.getBranches();
    if (!branches.empty()) {
        std::cout << "\nFirst few branches:" << std::endl;
        int count = 0;
        for (const auto& branch : branches) {
            if (count++ < 5) {
                std::cout << "  Branch " << branch.node1Id << " -> " << branch.node2Id 
                          << " (R=" << branch.resistance << " ohm, dir=" << branch.direction << ")" << std::endl;
            }
        }
        if (branches.size() > 5) {
            std::cout << "  ... and " << (branches.size() - 5) << " more branches" << std::endl;
        }
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Example completed successfully!" << std::endl;
    std::cout << "========================================" << std::endl;
    
    return 0;
}

