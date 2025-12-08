#include "PDNNetwork.h"
#include "PixelModel.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

/**
 * C++ version of PDN visualizer
 * Generates DOT format for Graphviz visualization
 */

void generateDOTFile(PDNNetwork& network, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    
    file << "digraph PDN_Network {\n";
    file << "  rankdir=TB;\n";
    file << "  node [shape=circle, style=filled, fillcolor=lightblue];\n";
    file << "  edge [fontsize=8];\n\n";
    
    // Get network dimensions
    int rows = network.getLayerRows(0);
    int cols = network.getLayerCols(0);
    
    // Create nodes
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            std::string nodeName = "N" + std::to_string(r) + "_" + std::to_string(c);
            file << "  " << nodeName << " [label=\"" << r << "," << c << "\"];\n";
        }
    }
    
    file << "\n";
    
    // Create horizontal edges (x-direction)
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols - 1; c++) {
            std::string node1 = "N" + std::to_string(r) + "_" + std::to_string(c);
            std::string node2 = "N" + std::to_string(r) + "_" + std::to_string(c + 1);
            
            PixelModel* model = network.getPixelModel(0, r, c);
            if (model) {
                float rx = model->getRx();
                file << "  " << node1 << " -> " << node2 
                     << " [label=\"Rx=" << std::fixed << std::setprecision(3) << rx << "\", color=blue];\n";
            } else {
                file << "  " << node1 << " -> " << node2 << " [color=gray];\n";
            }
        }
    }
    
    // Create vertical edges (y-direction)
    for (int r = 0; r < rows - 1; r++) {
        for (int c = 0; c < cols; c++) {
            std::string node1 = "N" + std::to_string(r) + "_" + std::to_string(c);
            std::string node2 = "N" + std::to_string(r + 1) + "_" + std::to_string(c);
            
            PixelModel* model = network.getPixelModel(0, r, c);
            if (model) {
                float ry = model->getRy();
                file << "  " << node1 << " -> " << node2 
                     << " [label=\"Ry=" << std::fixed << std::setprecision(3) << ry << "\", color=red];\n";
            } else {
                file << "  " << node1 << " -> " << node2 << " [color=gray];\n";
            }
        }
    }
    
    file << "}\n";
    file.close();
    std::cout << "DOT file generated: " << filename << std::endl;
    std::cout << "To visualize, run: dot -Tpng " << filename << " -o pdn_network.png" << std::endl;
}

void printNetworkASCII(PDNNetwork& network) {
    int rows = network.getLayerRows(0);
    int cols = network.getLayerCols(0);
    
    std::cout << "\n=== PDN Network ASCII Visualization ===\n\n";
    
    // Print nodes and horizontal connections
    for (int r = 0; r < rows; r++) {
        // Print nodes and horizontal branches
        for (int c = 0; c < cols; c++) {
            std::cout << "(" << r << "," << c << ")";
            if (c < cols - 1) {
                PixelModel* model = network.getPixelModel(0, r, c);
                if (model) {
                    std::cout << " --Rx=" << std::fixed << std::setprecision(3) 
                             << model->getRx() << "--> ";
                } else {
                    std::cout << " -----> ";
                }
            }
        }
        std::cout << "\n";
        
        // Print vertical branches
        if (r < rows - 1) {
            for (int c = 0; c < cols; c++) {
                PixelModel* model = network.getPixelModel(0, r, c);
                if (model) {
                    std::cout << "  | Ry=" << std::fixed << std::setprecision(3) 
                             << model->getRy();
                } else {
                    std::cout << "  |";
                }
                if (c < cols - 1) {
                    std::cout << "      ";
                }
            }
            std::cout << "\n";
        }
    }
    std::cout << "\n";
}

// Standalone visualization function
int visualizePDN(int rows = 5, int cols = 5) {
    PDNNetwork network(rows, cols, 1);
    
    // Create sample pixel models
    float baseRx = 0.1f, baseRy = 0.1f, baseRz = 0.5f;
    std::vector<std::vector<std::unique_ptr<PixelModel>>> pixelModels(rows);
    
    for (int r = 0; r < rows; r++) {
        pixelModels[r].resize(cols);
        for (int c = 0; c < cols; c++) {
            pixelModels[r][c] = std::make_unique<PixelModel>(
                "model_" + std::to_string(r) + "_" + std::to_string(c),
                baseRx, baseRy, baseRz, 1.0f
            );
            network.setPixelModel(0, r, c, pixelModels[r][c].get());
        }
    }
    
    network.buildNetwork();
    
    // Print ASCII visualization
    printNetworkASCII(network);
    
    // Generate DOT file
    generateDOTFile(network, "pdn_network.dot");
    
    return 0;
}

