#include "PDNNetwork.h"
#include "PixelModel.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>

/**
 * Sensitivity Analysis Tool for Pixel Model R Values
 * 
 * This tool sweeps R values (Rx, Ry, Rz) and analyzes their impact on:
 * - Total network resistance
 * - Path resistance
 * - IR drop (if voltage/current sources are set)
 */

struct SensitivityResult {
    float rx;
    float ry;
    float rz;
    double totalResistance;
    double avgBranchResistance;
    double maxBranchResistance;
    double minBranchResistance;
    double pathResistance;  // Resistance along a specific path
};

class SensitivityAnalyzer {
public:
    SensitivityAnalyzer(PDNNetwork* network) : _network(network) {}
    
    // Sweep a single resistance value (Rx, Ry, or Rz) while keeping others constant
    std::vector<SensitivityResult> sweepResistance(
        const std::string& resistanceType,  // "Rx", "Ry", or "Rz"
        float baseRx, float baseRy, float baseRz,
        float startRatio, float endRatio, int steps) {
        
        std::vector<SensitivityResult> results;
        float stepSize = (endRatio - startRatio) / (steps - 1);
        
        // Store original pixel models
        std::vector<std::vector<std::vector<PixelModel*>>> originalModels;
        savePixelModels(originalModels);
        
        for (int i = 0; i < steps; i++) {
            float ratio = startRatio + i * stepSize;
            
            // Create modified pixel models
            float rx = baseRx;
            float ry = baseRy;
            float rz = baseRz;
            
            if (resistanceType == "Rx") rx *= ratio;
            else if (resistanceType == "Ry") ry *= ratio;
            else if (resistanceType == "Rz") rz *= ratio;
            
            // Apply to all pixel models
            applyResistancesToNetwork(rx, ry, rz);
            
            // Rebuild network
            _network->buildNetwork();
            
            // Calculate metrics
            SensitivityResult result;
            result.rx = rx;
            result.ry = ry;
            result.rz = rz;
            calculateMetrics(result);
            
            results.push_back(result);
        }
        
        // Restore original models
        restorePixelModels(originalModels);
        _network->buildNetwork();
        
        return results;
    }
    
    // Sweep all resistances together (same ratio for all)
    std::vector<SensitivityResult> sweepAllResistances(
        float baseRx, float baseRy, float baseRz,
        float startRatio, float endRatio, int steps) {
        
        std::vector<SensitivityResult> results;
        float stepSize = (endRatio - startRatio) / (steps - 1);
        
        std::vector<std::vector<std::vector<PixelModel*>>> originalModels;
        savePixelModels(originalModels);
        
        for (int i = 0; i < steps; i++) {
            float ratio = startRatio + i * stepSize;
            
            float rx = baseRx * ratio;
            float ry = baseRy * ratio;
            float rz = baseRz * ratio;
            
            applyResistancesToNetwork(rx, ry, rz);
            _network->buildNetwork();
            
            SensitivityResult result;
            result.rx = rx;
            result.ry = ry;
            result.rz = rz;
            calculateMetrics(result);
            
            results.push_back(result);
        }
        
        restorePixelModels(originalModels);
        _network->buildNetwork();
        
        return results;
    }
    
    // Write results to CSV
    void writeResultsToCSV(const std::vector<SensitivityResult>& results, 
                          const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return;
        }
        
        // Header
        file << "Rx,Ry,Rz,TotalResistance,AvgBranchResistance,MaxBranchResistance,MinBranchResistance,PathResistance" << std::endl;
        
        // Data
        for (const auto& result : results) {
            file << std::fixed << std::setprecision(6)
                 << result.rx << ","
                 << result.ry << ","
                 << result.rz << ","
                 << result.totalResistance << ","
                 << result.avgBranchResistance << ","
                 << result.maxBranchResistance << ","
                 << result.minBranchResistance << ","
                 << result.pathResistance << std::endl;
        }
        
        file.close();
        std::cout << "Results written to " << filename << std::endl;
    }
    
    // Print sensitivity summary
    void printSensitivitySummary(const std::vector<SensitivityResult>& results,
                                const std::string& resistanceType) {
        if (results.empty()) return;
        
        std::cout << "\n=== Sensitivity Analysis: " << resistanceType << " ===" << std::endl;
        std::cout << "Steps: " << results.size() << std::endl;
        
        // Calculate sensitivity (change in metric / change in resistance)
        if (results.size() >= 2) {
            double resistanceChange = results.back().rx - results.front().rx;
            if (resistanceType == "Ry") {
                resistanceChange = results.back().ry - results.front().ry;
            } else if (resistanceType == "Rz") {
                resistanceChange = results.back().rz - results.front().rz;
            }
            
            double totalResChange = results.back().totalResistance - results.front().totalResistance;
            double sensitivity = (resistanceChange != 0) ? totalResChange / resistanceChange : 0.0;
            
            std::cout << "Resistance Range: [" 
                      << std::fixed << std::setprecision(6)
                      << results.front().rx << ", " << results.back().rx << "] Ω" << std::endl;
            std::cout << "Total Resistance Range: ["
                      << results.front().totalResistance << ", " 
                      << results.back().totalResistance << "] Ω" << std::endl;
            std::cout << "Sensitivity: " << std::fixed << std::setprecision(6) 
                      << sensitivity << " (ΔTotalR / Δ" << resistanceType << ")" << std::endl;
        }
        
        std::cout << "\nFirst result:" << std::endl;
        printResult(results.front());
        std::cout << "\nLast result:" << std::endl;
        printResult(results.back());
    }
    
private:
    PDNNetwork* _network;
    std::vector<std::vector<std::vector<PixelModel*>>> _savedModels;
    
    void savePixelModels(std::vector<std::vector<std::vector<PixelModel*>>>& /* saved */) {
        // Note: In a real implementation, you'd save the actual resistance values
        // For this simplified version, we'll just note that we need to restore
        _savedModels.clear();
    }
    
    void restorePixelModels(std::vector<std::vector<std::vector<PixelModel*>>>& /* saved */) {
        // Restore original resistance values
        // In this simplified version, we'll rebuild the network
    }
    
    void applyResistancesToNetwork(float rx, float ry, float rz) {
        // Apply resistances to all pixel models in the network
        // This is a simplified version that assumes uniform resistances
        // In reality, you'd iterate through all models and update their values
        for (int l = 0; l < _network->getLayerCount(); l++) {
            for (int r = 0; r < _network->getLayerRows(l); r++) {
                for (int c = 0; c < _network->getLayerCols(l); c++) {
                    PixelModel* model = _network->getPixelModel(l, r, c);
                    if (model) {
                        model->setResistances(rx, ry, rz);
                    }
                }
            }
        }
    }
    
    void calculateMetrics(SensitivityResult& result) {
        auto branches = _network->getBranches();
        
        if (branches.empty()) {
            result.totalResistance = 0.0;
            result.avgBranchResistance = 0.0;
            result.maxBranchResistance = 0.0;
            result.minBranchResistance = 0.0;
            return;
        }
        
        double total = 0.0;
        double min = branches[0].resistance;
        double max = branches[0].resistance;
        
        for (const auto& branch : branches) {
            total += branch.resistance;
            if (branch.resistance < min) min = branch.resistance;
            if (branch.resistance > max) max = branch.resistance;
        }
        
        result.totalResistance = total;
        result.avgBranchResistance = total / branches.size();
        result.maxBranchResistance = max;
        result.minBranchResistance = min;
        
        // Calculate path resistance (example: from corner to corner)
        // This is a simplified example - in reality you'd find the actual path
        // For now, we'll use a simple approximation: average resistance * path length
        int rows = _network->getLayerRows(0);
        int cols = _network->getLayerCols(0);
        int pathLength = rows + cols - 2;  // Manhattan distance
        result.pathResistance = result.avgBranchResistance * pathLength;
    }
    
    void printResult(const SensitivityResult& result) {
        std::cout << "  Rx: " << std::fixed << std::setprecision(6) << result.rx << " Ω" << std::endl;
        std::cout << "  Ry: " << result.ry << " Ω" << std::endl;
        std::cout << "  Rz: " << result.rz << " Ω" << std::endl;
        std::cout << "  Total Resistance: " << result.totalResistance << " Ω" << std::endl;
        std::cout << "  Avg Branch Resistance: " << result.avgBranchResistance << " Ω" << std::endl;
        std::cout << "  Min Branch Resistance: " << result.minBranchResistance << " Ω" << std::endl;
        std::cout << "  Max Branch Resistance: " << result.maxBranchResistance << " Ω" << std::endl;
    }
};

// Example usage
int main(int /* argc */, char* /* argv */[]) {
    // Create a simple test network
    int rows = 10, cols = 10, layers = 1;
    PDNNetwork network(rows, cols, layers);
    
    // Create pixel models with base resistance values
    float baseRx = 0.1f;  // 0.1 Ω
    float baseRy = 0.1f;  // 0.1 Ω
    float baseRz = 0.5f;  // 0.5 Ω
    
    // Store pixel models to manage lifetime
    std::vector<std::vector<std::unique_ptr<PixelModel>>> pixelModels(rows);
    
    // Initialize all grid points with pixel models
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
    network.printNetworkInfo();
    
    // Run sensitivity analysis
    SensitivityAnalyzer analyzer(&network);
    
    std::cout << "\n=== Running Sensitivity Analysis ===" << std::endl;
    
    // Sweep Rx from 0.5x to 2.0x
    auto rxResults = analyzer.sweepResistance("Rx", baseRx, baseRy, baseRz, 0.5f, 2.0f, 20);
    analyzer.printSensitivitySummary(rxResults, "Rx");
    analyzer.writeResultsToCSV(rxResults, "sensitivity_rx.csv");
    
    // Sweep Ry from 0.5x to 2.0x
    auto ryResults = analyzer.sweepResistance("Ry", baseRx, baseRy, baseRz, 0.5f, 2.0f, 20);
    analyzer.printSensitivitySummary(ryResults, "Ry");
    analyzer.writeResultsToCSV(ryResults, "sensitivity_ry.csv");
    
    // Sweep Rz from 0.5x to 2.0x
    auto rzResults = analyzer.sweepResistance("Rz", baseRx, baseRy, baseRz, 0.5f, 2.0f, 20);
    analyzer.printSensitivitySummary(rzResults, "Rz");
    analyzer.writeResultsToCSV(rzResults, "sensitivity_rz.csv");
    
    // Sweep all resistances together
    auto allResults = analyzer.sweepAllResistances(baseRx, baseRy, baseRz, 0.5f, 2.0f, 20);
    analyzer.printSensitivitySummary(allResults, "All");
    analyzer.writeResultsToCSV(allResults, "sensitivity_all.csv");
    
    return 0;
}

