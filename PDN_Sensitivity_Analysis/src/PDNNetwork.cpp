#include "PDNNetwork.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>

PDNNetwork::PDNNetwork(int rows, int cols, int layers)
    : _rows(rows), _cols(cols), _layers(layers), _nextNodeId(1) {
    // Initialize pixel model grid
    _pixelModels.resize(_layers);
    for (int l = 0; l < _layers; l++) {
        _pixelModels[l].resize(_rows);
        for (int r = 0; r < _rows; r++) {
            _pixelModels[l][r].resize(_cols, nullptr);
        }
    }
}

void PDNNetwork::setPixelModel(int layer, int row, int col, PixelModel* model) {
    if (layer >= 0 && layer < _layers && 
        row >= 0 && row < _rows && 
        col >= 0 && col < _cols) {
        _pixelModels[layer][row][col] = model;
    }
}

PixelModel* PDNNetwork::getPixelModel(int layer, int row, int col) const {
    if (layer >= 0 && layer < _layers && 
        row >= 0 && row < _rows && 
        col >= 0 && col < _cols) {
        return _pixelModels[layer][row][col];
    }
    return nullptr;
}

int PDNNetwork::getNodeId(int layer, int row, int col, const std::string& pgNet) {
    auto key = std::make_tuple(layer, row, col, pgNet);
    if (_nodeMap.find(key) == _nodeMap.end()) {
        int nodeId = _nextNodeId++;
        _nodeMap[key] = nodeId;
        _nodes.push_back(Node(nodeId, 0.0, 0.0, pgNet));
    }
    return _nodeMap[key];
}

PDNNetwork::Node* PDNNetwork::getNode(int layer, int row, int col, const std::string& pgNet) {
    auto key = std::make_tuple(layer, row, col, pgNet);
    if (_nodeMap.find(key) != _nodeMap.end()) {
        int nodeId = _nodeMap[key];
        for (auto& node : _nodes) {
            if (node.id == nodeId) {
                return &node;
            }
        }
    }
    return nullptr;
}

void PDNNetwork::setCurrentLoad(int layer, int row, int col, double current, const std::string& pgNet) {
    Node* node = getNode(layer, row, col, pgNet);
    if (node) {
        node->currentLoad = current;
    } else {
        int nodeId = getNodeId(layer, row, col, pgNet);
        for (auto& n : _nodes) {
            if (n.id == nodeId) {
                n.currentLoad = current;
                break;
            }
        }
    }
}

void PDNNetwork::setVoltageSource(int layer, int row, int col, double voltage, const std::string& pgNet, bool isPad) {
    Node* node = getNode(layer, row, col, pgNet);
    if (node) {
        node->voltage = voltage;
        if (isPad) {
            node->isPad = true;
        }
    } else {
        int nodeId = getNodeId(layer, row, col, pgNet);
        for (auto& n : _nodes) {
            if (n.id == nodeId) {
                n.voltage = voltage;
                if (isPad) {
                    n.isPad = true;
                }
                break;
            }
        }
    }
}

void PDNNetwork::addBranch(int node1, int node2, double resistance, const std::string& direction) {
    _branches.push_back(Branch(node1, node2, resistance, direction));
}

void PDNNetwork::buildNetwork() {
    _nodes.clear();
    _branches.clear();
    _nodeMap.clear();
    _nextNodeId = 1;
    
    // Create nodes for all grid points
    for (int l = 0; l < _layers; l++) {
        for (int r = 0; r < _rows; r++) {
            for (int c = 0; c < _cols; c++) {
                getNodeId(l, r, c, "VDD");
                getNodeId(l, r, c, "GND");
            }
        }
    }
    
    // Build branches from pixel models
    for (int l = 0; l < _layers; l++) {
        for (int r = 0; r < _rows; r++) {
            for (int c = 0; c < _cols; c++) {
                PixelModel* model = _pixelModels[l][r][c];
                if (!model) continue;
                
                // Horizontal (x-direction) connections
                if (c < _cols - 1) {
                    int node1 = getNodeId(l, r, c, "VDD");
                    int node2 = getNodeId(l, r, c + 1, "VDD");
                    addBranch(node1, node2, model->getRx(), "x");
                }
                
                // Vertical (y-direction) connections
                if (r < _rows - 1) {
                    int node1 = getNodeId(l, r, c, "VDD");
                    int node2 = getNodeId(l, r + 1, c, "VDD");
                    addBranch(node1, node2, model->getRy(), "y");
                }
                
                // Vertical between layers (z-direction)
                if (l < _layers - 1) {
                    int node1 = getNodeId(l, r, c, "VDD");
                    int node2 = getNodeId(l + 1, r, c, "VDD");
                    addBranch(node1, node2, model->getRz(), "z");
                }
            }
        }
    }
}

void PDNNetwork::buildNetworkWithPerBranchResistances(
    const std::map<std::tuple<int, int, int, std::string>, double>& branchResistances) {
    _nodes.clear();
    _branches.clear();
    _nodeMap.clear();
    _nextNodeId = 1;
    
    // Create nodes for all grid points
    for (int l = 0; l < _layers; l++) {
        for (int r = 0; r < _rows; r++) {
            for (int c = 0; c < _cols; c++) {
                getNodeId(l, r, c, "VDD");
                getNodeId(l, r, c, "GND");
            }
        }
    }
    
    // Build branches with per-branch resistances
    for (int l = 0; l < _layers; l++) {
        for (int r = 0; r < _rows; r++) {
            for (int c = 0; c < _cols; c++) {
                // Horizontal (x-direction) connections
                if (c < _cols - 1) {
                    auto key = std::make_tuple(l, r, c, std::string("x"));
                    auto it = branchResistances.find(key);
                    double resistance = 0.1;  // Default if not found
                    if (it != branchResistances.end()) {
                        resistance = it->second;
                    }
                    
                    int node1 = getNodeId(l, r, c, "VDD");
                    int node2 = getNodeId(l, r, c + 1, "VDD");
                    addBranch(node1, node2, resistance, "x");
                }
                
                // Vertical (y-direction) connections
                if (r < _rows - 1) {
                    auto key = std::make_tuple(l, r, c, std::string("y"));
                    auto it = branchResistances.find(key);
                    double resistance = 0.1;  // Default if not found
                    if (it != branchResistances.end()) {
                        resistance = it->second;
                    }
                    
                    int node1 = getNodeId(l, r, c, "VDD");
                    int node2 = getNodeId(l, r + 1, c, "VDD");
                    addBranch(node1, node2, resistance, "y");
                }
                
                // Vertical between layers (z-direction)
                if (l < _layers - 1) {
                    auto key = std::make_tuple(l, r, c, std::string("z"));
                    auto it = branchResistances.find(key);
                    double resistance = 0.5;  // Default if not found
                    if (it != branchResistances.end()) {
                        resistance = it->second;
                    }
                    
                    int node1 = getNodeId(l, r, c, "VDD");
                    int node2 = getNodeId(l + 1, r, c, "VDD");
                    addBranch(node1, node2, resistance, "z");
                }
            }
        }
    }
}

double PDNNetwork::calculatePathResistance(const std::vector<int>& nodeIds) const {
    if (nodeIds.size() < 2) return 0.0;
    
    double totalResistance = 0.0;
    for (size_t i = 0; i < nodeIds.size() - 1; i++) {
        int node1 = nodeIds[i];
        int node2 = nodeIds[i + 1];
        
        // Find branch connecting these nodes
        for (const auto& branch : _branches) {
            if ((branch.node1Id == node1 && branch.node2Id == node2) ||
                (branch.node1Id == node2 && branch.node2Id == node1)) {
                totalResistance += branch.resistance;
                break;
            }
        }
    }
    return totalResistance;
}

void PDNNetwork::printNetworkInfo() const {
    std::cout << "=== PDN Network Information ===" << std::endl;
    std::cout << "Grid: " << _rows << " x " << _cols << " x " << _layers << std::endl;
    std::cout << "Nodes: " << _nodes.size() << std::endl;
    std::cout << "Branches: " << _branches.size() << std::endl;
    std::cout << std::endl;
    
    // Print resistance statistics
    std::map<std::string, std::vector<double>> resistancesByDir;
    for (const auto& branch : _branches) {
        resistancesByDir[branch.direction].push_back(branch.resistance);
    }
    
    for (const auto& [dir, res] : resistancesByDir) {
        if (!res.empty()) {
            double sum = 0.0, min = res[0], max = res[0];
            for (double r : res) {
                sum += r;
                if (r < min) min = r;
                if (r > max) max = r;
            }
            double avg = sum / res.size();
            std::cout << "Direction " << dir << ": " << res.size() << " branches" << std::endl;
            std::cout << "  Min: " << std::fixed << std::setprecision(6) << min << " Ω" << std::endl;
            std::cout << "  Max: " << std::fixed << std::setprecision(6) << max << " Ω" << std::endl;
            std::cout << "  Avg: " << std::fixed << std::setprecision(6) << avg << " Ω" << std::endl;
        }
    }
}

void PDNNetwork::exportBranchDataForVisualization(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    
    // Header
    // Note: row is exported in VoltSpot coordinate convention (bottom-origin):
    // row_out = (rows - 1) - row_internal
    file << "layer,row,col,direction,resistance,node1,node2" << std::endl;
    
    // Build reverse mapping: node ID -> grid coordinates
    std::map<int, std::tuple<int, int, int>> nodeToGrid;
    for (const auto& [key, nodeId] : _nodeMap) {
        int layer = std::get<0>(key);
        int row = std::get<1>(key);
        int col = std::get<2>(key);
        std::string pgNet = std::get<3>(key);
        
        if (pgNet == "VDD") {
            nodeToGrid[nodeId] = std::make_tuple(layer, row, col);
        }
    }
    
    // Export branches
    for (const auto& branch : _branches) {
        auto it1 = nodeToGrid.find(branch.node1Id);
        auto it2 = nodeToGrid.find(branch.node2Id);
        
        if (it1 != nodeToGrid.end() && it2 != nodeToGrid.end()) {
            int layer1 = std::get<0>(it1->second);
            int row1 = std::get<1>(it1->second);
            int col1 = std::get<2>(it1->second);
            
            // Use the "from" node coordinates for branch location
            int branchLayer = layer1;
            int branchRow = (_rows - 1) - row1; // VoltSpot bottom-origin
            int branchCol = col1;
            
            file << branchLayer << "," << branchRow << "," << branchCol << ","
                 << branch.direction << "," << std::fixed << std::setprecision(8)
                 << branch.resistance << "," << branch.node1Id << "," << branch.node2Id << std::endl;
        }
    }
    
    file.close();
    std::cout << "Branch data exported to " << filename << std::endl;
}

void PDNNetwork::exportNodeDataForVisualization(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    
    // Header
    // Note: row is exported in VoltSpot coordinate convention (bottom-origin):
    // row_out = (rows - 1) - row_internal
    file << "layer,row,col,pgNet,nodeId,voltage,currentLoad,isPad" << std::endl;
    
    // Export all nodes with their coordinates
    for (const auto& [key, nodeId] : _nodeMap) {
        int layer = std::get<0>(key);
        int row = std::get<1>(key);
        int col = std::get<2>(key);
        std::string pgNet = std::get<3>(key);

        int row_out = (_rows - 1) - row; // VoltSpot bottom-origin
        
        // Find the node
        for (const auto& node : _nodes) {
            if (node.id == nodeId) {
                file << layer << "," << row_out << "," << col << "," << pgNet << ","
                     << nodeId << "," << std::fixed << std::setprecision(8)
                     << node.voltage << "," << node.currentLoad << "," << (node.isPad ? 1 : 0) << std::endl;
                break;
            }
        }
    }
    
    file.close();
    std::cout << "Node data exported to " << filename << std::endl;
}

