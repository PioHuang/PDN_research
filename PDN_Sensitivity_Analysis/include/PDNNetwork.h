#ifndef PDN_NETWORK_H
#define PDN_NETWORK_H

#include "PixelModel.h"
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <tuple>

/**
 * Simple PDN Network structure for sensitivity analysis
 * Represents a grid-based power delivery network
 */
class PDNNetwork {
public:
    struct Node {
        int id;
        double voltage;
        double currentLoad;
        std::string pgNetName;  // Power/Ground net name (e.g., "VDD", "GND")
        
        Node(int id = 0, double voltage = 0.0, double currentLoad = 0.0, const std::string& pgNet = "")
            : id(id), voltage(voltage), currentLoad(currentLoad), pgNetName(pgNet) {}
    };

    struct Branch {
        int node1Id;
        int node2Id;
        double resistance;
        std::string direction;  // "x", "y", or "z"
        
        Branch(int n1, int n2, double r, const std::string& dir = "")
            : node1Id(n1), node2Id(n2), resistance(r), direction(dir) {}
    };

    PDNNetwork(int rows, int cols, int layers = 1);
    
    // Grid setup
    void setPixelModel(int layer, int row, int col, PixelModel* model);
    PixelModel* getPixelModel(int layer, int row, int col) const;
    
    // Node operations
    Node* getNode(int layer, int row, int col, const std::string& pgNet = "VDD");
    void setCurrentLoad(int layer, int row, int col, double current, const std::string& pgNet = "VDD");
    void setVoltageSource(int layer, int row, int col, double voltage, const std::string& pgNet = "VDD");
    
    // Build network from pixel models
    void buildNetwork();
    
    // Get network statistics
    int getNodeCount() const { return _nodes.size(); }
    int getBranchCount() const { return _branches.size(); }
    std::vector<Branch> getBranches() const { return _branches; }
    std::vector<Node> getNodes() const { return _nodes; }
    int getLayerCount() const { return _layers; }
    int getLayerRows(int layer) const { 
        if (layer >= 0 && layer < _layers) return _rows; 
        return 0; 
    }
    int getLayerCols(int layer) const { 
        if (layer >= 0 && layer < _layers) return _cols; 
        return 0; 
    }
    
    // Calculate equivalent resistance for a path
    double calculatePathResistance(const std::vector<int>& nodeIds) const;
    
    // Print network info
    void printNetworkInfo() const;

private:
    int _rows, _cols, _layers;
    std::vector<Node> _nodes;
    std::vector<Branch> _branches;
    
    // Grid structure: [layer][row][col] -> PixelModel*
    std::vector<std::vector<std::vector<PixelModel*>>> _pixelModels;
    
    // Node mapping: (layer, row, col, pgNet) -> node index
    std::map<std::tuple<int, int, int, std::string>, int> _nodeMap;
    
    int _nextNodeId;
    
    int getNodeId(int layer, int row, int col, const std::string& pgNet);
    void addBranch(int node1, int node2, double resistance, const std::string& direction);
};

#endif // PDN_NETWORK_H

