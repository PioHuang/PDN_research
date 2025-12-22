#include "PDNLoadMapper.h"
#include "PDNNetwork.h"
#include "VoltSpotVirtualGrid.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <map>

namespace
{

    std::string trim(const std::string &s)
    {
        size_t start = 0;
        while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])))
        {
            start++;
        }
        size_t end = s.size();
        while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1])))
        {
            end--;
        }
        return s.substr(start, end - start);
    }

    std::string stripComment(const std::string &line)
    {
        const size_t pos = line.find('#');
        if (pos == std::string::npos)
        {
            return line;
        }
        return line.substr(0, pos);
    }

    std::vector<std::string> splitWS(const std::string &s)
    {
        std::istringstream iss(s);
        std::vector<std::string> out;
        std::string tok;
        while (iss >> tok)
        {
            out.push_back(tok);
        }
        return out;
    }

    double parseDouble(const std::string &s)
    {
        try
        {
            size_t idx = 0;
            const double v = std::stod(s, &idx);
            if (idx != s.size())
            {
                throw std::invalid_argument("trailing characters");
            }
            return v;
        }
        catch (const std::exception &)
        {
            throw std::invalid_argument("Invalid double: " + s);
        }
    }

    int parseInt(const std::string &s)
    {
        try
        {
            size_t idx = 0;
            const int v = std::stoi(s, &idx);
            if (idx != s.size())
            {
                throw std::invalid_argument("trailing characters");
            }
            return v;
        }
        catch (const std::exception &)
        {
            throw std::invalid_argument("Invalid int: " + s);
        }
    }

    struct FlpUnit
    {
        std::string name;
        double width_m;
        double height_m;
        double leftX_m;
        double bottomY_m;
    };

    std::vector<FlpUnit> parseFlpFile(const std::string &flpPath)
    {
        std::ifstream in(flpPath);
        if (!in.is_open())
        {
            throw std::runtime_error("Cannot open flp file: " + flpPath);
        }

        std::vector<FlpUnit> units;
        std::string line;
        while (std::getline(in, line))
        {
            line = trim(stripComment(line));
            if (line.empty())
                continue;

            const auto toks = splitWS(line);
            if (toks.size() < 5)
            {
                throw std::runtime_error("Invalid flp line (need 5 columns): " + line);
            }

            FlpUnit unit;
            unit.name = toks[0];
            unit.width_m = parseDouble(toks[1]);
            unit.height_m = parseDouble(toks[2]);
            unit.leftX_m = parseDouble(toks[3]);
            unit.bottomY_m = parseDouble(toks[4]);
            units.push_back(unit);
        }

        return units;
    }

    std::map<std::string, double> parsePtraceFile(const std::string &ptracePath)
    {
        std::ifstream in(ptracePath);
        if (!in.is_open())
        {
            throw std::runtime_error("Cannot open ptrace file: " + ptracePath);
        }

        std::string line;

        // Read unit names from first line
        std::getline(in, line);
        line = trim(line);
        if (line.empty())
        {
            throw std::runtime_error("Empty ptrace file (no unit names)");
        }

        const auto unitNames = splitWS(line);

        // Read first power trace line (or average if multiple lines)
        std::vector<double> powerSum(unitNames.size(), 0.0);
        int traceCount = 0;

        while (std::getline(in, line))
        {
            line = trim(line);
            if (line.empty())
                continue;

            const auto powerValues = splitWS(line);
            if (powerValues.size() != unitNames.size())
            {
                throw std::runtime_error("Power trace line has " + std::to_string(powerValues.size()) +
                                         " values, expected " + std::to_string(unitNames.size()));
            }

            for (size_t i = 0; i < powerValues.size(); i++)
            {
                powerSum[i] += parseDouble(powerValues[i]);
            }
            traceCount++;
        }

        // Average the power values
        std::map<std::string, double> unitPower;
        for (size_t i = 0; i < unitNames.size(); i++)
        {
            double avgPower = traceCount > 0 ? powerSum[i] / traceCount : powerSum[i];
            unitPower[unitNames[i]] = avgPower;
        }

        return unitPower;
    }

    int clamp(int value, int minVal, int maxVal)
    {
        return std::max(minVal, std::min(maxVal, value));
    }

    int physicalToGridX(double physicalX_m, double chipWidth_m, int gridCols)
    {
        if (gridCols <= 1)
            return 0;
        double normalized = physicalX_m / chipWidth_m;
        int gridX = static_cast<int>(std::round(normalized * (gridCols - 1)));
        return clamp(gridX, 0, gridCols - 1);
    }

    int physicalToGridY(double physicalY_m, double chipHeight_m, int gridRows)
    {
        if (gridRows <= 1)
            return 0;
        double normalized = physicalY_m / chipHeight_m;
        int gridY = static_cast<int>(std::round(normalized * (gridRows - 1)));
        return clamp(gridY, 0, gridRows - 1);
    }

} // namespace

void PDNLoadMapper::loadVoltageSources(PDNNetwork &network, const std::string &padlocPath,
                                       const VoltSpotVirtualGrid::Result &gridResult,
                                       double vddVoltage, double gndVoltage, int padlocFormat)
{
    std::ifstream in(padlocPath);
    if (!in.is_open())
    {
        throw std::runtime_error("Cannot open padloc file: " + padlocPath);
    }

    const int gridRows = gridResult.grid.rows;
    const int gridCols = gridResult.grid.cols;
    const double chipWidth_m = gridResult.chip.width_m;
    const double chipHeight_m = gridResult.chip.height_m;
    const double padPitch_m = gridResult.padPitch_m;
    // const int gridIntv = gridResult.gridIntv; // (unused)

    int padCount = 0;
    int invalidCount = 0;

    std::string line;
    while (std::getline(in, line))
    {
        line = trim(stripComment(line));
        if (line.empty())
            continue;

        const auto toks = splitWS(line);
        if (toks.size() < 3)
        {
            continue;
        }

        std::string padType = toks[0];
        int x, y;

        if (padlocFormat == 0 || padlocFormat == 1)
        {
            // Virtual grid coordinates (direct)
            x = parseInt(toks[1]);
            y = parseInt(toks[2]);
        }
        else
        {
            // Pad grid coordinates (need conversion)
            int padX = parseInt(toks[1]);
            int padY = parseInt(toks[2]);

            const int padCols = static_cast<int>(std::floor(chipWidth_m / padPitch_m));
            const int padRows = static_cast<int>(std::floor(chipHeight_m / padPitch_m));
            const int itvCol = (gridCols - 1) / (padCols - 1);
            const int itvRow = (gridRows - 1) / (padRows - 1);

            x = padX * itvCol;
            y = (padRows - padY - 1) * itvRow;
        }

        // Convert y coordinate (VoltSpot uses inverted y-axis)
        int gridY = gridRows - y - 1;
        int gridX = x;

        // Boundary check
        if (gridX < 0 || gridX >= gridCols || gridY < 0 || gridY >= gridRows)
        {
            std::cerr << "Warning: Pad location out of range: (" << x << ", " << y
                      << ") -> grid(" << gridX << ", " << gridY << ")" << std::endl;
            invalidCount++;
            continue;
        }

        // Set voltage source
        if (padType == "V" || padType == "v")
        {
            network.setVoltageSource(0, gridY, gridX, vddVoltage, "VDD", true);
            padCount++;
        }
        else if (padType == "G" || padType == "g")
        {
            network.setVoltageSource(0, gridY, gridX, gndVoltage, "GND", true);
            padCount++;
        }
        else
        {
            std::cerr << "Warning: Unknown pad type: " << padType << std::endl;
        }
    }

    std::cout << "Loaded " << padCount << " voltage sources from " << padlocPath << std::endl;
    if (invalidCount > 0)
    {
        std::cerr << "Warning: " << invalidCount << " invalid pad locations skipped" << std::endl;
    }
}

void PDNLoadMapper::loadCurrentLoads(PDNNetwork &network, const std::string &flpPath,
                                     const std::string &ptracePath,
                                     const VoltSpotVirtualGrid::Result &gridResult,
                                     double supplyVoltage)
{
    // Parse FLP file
    const auto flpUnits = parseFlpFile(flpPath);

    // Parse PTrace file
    const auto unitPower = parsePtraceFile(ptracePath);

    const int gridRows = gridResult.grid.rows;
    const int gridCols = gridResult.grid.cols;
    const double chipWidth_m = gridResult.chip.width_m;
    const double chipHeight_m = gridResult.chip.height_m;
    const double dx_m = gridResult.grid.dx_m;
    const double dy_m = gridResult.grid.dy_m;

    // Initialize grid current load map
    std::map<std::pair<int, int>, double> gridCurrentLoads;

    // Map each FLP unit to grid cells
    for (const auto &unit : flpUnits)
    {
        auto it = unitPower.find(unit.name);
        if (it == unitPower.end())
        {
            std::cerr << "Warning: Unit " << unit.name << " not found in ptrace file" << std::endl;
            continue;
        }

        double powerW = it->second;
        if (powerW <= 0.0)
            continue;

        // Calculate unit boundaries in grid coordinates
        double rightX_m = unit.leftX_m + unit.width_m;
        double topY_m = unit.bottomY_m + unit.height_m;

        int leftGridX = physicalToGridX(unit.leftX_m, chipWidth_m, gridCols);
        int rightGridX = physicalToGridX(rightX_m, chipWidth_m, gridCols);
        int bottomGridY = physicalToGridY(unit.bottomY_m, chipHeight_m, gridRows);
        int topGridY = physicalToGridY(topY_m, chipHeight_m, gridRows);

        // Ensure valid range
        leftGridX = clamp(leftGridX, 0, gridCols - 1);
        rightGridX = clamp(rightGridX, 0, gridCols - 1);
        bottomGridY = clamp(bottomGridY, 0, gridRows - 1);
        topGridY = clamp(topGridY, 0, gridRows - 1);

        // Calculate unit area
        double unitArea_m2 = unit.width_m * unit.height_m;
        if (unitArea_m2 <= 0.0)
            continue;

        // Distribute power to grid cells that overlap with this unit
        for (int gridY = bottomGridY; gridY <= topGridY; gridY++)
        {
            for (int gridX = leftGridX; gridX <= rightGridX; gridX++)
            {
                // Calculate grid cell boundaries in physical coordinates
                double cellLeftX_m = gridX * dx_m;
                double cellRightX_m = (gridX + 1) * dx_m;
                double cellBottomY_m = gridY * dy_m;
                double cellTopY_m = (gridY + 1) * dy_m;

                // Calculate overlap area
                double overlapLeftX = std::max(unit.leftX_m, cellLeftX_m);
                double overlapRightX = std::min(rightX_m, cellRightX_m);
                double overlapBottomY = std::max(unit.bottomY_m, cellBottomY_m);
                double overlapTopY = std::min(topY_m, cellTopY_m);

                if (overlapRightX > overlapLeftX && overlapTopY > overlapBottomY)
                {
                    double overlapArea_m2 = (overlapRightX - overlapLeftX) * (overlapTopY - overlapBottomY);
                    double powerFraction = overlapArea_m2 / unitArea_m2;
                    double cellPowerW = powerW * powerFraction;
                    double cellCurrentA = cellPowerW / supplyVoltage;

                    // VoltSpot internal grid indexing uses row 0 at TOP, while physical Y=0 is at BOTTOM.
                    // Convert physical-grid row (bottom-origin) to internal row (top-origin).
                    int internalRow = (gridRows - 1) - gridY;
                    auto key = std::make_pair(internalRow, gridX);
                    gridCurrentLoads[key] += cellCurrentA;
                }
            }
        }
    }

    // Apply current loads to network
    int loadCount = 0;
    for (const auto &[key, current] : gridCurrentLoads)
    {
        int gridY = key.first;
        int gridX = key.second;
        // VoltSpot models load current between VDD and GND rails:
        // VDD gets -I in RHS, GND gets +I in RHS. Store the same magnitude on both.
        network.setCurrentLoad(0, gridY, gridX, current, "VDD");
        network.setCurrentLoad(0, gridY, gridX, current, "GND");
        loadCount++;
    }

    std::cout << "Loaded current loads from " << flpPath << " and " << ptracePath << std::endl;
    std::cout << "  Total units: " << flpUnits.size() << std::endl;
    std::cout << "  Grid cells with load: " << loadCount << std::endl;
}
