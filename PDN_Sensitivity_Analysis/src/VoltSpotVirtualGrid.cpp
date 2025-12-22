#include "VoltSpotVirtualGrid.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

std::string trim(const std::string& s) {
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
        start++;
    }
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
        end--;
    }
    return s.substr(start, end - start);
}

std::string stripComment(const std::string& line) {
    const size_t pos = line.find('#');
    if (pos == std::string::npos) {
        return line;
    }
    return line.substr(0, pos);
}

std::vector<std::string> splitWS(const std::string& s) {
    std::istringstream iss(s);
    std::vector<std::string> out;
    std::string tok;
    while (iss >> tok) {
        out.push_back(tok);
    }
    return out;
}

double parseDouble(const std::string& s) {
    try {
        size_t idx = 0;
        const double v = std::stod(s, &idx);
        if (idx != s.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return v;
    } catch (const std::exception&) {
        throw std::invalid_argument("Invalid double: " + s);
    }
}

int parseInt(const std::string& s) {
    try {
        size_t idx = 0;
        const int v = std::stoi(s, &idx);
        if (idx != s.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return v;
    } catch (const std::exception&) {
        throw std::invalid_argument("Invalid int: " + s);
    }
}

bool isNullToken(const std::string& s) {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return (t == "(null)" || t == "null");
}

VoltSpotVirtualGrid::Direction parseDirection01(const std::string& s) {
    const int d = parseInt(s);
    if (d == 0) return VoltSpotVirtualGrid::Direction::X;
    if (d == 1) return VoltSpotVirtualGrid::Direction::Y;
    throw std::invalid_argument("Direction must be 0 (X) or 1 (Y), got: " + s);
}

} // namespace

VoltSpotVirtualGrid::ChipSize VoltSpotVirtualGrid::parseVoltSpotFlpChipSize(const std::string& flpPath) {
    std::ifstream in(flpPath);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open flp file: " + flpPath);
    }

    double xmax = 0.0;
    double ymax = 0.0;

    std::string line;
    while (std::getline(in, line)) {
        line = trim(stripComment(line));
        if (line.empty()) continue;

        // VoltSpot flp format:
        // <name> <width> <height> <left_x> <bottom_y>   (meters)
        const auto toks = splitWS(line);
        if (toks.size() < 5) {
            throw std::runtime_error("Invalid flp line (need 5 columns): " + line);
        }

        const double w = parseDouble(toks[1]);
        const double h = parseDouble(toks[2]);
        const double x = parseDouble(toks[3]);
        const double y = parseDouble(toks[4]);

        xmax = std::max(xmax, x + w);
        ymax = std::max(ymax, y + h);
    }

    if (xmax <= 0.0 || ymax <= 0.0) {
        throw std::runtime_error("Computed non-positive chip size from flp: " + flpPath);
    }

    return ChipSize{.width_m = xmax, .height_m = ymax};
}

std::unordered_map<std::string, std::string> VoltSpotVirtualGrid::parseVoltSpotConfigFile(const std::string& configPath) {
    std::ifstream in(configPath);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config file: " + configPath);
    }

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(stripComment(line));
        if (line.empty()) continue;

        const auto toks = splitWS(line);
        if (toks.empty()) continue;

        // VoltSpot config format uses "-key value".
        if (toks[0].size() < 2 || toks[0][0] != '-') continue;
        const std::string key = toks[0].substr(1);
        if (toks.size() < 2) continue;
        kv[key] = toks[1];
    }
    return kv;
}

VoltSpotVirtualGrid::Config VoltSpotVirtualGrid::configFromMap(const std::unordered_map<std::string, std::string>& kv) {
    Config cfg;

    const auto itPitch = kv.find("PDN_padpitch");
    if (itPitch != kv.end()) {
        cfg.padPitch_m = parseDouble(itPitch->second);
    }

    const auto itIntv = kv.find("PDN_grid_intv");
    if (itIntv != kv.end()) {
        cfg.gridIntv = parseInt(itIntv->second);
    }

    const auto itMlcf = kv.find("mlayer_spec_file");
    if (itMlcf != kv.end() && !isNullToken(itMlcf->second)) {
        cfg.mlcfPath = itMlcf->second;
    }

    return cfg;
}

std::vector<VoltSpotVirtualGrid::MetalLayerSpec> VoltSpotVirtualGrid::parseMlcf(const std::string& mlcfPath) {
    std::ifstream in(mlcfPath);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open mlcf file: " + mlcfPath);
    }

    // mlcf file format (6 lines per layer):
    // layerId, pitch, width, thick, rho, direction(0/1)
    std::vector<std::string> tokens;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(stripComment(line));
        if (line.empty()) continue;
        const auto toks = splitWS(line);
        tokens.insert(tokens.end(), toks.begin(), toks.end());
    }

    if (tokens.size() % 6 != 0) {
        throw std::runtime_error("Invalid mlcf token count; expected multiple of 6: " + mlcfPath);
    }

    const int n = static_cast<int>(tokens.size() / 6);
    if (n <= 0) {
        throw std::runtime_error("Empty mlcf: " + mlcfPath);
    }

    std::vector<MetalLayerSpec> layers;
    layers.resize(static_cast<size_t>(n));
    std::vector<bool> seen(static_cast<size_t>(n), false);

    for (int k = 0; k < n; k++) {
        const int base = k * 6;
        const int layerId = parseInt(tokens[base + 0]);
        if (layerId < 0 || layerId >= n) {
            throw std::runtime_error("Invalid mlcf layerId: " + std::to_string(layerId));
        }
        if (seen[static_cast<size_t>(layerId)]) {
            throw std::runtime_error("Duplicate mlcf layerId: " + std::to_string(layerId));
        }
        seen[static_cast<size_t>(layerId)] = true;

        MetalLayerSpec spec;
        spec.layerId = layerId;
        spec.pitch_m = parseDouble(tokens[base + 1]);
        spec.width_m = parseDouble(tokens[base + 2]);
        spec.thickness_m = parseDouble(tokens[base + 3]);
        spec.rho_ohm_m = parseDouble(tokens[base + 4]);
        spec.direction = parseDirection01(tokens[base + 5]);

        layers[static_cast<size_t>(layerId)] = spec;
    }

    // VoltSpot enforces: even number of layers, layer0 X, alternating direction.
    if (n % 2 != 0) {
        throw std::runtime_error("VoltSpot requires even number of metal layers in mlcf.");
    }
    for (int i = 0; i < n; i++) {
        const Direction expected = (i % 2 == 0) ? Direction::X : Direction::Y;
        if (layers[static_cast<size_t>(i)].direction != expected) {
            throw std::runtime_error("VoltSpot expects metal layer directions: 0=X,1=Y,2=X,...; mismatch at layer " +
                                     std::to_string(i));
        }
    }

    return layers;
}

std::vector<VoltSpotVirtualGrid::MetalLayerSpec> VoltSpotVirtualGrid::defaultMetalStack() {
    // Matches VoltSpot populate_default_mlayers() (PDN_sim.c).
    std::vector<MetalLayerSpec> layers;
    layers.push_back(MetalLayerSpec{.layerId = 0, .pitch_m = 30e-6, .width_m = 10e-6, .thickness_m = 3.5e-6, .rho_ohm_m = 1.68e-8, .direction = Direction::X});
    layers.push_back(MetalLayerSpec{.layerId = 1, .pitch_m = 30e-6, .width_m = 10e-6, .thickness_m = 3.5e-6, .rho_ohm_m = 1.68e-8, .direction = Direction::Y});
    layers.push_back(MetalLayerSpec{.layerId = 2, .pitch_m = 810e-9, .width_m = 400e-9, .thickness_m = 720e-9, .rho_ohm_m = 1.68e-8, .direction = Direction::X});
    layers.push_back(MetalLayerSpec{.layerId = 3, .pitch_m = 810e-9, .width_m = 400e-9, .thickness_m = 720e-9, .rho_ohm_m = 1.68e-8, .direction = Direction::Y});
    return layers;
}

VoltSpotVirtualGrid::Result VoltSpotVirtualGrid::build(const ChipSize& chip, const Config& cfg) {
    if (chip.width_m <= 0.0 || chip.height_m <= 0.0) {
        throw std::invalid_argument("Chip size must be positive");
    }
    if (cfg.padPitch_m <= 0.0) {
        throw std::invalid_argument("PDN_padpitch must be > 0");
    }
    if (cfg.gridIntv < 1) {
        throw std::invalid_argument("PDN_grid_intv must be >= 1");
    }

    Result out;
    out.chip = chip;
    out.padPitch_m = cfg.padPitch_m;
    out.gridIntv = cfg.gridIntv;

    // VoltSpot sizing:
    // pad_grid_col = floor(width / padPitch)
    // pad_grid_row = floor(height / padPitch)
    // cols = gridIntv * (pad_grid_col - 1) + 1
    // rows = gridIntv * (pad_grid_row - 1) + 1
    const int padCols = static_cast<int>(std::floor(chip.width_m / cfg.padPitch_m));
    const int padRows = static_cast<int>(std::floor(chip.height_m / cfg.padPitch_m));
    if (padCols < 1 || padRows < 1) {
        throw std::runtime_error("PDN_padpitch larger than chip size (pad grid would be < 1).");
    }

    out.grid.cols = cfg.gridIntv * (padCols - 1) + 1;
    out.grid.rows = cfg.gridIntv * (padRows - 1) + 1;
    if (out.grid.cols < 2 || out.grid.rows < 2) {
        throw std::runtime_error("Virtual grid too small; increase chip size or decrease pad pitch.");
    }

    out.grid.dx_m = chip.width_m / static_cast<double>(out.grid.cols - 1);
    out.grid.dy_m = chip.height_m / static_cast<double>(out.grid.rows - 1);

    if (!cfg.mlcfPath.empty()) {
        out.metal = parseMlcf(cfg.mlcfPath);
    } else {
        out.metal = defaultMetalStack();
    }

    // VoltSpot per-layer edge resistance model (PDN_sim.c populate_R_model_PDN):
    // For each physical metal layer:
    // metal_grid_col = floor(cw/(2*pitch)), metal_grid_row = floor(ch/(2*pitch))
    // if X: r = (1/(nc-1))*(nr/metal_grid_row) * (rho*cw/(w*t))
    // if Y: r = (1/(nr-1))*(nc/metal_grid_col) * (rho*ch/(w*t))
    const double cw = chip.width_m;
    const double ch = chip.height_m;
    const double nr = static_cast<double>(out.grid.rows);
    const double nc = static_cast<double>(out.grid.cols);

    out.r.perMetal_edgeR_ohm.resize(out.metal.size(), 0.0);
    double gx = 0.0;
    double gy = 0.0;

    for (const auto& layer : out.metal) {
        if (layer.pitch_m <= 0 || layer.width_m <= 0 || layer.thickness_m <= 0 || layer.rho_ohm_m <= 0) {
            throw std::runtime_error("Invalid metal parameters in mlcf (must be positive).");
        }

        const int metal_grid_col = static_cast<int>(std::floor(cw / (2.0 * layer.pitch_m)));
        const int metal_grid_row = static_cast<int>(std::floor(ch / (2.0 * layer.pitch_m)));
        if (metal_grid_col <= 0 || metal_grid_row <= 0) {
            throw std::runtime_error("Metal pitch too large for chip dimensions; metal_grid_col/row becomes 0.");
        }

        const double baseR_x = layer.rho_ohm_m * cw / (layer.width_m * layer.thickness_m);
        const double baseR_y = layer.rho_ohm_m * ch / (layer.width_m * layer.thickness_m);

        double r_edge = 0.0;
        if (layer.direction == Direction::X) {
            r_edge = (1.0 / (nc - 1.0)) * (nr / static_cast<double>(metal_grid_row)) * baseR_x;
            gx += 1.0 / r_edge;
        } else {
            r_edge = (1.0 / (nr - 1.0)) * (nc / static_cast<double>(metal_grid_col)) * baseR_y;
            gy += 1.0 / r_edge;
        }
        out.r.perMetal_edgeR_ohm[static_cast<size_t>(layer.layerId)] = r_edge;
    }

    if (gx <= 0.0 || gy <= 0.0) {
        throw std::runtime_error("Failed to compute merged rx/ry (no layers in X or Y?)");
    }
    out.r.rx_edge_ohm = 1.0 / gx;
    out.r.ry_edge_ohm = 1.0 / gy;

    return out;
}

