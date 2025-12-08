#ifndef PIXEL_MODEL_H
#define PIXEL_MODEL_H

#include <string>
#include <tuple>
#include <stdexcept>

/**
 * PixelModel represents the resistance model for a pixel in the PDN grid
 * It contains three resistance values: Rx, Ry, Rz (horizontal, vertical, vertical between layers)
 */
class PixelModel {
public:
    PixelModel(const std::string& name = "", float rx = 0.0f, float ry = 0.0f, float rz = 0.0f, float voltage = 0.0f)
        : _name(name), _rx(rx), _ry(ry), _rz(rz), _voltage(voltage) {}

    // Getters
    std::string getName() const { return _name; }
    float getRx() const { return _rx; }
    float getRy() const { return _ry; }
    float getRz() const { return _rz; }
    float getVoltage() const { return _voltage; }
    std::tuple<float, float, float> getResistances() const { return std::make_tuple(_rx, _ry, _rz); }

    // Setters
    void setName(const std::string& name) { _name = name; }
    void setRx(float rx) { 
        if (rx < 0) throw std::invalid_argument("Rx must be positive");
        _rx = rx; 
    }
    void setRy(float ry) { 
        if (ry < 0) throw std::invalid_argument("Ry must be positive");
        _ry = ry; 
    }
    void setRz(float rz) { 
        if (rz < 0) throw std::invalid_argument("Rz must be positive");
        _rz = rz; 
    }
    void setVoltage(float voltage) { 
        if (voltage < 0) throw std::invalid_argument("Voltage must be positive");
        _voltage = voltage; 
    }
    void setResistances(float rx, float ry, float rz) {
        setRx(rx);
        setRy(ry);
        setRz(rz);
    }

    // Scale resistances by a ratio (for sensitivity analysis)
    void scaleResistances(float ratio) {
        _rx *= ratio;
        _ry *= ratio;
        _rz *= ratio;
    }

    void print() const;

private:
    std::string _name;
    float _rx;  // Horizontal resistance
    float _ry;  // Vertical resistance
    float _rz;  // Vertical resistance between layers
    float _voltage;
};

#endif // PIXEL_MODEL_H

