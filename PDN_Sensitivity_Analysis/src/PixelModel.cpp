#include "PixelModel.h"
#include <iostream>
#include <iomanip>

void PixelModel::print() const {
    std::cout << "Pixel Model: " << _name << std::endl;
    std::cout << "  Rx: " << std::fixed << std::setprecision(6) << _rx << " Ω" << std::endl;
    std::cout << "  Ry: " << std::fixed << std::setprecision(6) << _ry << " Ω" << std::endl;
    std::cout << "  Rz: " << std::fixed << std::setprecision(6) << _rz << " Ω" << std::endl;
    std::cout << "  Voltage: " << std::fixed << std::setprecision(3) << _voltage << " V" << std::endl;
}

