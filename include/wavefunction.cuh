//
// Created by mattw on 07/02/2022.
//

#ifndef BECUDA_WAVEFUNCTION_H
#define BECUDA_WAVEFUNCTION_H

#include <cufft.h>
#include <string>
#include <iostream>
#include <random>
#include <chrono>
#include "grid.cuh"

// Abstract base class
class Wavefunction
{

};

class Wavefunction2D : Wavefunction
{
private:
    // FFT plan
    cufftHandle fftPlan{};

    // Initial State functions
    void setPolarInitialState() const;

public:
    // --------------------
    // Variables
    // --------------------

    // Reference to grid object
    Grid2D &grid;

    // Wavefunction components
    cufftComplex *plusComponent{}, *zeroComponent{}, *minusComponent{};
    cufftComplex *plusFourierComponent{}, *zeroFourierComponent{}, *minusFourierComponent{};  // K-space versions

    // --------------------
    // Methods
    // --------------------

    // Constructor
    explicit Wavefunction2D(Grid2D &grid);

    // Destructor
    ~Wavefunction2D();

    // FFT-related functions
    void generateFFTPlans();

    // Initial state functions
    void setInitialState(const std::string &groundState) const;
    void addNoise(std::string const &components, float mean, float stddev) const;
};

#endif //BECUDA_WAVEFUNCTION_H
