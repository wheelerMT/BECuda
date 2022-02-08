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
    cufftHandle m_FFTPlan{};

public:
    // --------------------
    // Variables
    // --------------------

    // Reference to grid object
    Grid2D &grid;

    // Wavefunction components
    cufftComplex *plus{}, *zero{}, *minus{};
    cufftComplex *plus_k{}, *zero_k{}, *minus_k{};  // K-space versions

    // --------------------
    // Methods
    // --------------------

    // Constructor
    explicit Wavefunction2D(Grid2D &grid);

    // Destructor
    ~Wavefunction2D();

    // FFT-related functions
    void generateFFTPlans();

    void fft();

    void ifft();

    // Initial state functions
    void setInitialState(const std::string &gsPhase) const;
    void add_noise(std::string const &components, float mean, float stddev) const;
};

#endif //BECUDA_WAVEFUNCTION_H
