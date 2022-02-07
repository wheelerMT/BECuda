//
// Created by mattw on 07/02/2022.
//

#ifndef BECUDA_WAVEFUNCTION_H
#define BECUDA_WAVEFUNCTION_H

#include <cufft.h>
#include <string>
#include <iostream>
#include "grid.cuh"

// Abstract base class
class Wavefunction
{

};

class Wavefunction2D : Wavefunction
{
private:
    // FFT plans for each wavefunction component
    cufftHandle m_planPlus{}, m_planZero{}, m_planMinus{};

public:
    // Constructor
    explicit Wavefunction2D(Grid2D &grid);

    // Destructor
    ~Wavefunction2D();

    // FFT-related functions
    void generateFFTPlans();
    void executeFFT();
    void executeIFFT();

    // Initial state functions
    void setInitialState(const std::string &gsPhase);

    // Reference to grid object
    Grid2D &grid;

    // Wavefunction components
    cufftComplex *plus{}, *zero{}, *minus{};
    cufftComplex *plus_k{}, *zero_k{}, *minus_k{};  // K-space versions

};

#endif //BECUDA_WAVEFUNCTION_H
