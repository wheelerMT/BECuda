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
    cufftHandle fftPlan{};

    void setPolarInitialState() const;

public:
    Grid2D grid;

    cufftComplex *plusComponent{}, *zeroComponent{}, *minusComponent{};
    cufftComplex *plusFourierComponent{}, *zeroFourierComponent{}, *minusFourierComponent{};  // K-space versions

    explicit Wavefunction2D(Grid2D &grid);

    ~Wavefunction2D();

    void generateFFTPlans();

    void setInitialState(const std::string &groundState) const;

    void addNoiseToComponents(std::string const &components, float mean, float stddev) const;
};

#endif //BECUDA_WAVEFUNCTION_H
