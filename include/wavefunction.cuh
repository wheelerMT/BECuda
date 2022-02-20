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

    void generateFFTPlans();

    void setPolarInitialState() const;

public:
    Grid2D grid;
    double *trappingPotential{};

    // Device component variables
    cufftComplex *plusComponent{}, *zeroComponent{}, *minusComponent{};
    cufftComplex *plusFourierComponent{}, *zeroFourierComponent{}, *minusFourierComponent{};

    explicit Wavefunction2D(Grid2D &grid);

    ~Wavefunction2D();

    void setTrappingPotential(const double *trappingPotential) const;

    void setInitialState(const std::string &groundState) const;

    void addNoiseToComponents(std::string const &components, float mean, float stddev) const;

    void fft() const;

    void ifft() const;

};

#endif //BECUDA_WAVEFUNCTION_H
