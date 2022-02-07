//
// Created by mattw on 07/02/2022.
//

#include "wavefunction.cuh"


Wavefunction2D::Wavefunction2D(Grid2D &grid) : grid{grid}
{
    // Allocate arrays
    plus = new cufftComplex[grid.nx * grid.ny];
    zero = new cufftComplex[grid.nx * grid.ny];
    minus = new cufftComplex[grid.nx * grid.ny];
    plus_k = new cufftComplex[grid.nx * grid.ny];
    zero_k = new cufftComplex[grid.nx * grid.ny];
    minus_k = new cufftComplex[grid.nx * grid.ny];

    // Initialise FFT plans
    generateFFTPlans();
}

Wavefunction2D::~Wavefunction2D()
{
    // This needs to call appropriate functions
    // to de-allocate arrays on device memory
}

void Wavefunction2D::generateFFTPlans()
{
    // Generate CUDA FFT plans for each component check for errors
    if (cufftPlan2d(&m_planPlus, grid.nx, grid.ny, CUFFT_C2C) != CUFFT_SUCCESS)
    {
        std::cerr << "CUFFT error: Plan creation failed...\n";
        return;
    };
    if (cufftPlan2d(&m_planZero, grid.nx, grid.ny, CUFFT_C2C) != CUFFT_SUCCESS)
    {
        std::cerr << "CUFFT error: Plan creation failed...\n";
        return;
    };
    if (cufftPlan2d(&m_planMinus, grid.nx, grid.ny, CUFFT_C2C) != CUFFT_SUCCESS)
    {
        std::cerr << "CUFFT error: Plan creation failed...\n";
        return;
    };
}

void Wavefunction2D::setInitialState(const std::string &gsPhase)
{
    if (gsPhase == "polar")
    {
        for (int i = 0; i < grid.nx; ++i)
        {
            for (int j = 0; j < grid.ny; ++j)
            {
                plus[j + i * grid.nx] = {0., 0.};
                zero[j + i * grid.nx] = {1., 0.};
                minus[j + i * grid.nx] = {0., 0.};
            }
        }
    }
}
