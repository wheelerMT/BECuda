//
// Created by mattw on 07/02/2022.
//

#include "wavefunction.cuh"


Wavefunction2D::Wavefunction2D(Grid2D &grid) : grid{grid}
{
    // Allocate arrays on device
    cudaMalloc(&plus, grid.nx * grid.ny * sizeof(cufftComplex));
    cudaMalloc(&zero, grid.nx * grid.ny * sizeof(cufftComplex));
    cudaMalloc(&minus, grid.nx * grid.ny * sizeof(cufftComplex));
    cudaMalloc(&plus_k, grid.nx * grid.ny * sizeof(cufftComplex));
    cudaMalloc(&zero_k, grid.nx * grid.ny * sizeof(cufftComplex));
    cudaMalloc(&minus_k, grid.nx * grid.ny * sizeof(cufftComplex));

    // Initialise FFT plans
    generateFFTPlans();
}

Wavefunction2D::~Wavefunction2D()
{
    // Free device memory
    cudaFree(plus);
    cudaFree(zero);
    cudaFree(minus);
    cudaFree(plus_k);
    cudaFree(zero_k);
    cudaFree(minus_k);
}

void Wavefunction2D::generateFFTPlans()
{
    // Generate CUDA FFT plans for each component
    cufftPlan2d(&m_FFTPlan, grid.nx, grid.ny, CUFFT_C2C);

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
