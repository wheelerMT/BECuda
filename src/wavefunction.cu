//
// Created by mattw on 07/02/2022.
//

#include "wavefunction.cuh"


Wavefunction2D::Wavefunction2D(Grid2D &grid) : grid{grid}
{
    // Allocate arrays
    plusComponent = new cufftComplex[grid.nx * grid.ny]{};
    zeroComponent = new cufftComplex[grid.nx * grid.ny]{};
    minusComponent = new cufftComplex[grid.nx * grid.ny]{};
    plusFourierComponent = new cufftComplex[grid.nx * grid.ny]{};
    zeroFourierComponent = new cufftComplex[grid.nx * grid.ny]{};
    minusFourierComponent = new cufftComplex[grid.nx * grid.ny]{};
}

Wavefunction2D::~Wavefunction2D()
{
    // Free device memory
    cudaFree(plusComponent);
    cudaFree(zeroComponent);
    cudaFree(minusComponent);
    cudaFree(plusFourierComponent);
    cudaFree(zeroFourierComponent);
    cudaFree(minusFourierComponent);
}

void Wavefunction2D::generateFFTPlans()
{
    // Generate CUDA FFT plans for each component
    cufftPlan2d(&fftPlan, grid.nx, grid.ny, CUFFT_C2C);

}

void Wavefunction2D::setInitialState(const std::string &groundState) const
{
    if (groundState == "polar")
    {
        setPolarInitialState();
    }

    // Can add more ground states as needed
}

void Wavefunction2D::setPolarInitialState() const
{
    for (int i = 0; i < grid.nx; ++i)
    {
        for (int j = 0; j < grid.ny; ++j)
        {
            plusComponent[j + i * grid.nx] = {0., 0.};
            zeroComponent[j + i * grid.nx] = {1., 0.};
            minusComponent[j + i * grid.nx] = {0., 0.};
        }
    }
}

void Wavefunction2D::addNoise(const std::string &components, float mean, float stddev) const
{
    // Construct random generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator{seed};
    std::normal_distribution<float> norm_dist{mean, stddev};

    if (components == "outer")
    {
        std::cout << "Adding noise...\n";

        // Add noise to outer components
        for (int i = 0; i < grid.nx; i++)
        {
            for (int j = 0; j < grid.ny; j++)
            {
                plusComponent[j + i * grid.nx].x += norm_dist(generator);
                plusComponent[j + i * grid.nx].y += norm_dist(generator);
                minusComponent[j + i * grid.nx].x += norm_dist(generator);
                minusComponent[j + i * grid.nx].y += norm_dist(generator);
            }
        }
    }

    // Add other component combinations as needed
}
