//
// Created by mattw on 07/02/2022.
//

#include "wavefunction.cuh"


Wavefunction2D::Wavefunction2D(Grid2D &grid) : grid{grid}
{
    // Allocate arrays
    plus = new cufftComplex[grid.nx * grid.ny]{};
    zero = new cufftComplex[grid.nx * grid.ny]{};
    minus = new cufftComplex[grid.nx * grid.ny]{};
    plus_k = new cufftComplex[grid.nx * grid.ny]{};
    zero_k = new cufftComplex[grid.nx * grid.ny]{};
    minus_k = new cufftComplex[grid.nx * grid.ny]{};

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
            plus[j + i * grid.nx] = {0., 0.};
            zero[j + i * grid.nx] = {1., 0.};
            minus[j + i * grid.nx] = {0., 0.};
        }
    }
}

void Wavefunction2D::add_noise(const std::string &components, float mean, float stddev) const
{
    // Construct random generator
    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator{seed1};
    std::normal_distribution<float> norm_dist{mean, stddev};

    if (components == "outer")
    {
        std::cout << "Adding noise...\n";

        // Add noise to outer components
        for (int i = 0; i < grid.nx; i++)
        {
            for (int j = 0; j < grid.ny; j++)
            {
                plus[j + i * grid.nx].x += norm_dist(generator);
                plus[j + i * grid.nx].y += norm_dist(generator);
                minus[j + i * grid.nx].x += norm_dist(generator);
                minus[j + i * grid.nx].y += norm_dist(generator);
            }
        }
    }

    // Add other component combinations as needed
}
