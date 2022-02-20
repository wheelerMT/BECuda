#include "wavefunction.cuh"


Wavefunction2D::Wavefunction2D(Grid2D &grid) : grid{grid}
{
    cudaMallocManaged(&plusComponent, grid.xNumGridPts * grid.yNumGridPts * sizeof(cufftComplex));
    cudaMallocManaged(&zeroComponent, grid.xNumGridPts * grid.yNumGridPts * sizeof(cufftComplex));
    cudaMallocManaged(&minusComponent, grid.xNumGridPts * grid.yNumGridPts * sizeof(cufftComplex));
    cudaMallocManaged(&plusFourierComponent, grid.xNumGridPts * grid.yNumGridPts * sizeof(cufftComplex));
    cudaMallocManaged(&zeroFourierComponent, grid.xNumGridPts * grid.yNumGridPts * sizeof(cufftComplex));
    cudaMallocManaged(&minusFourierComponent, grid.xNumGridPts * grid.yNumGridPts * sizeof(cufftComplex));

    trappingPotential = new double[grid.xNumGridPts * grid.yNumGridPts]{};

    generateFFTPlans();
}

void Wavefunction2D::generateFFTPlans()
{
    cufftPlan2d(&fftPlan, grid.xNumGridPts, grid.yNumGridPts, CUFFT_C2C);
}

Wavefunction2D::~Wavefunction2D()
{
    cufftDestroy(fftPlan);

    cudaFree(plusComponent);
    cudaFree(zeroComponent);
    cudaFree(minusComponent);
    cudaFree(plusFourierComponent);
    cudaFree(zeroFourierComponent);
    cudaFree(minusFourierComponent);

    delete[] trappingPotential;

}

void Wavefunction2D::setTrappingPotential(const double *newTrappingPotential) const
{
    for (int i = 0; i < grid.xNumGridPts; ++i)
    {
        for (int j = 0; j < grid.xNumGridPts; ++j)
        {
            trappingPotential[j + i * grid.yNumGridPts] = newTrappingPotential[j + i * grid.yNumGridPts];
        }
    }
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
    for (int i = 0; i < grid.xNumGridPts; ++i)
    {
        for (int j = 0; j < grid.yNumGridPts; ++j)
        {
            plusComponent[j + i * grid.xNumGridPts] = {0., 0.};
            zeroComponent[j + i * grid.xNumGridPts] = {1., 0.};
            minusComponent[j + i * grid.xNumGridPts] = {0., 0.};
        }
    }
}

void Wavefunction2D::addNoiseToComponents(const std::string &components, float mean, float stddev) const
{
    // Construct random generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator{seed};
    std::normal_distribution<float> norm_dist{mean, stddev};

    if (components == "outer")
    {
        std::cout << "Adding noise...\n";
        for (int i = 0; i < grid.xNumGridPts; i++)
        {
            for (int j = 0; j < grid.yNumGridPts; j++)
            {
                plusComponent[j + i * grid.xNumGridPts].x += norm_dist(generator);
                plusComponent[j + i * grid.xNumGridPts].y += norm_dist(generator);
                minusComponent[j + i * grid.xNumGridPts].x += norm_dist(generator);
                minusComponent[j + i * grid.xNumGridPts].y += norm_dist(generator);
            }
        }
    }

    // Add other component combinations as needed
}

void Wavefunction2D::fft() const
{
    cufftExecC2C(fftPlan, plusComponent, plusFourierComponent, CUFFT_FORWARD);
    cufftExecC2C(fftPlan, zeroComponent, zeroFourierComponent, CUFFT_FORWARD);
    cufftExecC2C(fftPlan, minusComponent, minusFourierComponent, CUFFT_FORWARD);
    cudaDeviceSynchronize();
}

void Wavefunction2D::ifft() const
{
    cufftExecC2C(fftPlan, plusFourierComponent, plusComponent, CUFFT_INVERSE);
    cufftExecC2C(fftPlan, zeroFourierComponent, zeroComponent, CUFFT_INVERSE);
    cufftExecC2C(fftPlan, minusFourierComponent, minusComponent, CUFFT_INVERSE);
    cudaDeviceSynchronize();
}