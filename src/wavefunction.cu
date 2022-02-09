#include "wavefunction.cuh"


Wavefunction2D::Wavefunction2D(Grid2D &grid) : grid{grid}
{
    plusComponent = new cufftComplex[grid.xNumGridPts * grid.yNumGridPts]{};
    zeroComponent = new cufftComplex[grid.xNumGridPts * grid.yNumGridPts]{};
    minusComponent = new cufftComplex[grid.xNumGridPts * grid.yNumGridPts]{};
    plusFourierComponent = new cufftComplex[grid.xNumGridPts * grid.yNumGridPts]{};
    zeroFourierComponent = new cufftComplex[grid.xNumGridPts * grid.yNumGridPts]{};
    minusFourierComponent = new cufftComplex[grid.xNumGridPts * grid.yNumGridPts]{};
}

Wavefunction2D::~Wavefunction2D()
{
    cudaFree(plusComponent);
    cudaFree(zeroComponent);
    cudaFree(minusComponent);
    cudaFree(plusFourierComponent);
    cudaFree(zeroFourierComponent);
    cudaFree(minusFourierComponent);
}

void Wavefunction2D::generateFFTPlans()
{
    cufftPlan2d(&fftPlan, grid.xNumGridPts, grid.yNumGridPts, CUFFT_C2C);
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

void Wavefunction2D::addNoise(const std::string &components, float mean, float stddev) const
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
