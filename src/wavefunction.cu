#include "wavefunction.cuh"


Wavefunction2D::Wavefunction2D(Grid2D &grid) : grid{grid}
{
    h_plusComponent = new cufftComplex[grid.xNumGridPts * grid.yNumGridPts]{};
    h_zeroComponent = new cufftComplex[grid.xNumGridPts * grid.yNumGridPts]{};
    h_minusComponent = new cufftComplex[grid.xNumGridPts * grid.yNumGridPts]{};

    trappingPotential = new double[grid.xNumGridPts * grid.yNumGridPts] {};
}

Wavefunction2D::~Wavefunction2D()
{
    delete[] h_plusComponent;
    delete[] h_zeroComponent;
    delete[] h_minusComponent;

    cudaFree(plusComponent);
    cudaFree(zeroComponent);
    cudaFree(minusComponent);
    cudaFree(plusFourierComponent);
    cudaFree(zeroFourierComponent);
    cudaFree(minusFourierComponent);
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
            h_plusComponent[j + i * grid.xNumGridPts] = {0., 0.};
            h_zeroComponent[j + i * grid.xNumGridPts] = {1., 0.};
            h_minusComponent[j + i * grid.xNumGridPts] = {0., 0.};
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
                h_plusComponent[j + i * grid.xNumGridPts].x += norm_dist(generator);
                h_plusComponent[j + i * grid.xNumGridPts].y += norm_dist(generator);
                h_minusComponent[j + i * grid.xNumGridPts].x += norm_dist(generator);
                h_minusComponent[j + i * grid.xNumGridPts].y += norm_dist(generator);
            }
        }
    }

    // Add other component combinations as needed
}
