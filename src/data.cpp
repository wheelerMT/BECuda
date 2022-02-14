//
// Created by mattw on 10/02/2022.
//

#include "data.h"

// Telling HDF5 how to create a cufftComplex type
HighFive::CompoundType createCompoundCufftComplex()
{
    return {{"real", HighFive::AtomicType<float>{}},
            {"imaginary", HighFive::AtomicType<float>{}}};
}

HIGHFIVE_REGISTER_TYPE(cufftComplex, createCompoundCufftComplex)


Spin1Parameters::Spin1Parameters(double c0, double c2, double p, double q, int nt, double dt)
        : spinIndependentInteraction{c0}, spinDependentInteraction{c2}, linearZeeman{p}, quadraticZeeman{q},
          numOfTimeSteps{nt}, timeStep{dt}
{
}

Spin1DataManager2D::Spin1DataManager2D(const std::string &filename, const Spin1Parameters &params, const Grid2D &grid)
        : filename{filename},
          file{filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate}
{
    saveParametersData(params, grid);

    generateWavefunctionDatasets(grid);
}

void Spin1DataManager2D::saveParametersData(const Spin1Parameters &params, const Grid2D &grid)
{
    // Save condensate and time parameters to file
    file.createDataSet("/parameters/c0", params.spinIndependentInteraction);
    file.createDataSet("/parameters/c2", params.spinDependentInteraction);
    file.createDataSet("/time/nt", params.numOfTimeSteps);
    file.createDataSet("/time/dt", params.timeStep);

    // Save grid parameters to file
    file.createDataSet("/grid/nx", grid.xNumGridPts);
    file.createDataSet("/grid/ny", grid.yNumGridPts);
    file.createDataSet("/grid/dx", grid.xGridSpacing);
    file.createDataSet("/grid/dy", grid.yGridSpacing);

}

void Spin1DataManager2D::generateWavefunctionDatasets(const Grid2D &grid)
{
    // Define dataspaces with arbitrary length of last dimension
    HighFive::DataSpace ds_plus = HighFive::DataSpace(
            {static_cast<unsigned long long>(grid.xNumGridPts * grid.yNumGridPts), 1},
            {static_cast<unsigned long long>(grid.xNumGridPts * grid.yNumGridPts), HighFive::DataSpace::UNLIMITED});
    HighFive::DataSpace ds_zero = HighFive::DataSpace(
            {static_cast<unsigned long long>(grid.xNumGridPts * grid.yNumGridPts), 1},
            {static_cast<unsigned long long>(grid.xNumGridPts * grid.yNumGridPts), HighFive::DataSpace::UNLIMITED});
    HighFive::DataSpace ds_minus = HighFive::DataSpace(
            {static_cast<unsigned long long>(grid.xNumGridPts * grid.yNumGridPts), 1},
            {static_cast<unsigned long long>(grid.xNumGridPts * grid.yNumGridPts), HighFive::DataSpace::UNLIMITED});

    // Use chunking
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(
            std::vector<hsize_t>{static_cast<unsigned long long>((grid.xNumGridPts * grid.yNumGridPts) / 4), 1}));

    // Create cufftComplex type for HDF5
    auto t_cufftComplex = createCompoundCufftComplex();
    t_cufftComplex.commit(file, "cufftComplex");

    // Create wavefunction datasets
    file.createDataSet("/wavefunction/psi_plus", ds_plus,
                       t_cufftComplex, props);
    file.createDataSet("/wavefunction/psi_zero", ds_zero,
                       t_cufftComplex, props);
    file.createDataSet("/wavefunction/psi_minus", ds_minus,
                       t_cufftComplex, props);
}

void Spin1DataManager2D::saveWavefunctionData(Wavefunction2D &wavefunction)
{
    // Load in datasets
    HighFive::DataSet dsPlus = file.getDataSet("/wavefunction/psi_plus");
    HighFive::DataSet dsZero = file.getDataSet("/wavefunction/psi_zero");
    HighFive::DataSet dsMinus = file.getDataSet("/wavefunction/psi_minus");

    // Resize datasets
    dsPlus.resize({static_cast<unsigned long long>(wavefunction.grid.xNumGridPts * wavefunction.grid.yNumGridPts),
                   saveIndex + 1});
    dsZero.resize({static_cast<unsigned long long>(wavefunction.grid.xNumGridPts * wavefunction.grid.yNumGridPts),
                   saveIndex + 1});
    dsMinus.resize({static_cast<unsigned long long>(wavefunction.grid.xNumGridPts * wavefunction.grid.yNumGridPts),
                    saveIndex + 1});

    // FFT so we update real-space arrays
//    wavefunction.ifft();

    // Save new wavefunction data
    dsPlus.select({0, saveIndex},
                  {static_cast<unsigned long long>(wavefunction.grid.xNumGridPts * wavefunction.grid.yNumGridPts),
                   1}).write(wavefunction.h_plusComponent);
    dsZero.select({0, saveIndex},
                  {static_cast<unsigned long long>(wavefunction.grid.xNumGridPts * wavefunction.grid.yNumGridPts),
                   1}).write(wavefunction.h_zeroComponent);
    dsMinus.select({0, saveIndex},
                   {static_cast<unsigned long long>(wavefunction.grid.xNumGridPts * wavefunction.grid.yNumGridPts),
                    1}).write(wavefunction.h_minusComponent);

    saveIndex += 1;
}
