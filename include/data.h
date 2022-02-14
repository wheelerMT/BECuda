//
// Created by mattw on 10/02/2022.
//

#ifndef BECUDA_DATA_H
#define BECUDA_DATA_H

#include "highfive/H5File.hpp"
#include "grid.cuh"
#include "wavefunction.cuh"

struct Spin1Parameters
{
    Spin1Parameters() = default;
    Spin1Parameters(double c0, double c2, double p, double q, int nt, double dt);

    double spinIndependentInteraction;
    double spinDependentInteraction;
    double linearZeeman;
    double quadraticZeeman;

    int numOfTimeSteps;
    double timeStep;
};


class Spin1DataManager2D
{
private:
    unsigned int saveIndex{0};

    void saveParametersData(const Spin1Parameters &params, const Grid2D &grid);

    void generateWavefunctionDatasets(const Grid2D &grid);

public:
    // Constructor
    Spin1DataManager2D(const std::string &filename, const Spin1Parameters &params, const Grid2D &grid);

    void saveWavefunctionData(Wavefunction2D &wavefunction);

    std::string filename;
    HighFive::File file;
};

#endif //BECUDA_DATA_H
