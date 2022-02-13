//
// Created by mattw on 10/02/2022.
//

#ifndef BECUDA_DATA_H
#define BECUDA_DATA_H

#include "highfive/H5File.hpp"

struct Spin1Parameters
{
    Spin1Parameters(double c0, double c2, double p, double q, int nt, double dt);

    double spinIndependentInteraction;
    double spinDependentInteraction;
    double linearZeeman;
    double quadraticZeeman;

    int numOfTimeSteps;
    double timeStep;
};


class DataManager
{

};

#endif //BECUDA_DATA_H
