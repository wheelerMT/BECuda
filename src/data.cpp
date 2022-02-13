//
// Created by mattw on 10/02/2022.
//

#include "data.h"

Spin1Parameters::Spin1Parameters(double c0, double c2, double p, double q, int nt, double dt)
        : spinIndependentInteraction{c0}, spinDependentInteraction{c2}, linearZeeman{p}, quadraticZeeman{q},
        numOfTimeSteps{nt}, timeStep{dt}
{
}
