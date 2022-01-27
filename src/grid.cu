//
// Created by mattw on 27/01/2022.
//

#include "grid.cuh"

void Grid2D::constructGridParameters()
{

}

void Grid2D::constructGrids()
{

}

void Grid2D::fftshift()
{

}

Grid2D::Grid2D(unsigned int nx, unsigned int ny, double dx, double dy)
        : nx{nx}, ny{ny}, dx{dx}, dy{dy}
{
    constructGridParameters();
    constructGrids();
}

