//
// Created by mattw on 27/01/2022.
//

#include <vector>
#include "constants.h"
#include "grid.cuh"

void Grid2D::constructGridParameters()
{
    // K-space grid spacing
    dkx = PI / (nx / 2. * dx);
    dky = PI / (ny / 2. * dy);

    // Set length of sides of box
    lenX = nx * dx;
    lenY = ny * dy;
}

void Grid2D::constructGrids()
{
    // Allocate meshgrid arrays
    X = new double[nx * ny];
    Y = new double[nx * ny];
    Kx = new double[nx * ny];
    Ky = new double[nx * ny];
    K = new double[nx * ny];

    // Construct grids
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            X[j + i * ny] = (j - nx / 2.) * dx;
            Kx[j + i * ny] = (j - ny / 2.) * dkx;
            Y[j + i * ny] = (j - nx / 2.) * dy;
            Ky[j + i * ny] = (j - ny / 2.) * dky;
            K[j + i * ny] = std::pow(Kx[j + i * ny], 2) + std::pow(Ky[j + i * ny], 2);
        }
    }
}

void Grid2D::fftshift() const
{
    /*
        Shifts the zero-frequency component to the center
        of the spectrum.
    */

    // Make a copy of K-space arrays
    std::vector<double> kxCopy(nx * ny);
    std::vector<double> kyCopy(nx * ny);

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            kxCopy[j + i * nx] = Kx[j + i * nx];
            kyCopy[j + i * nx] = Ky[j + i * nx];
        }
    }

    // Reverse each row
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            if (j < nx / 2)
            {
                Kx[j + i * ny] = kxCopy[nx / 2 + j + i * ny];
                Ky[j + i * ny] = kyCopy[nx / 2 + j + i * ny];
            } else if (j >= nx / 2)
            {
                Kx[j + i * ny] = kxCopy[j - nx / 2 + i * ny];
                Ky[j + i * ny] = kyCopy[j - nx / 2 + i * ny];
            }

        }
    }

    // Update array copies
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            kxCopy[j + i * nx] = Kx[j + i * nx];
            kyCopy[j + i * nx] = Ky[j + i * nx];
        }
    }

    // Reverse each column
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            if (j < nx / 2)
            {
                Kx[i + j * nx] = kxCopy[(nx / 2 + j) * nx + i];
                Ky[i + j * nx] = kyCopy[(nx / 2 + j) * nx + i];
            } else if (j >= nx / 2)
            {
                Kx[i + j * nx] = kxCopy[(j - nx / 2) * nx + i];
                Ky[i + j * nx] = kyCopy[(j - nx / 2) * nx + i];
            }
        }
    }

    // Re-update wavenumber, k
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < nx; ++j)
        {
            K[j + i * nx] = std::pow(Kx[j + i * nx], 2) + std::pow(Ky[j + i * nx], 2);
        }
    }
}

Grid2D::Grid2D(unsigned int nx, unsigned int ny, double dx, double dy)
        : nx{nx}, ny{ny}, dx{dx}, dy{dy}
{
    constructGridParameters();
    constructGrids();
}

Grid2D::~Grid2D()
{
    // Delete dynamically allocated arrays
    delete[] X;
    delete[] Y;
    delete[] Kx;
    delete[] Ky;
    delete[] K;
}
