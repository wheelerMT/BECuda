#include <vector>
#include "constants.h"
#include "grid.cuh"

void Grid2D::constructGridParameters()
{
    // K-space grid spacing
    xFourierGridSpacing = PI / (xNumGridPts / 2. * xGridSpacing);
    yFourierGridSpacing = PI / (yNumGridPts / 2. * yGridSpacing);

    // Set length of sides of box
    xLengthOfBox = xNumGridPts * xGridSpacing;
    yLengthOfBox = yNumGridPts * yGridSpacing;
}

void Grid2D::constructGrids()
{
    // Allocate meshgrid arrays
    xMesh = new double[xNumGridPts * yNumGridPts];
    yMesh = new double[xNumGridPts * yNumGridPts];
    xFourierMesh = new double[xNumGridPts * yNumGridPts];
    yFourierMesh = new double[xNumGridPts * yNumGridPts];
    wavenumberMesh = new double[xNumGridPts * yNumGridPts];

    // Construct grids
    for (int i = 0; i < xNumGridPts; ++i)
    {
        for (int j = 0; j < yNumGridPts; ++j)
        {
            xMesh[j + i * yNumGridPts] = (j - xNumGridPts / 2.) * xGridSpacing;
            xFourierMesh[j + i * yNumGridPts] = (j - yNumGridPts / 2.) * xFourierGridSpacing;
            yMesh[j + i * yNumGridPts] = (j - xNumGridPts / 2.) * yGridSpacing;
            yFourierMesh[j + i * yNumGridPts] = (j - yNumGridPts / 2.) * yFourierGridSpacing;
            wavenumberMesh[j + i * yNumGridPts] = std::pow(xFourierMesh[j + i * yNumGridPts], 2)
                                                  + std::pow(yFourierMesh[j + i * yNumGridPts], 2);
        }
    }
}

void Grid2D::fftshift() const
{
    /*
    * Shifts the zero-frequency component to the center
    * of the spectrum.
    */

    std::vector<double> xFourierMeshCopy(xNumGridPts * yNumGridPts);
    std::vector<double> yFourierMeshCopy(xNumGridPts * yNumGridPts);

    for (int i = 0; i < xNumGridPts; ++i)
    {
        for (int j = 0; j < yNumGridPts; ++j)
        {
            xFourierMeshCopy[j + i * xNumGridPts] = xFourierMesh[j + i * xNumGridPts];
            yFourierMeshCopy[j + i * xNumGridPts] = yFourierMesh[j + i * xNumGridPts];
        }
    }

    // Reverse each row
    for (int i = 0; i < xNumGridPts; ++i)
    {
        for (int j = 0; j < xNumGridPts; ++j)
        {
            if (j < xNumGridPts / 2)
            {
                xFourierMesh[j + i * yNumGridPts] = xFourierMeshCopy[xNumGridPts / 2 + j + i * yNumGridPts];
                yFourierMesh[j + i * yNumGridPts] = yFourierMeshCopy[xNumGridPts / 2 + j + i * yNumGridPts];
            } else if (j >= xNumGridPts / 2)
            {
                xFourierMesh[j + i * yNumGridPts] = xFourierMeshCopy[j - xNumGridPts / 2 + i * yNumGridPts];
                yFourierMesh[j + i * yNumGridPts] = yFourierMeshCopy[j - xNumGridPts / 2 + i * yNumGridPts];
            }

        }
    }

    for (int i = 0; i < xNumGridPts; ++i)
    {
        for (int j = 0; j < yNumGridPts; ++j)
        {
            xFourierMeshCopy[j + i * xNumGridPts] = xFourierMesh[j + i * xNumGridPts];
            yFourierMeshCopy[j + i * xNumGridPts] = yFourierMesh[j + i * xNumGridPts];
        }
    }

    // Reverse each column
    for (int i = 0; i < xNumGridPts; ++i)
    {
        for (int j = 0; j < xNumGridPts; ++j)
        {
            if (j < xNumGridPts / 2)
            {
                xFourierMesh[i + j * xNumGridPts] = xFourierMeshCopy[(xNumGridPts / 2 + j) * xNumGridPts + i];
                yFourierMesh[i + j * xNumGridPts] = yFourierMeshCopy[(xNumGridPts / 2 + j) * xNumGridPts + i];
            } else if (j >= xNumGridPts / 2)
            {
                xFourierMesh[i + j * xNumGridPts] = xFourierMeshCopy[(j - xNumGridPts / 2) * xNumGridPts + i];
                yFourierMesh[i + j * xNumGridPts] = yFourierMeshCopy[(j - xNumGridPts / 2) * xNumGridPts + i];
            }
        }
    }

    // Re-update wavenumber, k
    for (int i = 0; i < xNumGridPts; ++i)
    {
        for (int j = 0; j < xNumGridPts; ++j)
        {
            wavenumberMesh[j + i * xNumGridPts] =
                    std::pow(xFourierMesh[j + i * xNumGridPts], 2) + std::pow(yFourierMesh[j + i * xNumGridPts], 2);
        }
    }
}

Grid2D::Grid2D(int xNumGridPts, int yNumGridPts, double xGridSpacing, double yGridSpacing)
        : xNumGridPts{xNumGridPts}, yNumGridPts{yNumGridPts},
          xGridSpacing{xGridSpacing}, yGridSpacing{yGridSpacing}
{
    constructGridParameters();
    constructGrids();
}

Grid2D::Grid2D(const Grid2D &grid)
        : xNumGridPts{grid.xNumGridPts}, yNumGridPts{grid.yNumGridPts},
          xGridSpacing{grid.xGridSpacing}, yGridSpacing{grid.yGridSpacing}
{
    constructGridParameters();
    constructGrids();
}

Grid2D::~Grid2D()
{
    delete[] xMesh;
    delete[] yMesh;
    delete[] xFourierMesh;
    delete[] yFourierMesh;
    delete[] wavenumberMesh;
}
