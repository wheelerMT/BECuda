//
// Created by mattw on 27/01/2022.
//

#ifndef BECUDA_GRID_CUH
#define BECUDA_GRID_CUH

// Abstract base class
class Grid
{

};

class Grid2D : Grid
{
private:
    void constructGridParameters();

    void constructGrids();

    void fftshift();

public:
    // Grid points
    const unsigned int nx{};
    const unsigned int ny{};

    // Grid spacing
    double dx{};
    double dy{};
    double dkx{};
    double dky{};

    // Grid lengths
    double len_x{};
    double len_y{};

    // Grids
    doubleArray_t X{};
    doubleArray_t Y{};
    doubleArray_t Kx{};
    doubleArray_t Ky{};
    doubleArray_t K{};

    // Constructors
    Grid2D(unsigned int nx, unsigned int ny, double dx, double dy);

};

#endif //BECUDA_GRID_CUH
