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

public:
    // --------------------
    // Variables
    // --------------------

    // Grid points
    const int nx{};
    const int ny{};

    // Grid spacing
    double dx{};
    double dy{};
    double dkx{};
    double dky{};

    // Grid lengths
    double lenX{};
    double lenY{};

    // Grids
    double *X{};
    double *Y{};
    double *Kx{};
    double *Ky{};
    double *K{};  // Square of wave number, |k|^2 = kx^2 + ky^2

    // --------------------
    // Methods
    // --------------------

    // Constructors
    Grid2D(int nx, int ny, double dx, double dy);

    // Destructor
    ~Grid2D();

    // FFT functions
    void fftshift() const;

};

#endif //BECUDA_GRID_CUH
