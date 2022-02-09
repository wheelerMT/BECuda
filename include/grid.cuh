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

    // Grid points
    const int xNumGridPts{};
    const int yNumGridPts{};

    // Grid spacing
    double xGridSpacing{};
    double yGridSpacing{};
    double xFourierGridSpacing{};
    double yFourierGridSpacing{};

    // Grid lengths
    double xLengthOfBox{};
    double yLengthOfBox{};

    // Grids
    double *xMesh{};
    double *yMesh{};
    double *xFourierMesh{};
    double *yFourierMesh{};
    double *wavenumberMesh{};  // Square of wave number, |k|^2 = kx^2 + ky^2

    // --------------------
    // Methods
    // --------------------

    // Constructors
    Grid2D(int xNumGridPts, int yNumGridPts, double xGridSpacing, double yGridSpacing);

    // Destructor
    ~Grid2D();

    // FFT functions
    void fftshift() const;

};

#endif //BECUDA_GRID_CUH
