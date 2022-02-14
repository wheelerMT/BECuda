import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np

"""Visualises the component densities for each frame of data."""


def convert_to_complex(wavefunction_dataset):
    wavefunction = np.empty((nx, ny, wavefunction_dataset.shape[-1]), dtype='complex64')
    for frame in range(wavefunction_dataset.shape[-1]):
        for i in range(nx):
            for j in range(ny):
                wavefunction[i, j, frame] = \
                    complex(wavefunction_dataset[i, j, frame][0], wavefunction_dataset[i, j, frame][1])

    return wavefunction


filename = 'test'
data_file = h5py.File(f'../data/{filename}.hdf5', 'r')

nx, ny = data_file['grid/nx'][...], data_file['grid/ny'][...]
dx, dy = data_file['grid/dx'][...], data_file['grid/dy'][...]
X, Y = np.meshgrid(np.arange(-nx // 2, nx // 2) * dx, np.arange(-ny // 2, ny // 2) * dy)

psi_plus_data = data_file['wavefunction/psi_plus']
psi_zero_data = data_file['wavefunction/psi_zero']
psi_minus_data = data_file['wavefunction/psi_minus']

psi_plus_data = np.reshape(psi_plus_data, (nx, ny, psi_plus_data.shape[-1]))
psi_zero_data = np.reshape(psi_zero_data, (nx, ny, psi_zero_data.shape[-1]))
psi_minus_data = np.reshape(psi_minus_data, (nx, ny, psi_minus_data.shape[-1]))

psi_plus = convert_to_complex(psi_plus_data)
psi_zero = convert_to_complex(psi_zero_data)
psi_minus = convert_to_complex(psi_minus_data)

fig = plt.figure(figsize=(12, 4))
plotting_grid = AxesGrid(fig, 111, nrows_ncols=(1, 3), cbar_mode='single', cbar_location='right',
                         axes_pad=0.1, cbar_pad=0.1)

for axis in plotting_grid:
    if axis == plotting_grid[0]:
        axis.set_ylabel(r'$\tilde{y}$')
        axis.set_title(r'$|\psi_+|^2$')
    if axis == plotting_grid[1]:
        axis.set_title(r'$|\psi_0|^2$')
    if axis == plotting_grid[2]:
        axis.set_title(r'$|\psi_-|^2$')
    axis.set_xlabel(r'$\tilde{x}$')
    axis.set_aspect('equal')

# Do initial plot
psi_plus_plot = plotting_grid[0].pcolormesh(X, Y, abs(psi_plus[:, :, 0]) ** 2, vmin=0, vmax=1)
psi_zero_plot = plotting_grid[1].pcolormesh(X, Y, abs(psi_zero[:, :, 0]) ** 2, vmin=0, vmax=1)
psi_minus_plot = plotting_grid[2].pcolormesh(X, Y, abs(psi_minus[:, :, 0]) ** 2, vmin=0, vmax=1)

# Add colour bar to last plotting axis
cbar = plotting_grid[-1].cax.colorbar(psi_minus_plot)

# Update plots
for k in range(psi_plus.shape[-1]):
    psi_plus_plot = plotting_grid[0].pcolormesh(X, Y, abs(psi_plus[:, :, k]) ** 2)
    psi_zero_plot = plotting_grid[1].pcolormesh(X, Y, abs(psi_zero[:, :, k]) ** 2)
    psi_minus_plot = plotting_grid[2].pcolormesh(X, Y, abs(psi_minus[:, :, k]) ** 2)
    plt.pause(0.05)
plt.show()
