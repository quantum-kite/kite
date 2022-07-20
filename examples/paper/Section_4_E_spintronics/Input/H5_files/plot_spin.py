import numpy as np
import matplotlib.pyplot as mpl
import h5py
import matplotlib.gridspec as gridspec

from pybinding.constants import hbar


def read_component(path, file, component):
    name = path + 'gaussian_wavepacket_' + file + '.h5'
    file_input = h5py.File(name, 'r+')

    scaled_time = float(file_input['Calculation']['gaussian_wave_packet']['timestep'].value)
    scale_a = float(file_input['EnergyScale'].value)
    deltaT = scaled_time * hbar / scale_a

    # get the S data from the hdf5
    spin = np.array(file_input['Calculation']['gaussian_wave_packet'][component]).flatten()

    num_points = spin.shape[0]
    timesteps = np.linspace(0, deltaT * num_points, num_points)

    return spin, timesteps


# ************************************
# keep these definitions for kite website
import seaborn as sns

mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['savefig.dpi'] = 100
sns.set_style("white")
# Kite color scheme
colors = ["dusty purple", "faded green", "windows blue", "amber", "greyish"]
current_palette = sns.xkcd_palette(colors)
sns.set_palette(current_palette)
sns.set_style("ticks")
sns.set_context("talk", font_scale=1)

fig = mpl.figure(figsize=(16, 12))
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.1, right=0.95, bottom=0.15, top=0.9, wspace=0.1, hspace=0.2)
lambda_block = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[:, :], wspace=0.2, hspace=0.05)
# Reconstruct the DOS function
# ax0 = plt.gca()
ax0 = mpl.subplot(lambda_block[0])
ax1 = mpl.subplot(lambda_block[1])

# read Sx_and
sx_and, timesteps_sx_and = read_component('anderson_disorder/', 'Sx', 'Sx')
# read Sz_and
sz_and, timesteps_sz_and = read_component('anderson_disorder/', 'Sz', 'Sz')

# read Sx and + resonant
sx_and_res, timesteps_sx_and_res = read_component('resonant_scatterers/resonant_concentration/', 'Sx', 'Sx')
# read Sz and + resonant
sz_and_res, timesteps_sz_and_res = read_component('resonant_scatterers/resonant_concentration/', 'Sz', 'Sz')

ax0.plot(timesteps_sx_and * 1e12, sx_and.real, c='C0', label='Sx Anderson', linestyle='--')
ax0.plot(timesteps_sx_and_res * 1e12, sx_and_res.real, c='C2', label='Sx Resonant')

np.savetxt('sx_anderson.dat', np.column_stack((timesteps_sx_and * 1e12, sx_and.real)))
np.savetxt('sz_anderson.dat', np.column_stack((timesteps_sz_and * 1e12, sz_and.real)))
np.savetxt('sx_magnetic_anderson.dat', np.column_stack((timesteps_sx_and_res * 1e12, sx_and_res.real)))
np.savetxt('sz_magnetic_anderson.dat', np.column_stack((timesteps_sz_and_res * 1e12, sz_and_res.real)))

ax1.plot(timesteps_sz_and * 1e12, sz_and.real, c='C1', label='Sz Anderson', linestyle='--')
ax1.plot(timesteps_sz_and_res * 1e12, sz_and_res.real, c='C3', label='Sz Resonant')

ax0.text(0.3, 0.95, 'Spin precession of Sx polarization', horizontalalignment='center',
         verticalalignment='center', transform=ax0.transAxes)
ax1.text(0.3, 0.95, 'Spin precession of Sz polarization', horizontalalignment='center',
         verticalalignment='center', transform=ax1.transAxes)

mpl.sca(ax0)
mpl.xlabel('t (ps)')
mpl.legend()

mpl.ylabel(r'Spin precession')
mpl.sca(ax1)
mpl.legend()
mpl.xlabel('t (ps)')

mpl.savefig('spin_precession_comparison.png', dpi=320)
mpl.show()
