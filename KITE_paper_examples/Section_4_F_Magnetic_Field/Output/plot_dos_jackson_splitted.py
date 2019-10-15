import h5py
import numpy as np
import matplotlib.pyplot as mpl

from matplotlib import ticker
from pybinding.constants import e, hbar


def get_moments_and_scales(file_name):
    file_input = h5py.File(file_name, 'r+')

    # get the DOS moments data from the hdf5
    moments = np.array(file_input['Calculation']['dos']['MU']).flatten()

    # get the energy scale, b is 0 at the moment
    a, b = file_input['EnergyScale'].value, file_input['EnergyShift'].value

    return moments, a, b


def get_size(file_name):
    file_input = h5py.File(file_name, 'r+')

    l_size = np.array(file_input['L']).flatten()
    l1, l2 = l_size[0], l_size[1]
    num_orbitals = file_input['NOrbitals'].value

    return num_orbitals, l1, l2


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
sns.set_context("talk", font_scale=1.3)

# ************************************
# read h5 just to know the number of moments, enters the name of the DOS file
file_name = 'phmag.h5'
moments_KITE, a_scale, b_scale = get_moments_and_scales(file_name)
num_orbitals, _, _ = get_size(file_name)

# this was the grid where the DOS was evaluated
num_points = 5000
energy1 = np.linspace(0.335, 0.37, num_points)
energy2 = np.linspace(-1.17, -1.2, num_points)

# ************************************
# define figure
fig = mpl.figure(figsize=(16, 8))
(ax1, ax2) = fig.subplots(1, 2, sharey='row')

# Plot the DOS function
select_part_mom = int(moments_KITE.shape[0])
for idx_e, axis_s in enumerate([ax2, ax1]):

    dos_jackson = np.loadtxt('dos_Jackson_broad_nM_{}_part_{}.dat'.format(select_part_mom, idx_e))
    axis_s.plot(dos_jackson[:, 0], dos_jackson[:, 1], c='C0', linewidth=2)

# zoom-in / limit the view to different portions of the data
ax1.set_xlim(energy2.min(), energy2.max())  # outliers only
ax1.set_ylim(0, 1.)  # outliers only
ax2.set_xlim(energy1.min(), energy1.max())  # most of the data
ax2.set_ylim(0, 1.)  # outliers only

# hide the spines between ax1 and ax2
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax1.yaxis.tick_left()
ax1.tick_params(labelright=False)
ax2.yaxis.tick_right()

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax2.plot((-d, +d), (-d, +d), **kwargs)

# What's cool about this is that now if we vary the distance between
# ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
# the diagonal lines will move accordingly, and stay right at the tips
# of the spines they are 'breaking'

me = 9.10938356 * 10 ** (-31)  # [kg]
B_field = 7.93532
idx_n = np.arange(100)
Ev = -1.18
Ec = Ev + 1.52

omega_e = e * B_field / me
omega_c = 2.657 * omega_e
omega_v = 2.182 * omega_e

energy_ll_c = Ec + (idx_n + 0.5) * hbar * omega_c
energy_ll_v = Ev - (idx_n + 0.5) * hbar * omega_v

for ellc, ellv in zip(energy_ll_c, energy_ll_v):
    ax1.axvline(x=ellv, linestyle='--', c='k', alpha=0.5)
    ax2.axvline(x=ellc, linestyle='--', c='k', alpha=0.5)

ax1.set_xlabel('E (eV)')
ax1.xaxis.set_label_coords(1.08, -0.08)
ax1.get_xaxis().set_major_formatter(ticker.StrMethodFormatter('{x:,.2f}'))
ax2.get_xaxis().set_major_formatter(ticker.StrMethodFormatter('{x:,.2f}'))
ax1.set_ylabel(r'${\rm states / (atom \cdot eV)}$')

mpl.savefig('phosphorene_DOS_magnetic.png', dpi=640)
mpl.show()
