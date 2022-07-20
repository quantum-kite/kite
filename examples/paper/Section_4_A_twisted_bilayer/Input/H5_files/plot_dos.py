import h5py
import matplotlib.gridspec as gridspec
import numpy as np
import copy
import matplotlib.pyplot as mpl

from pybinding.repository import graphene


def lorentz_kernel(num_moments, broadening):
    """Define a Lorentz damping kernel

    Parameters
    ----------
    num_moments : int
        Number of expansion polynomials.
    broadening : float
        Defines the resolution.
    """
    lambda_b = num_moments * broadening
    ns = np.arange(num_moments)
    return np.sinh(lambda_b * (1 - ns / num_moments)) / np.sinh(lambda_b)


def jackson_kernel(num_moments):
    """Define a Jackson damping kernel

    Parameters
    ----------
    num_moments : int
        Number of expansion polynomials.
    """

    ns = np.arange(num_moments)
    f = np.pi / (num_moments + 1)
    return ((num_moments + 1 - ns) * np.cos(ns * f) + np.sin(ns * f) / np.tan(f)) / (num_moments + 1)


def reconstruct_dos(moments, energy, scale_factors, **kwargs):
    """Reconstruct Density of States

    Parameters
    ----------
    moments : np.array
        Array of moments
    energy : np.array
        Array of energy points for evaluating the DOS
    scale_factors : tuple or list
        Two scaling factors
    **kwargs :
        Broadening
    """
    moments0 = copy.deepcopy(moments)
    broadening = kwargs.get('broadening', None)

    num_moments = moments0.shape[0]

    a, b = scale_factors

    # apply kernel
    if broadening:
        kernel = lorentz_kernel(num_moments, broadening)
    else:
        kernel = jackson_kernel(num_moments)

    moments0 *= kernel

    scaled_energy = (energy - b) / a
    ns = np.arange(num_moments)
    k = 2 / (a * np.pi * np.sqrt(1 - scaled_energy ** 2))
    chebyshev = np.cos(ns * np.arccos(scaled_energy[:, np.newaxis])) / (1 + np.uint32(ns == 0))

    dos = k * np.sum(moments0.real * chebyshev, axis=1)
    return dos


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


def monolayer_density(energy):
    # define the density of monolayer graphene
    lattice_gr = graphene.monolayer()
    l1gr, l2gr = lattice_gr.vectors
    l1gr, l2gr = l1gr[0:2], l2gr[0:2]

    Ac = np.linalg.norm(np.cross(l1gr, l2gr))
    vf = 3 / (2) * abs(graphene.t) * graphene.a_cc  #: [eV*nm] Fermi velocity

    dos_slg = 1 * Ac / np.pi * np.abs(energy) / vf ** 2

    return dos_slg


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

# ************************************
# read unrelaxed
file_name = 'unrelaxed/dos_tblg_1.050_12000_moments_more_neigh.h5'
moments_unrlx, a_unrlx, b_unrlx = get_moments_and_scales(file_name)
num_orbitals_unrlx, l1_unrlx, l2_unrlx = get_size(file_name)

# read relaxed
file_name = 'relaxed/dos_relaxed_tblg_1.050_12000_moments_more_neigh.h5'
moments_rlx, a_rlx, b_rlx = get_moments_and_scales(file_name)
num_orbitals_rlx, l1_rlx, l2_rlx = get_size(file_name)

# define energy grid
num_points = 3000
energy = np.linspace(-1.02, 1.02, num_points)

# ************************************
# define figure
fig = mpl.figure(figsize=(16, 12))
fig.suptitle(r'$\theta=1.05^\circ$' + '\n' + 'num orbitals {}'.format(num_orbitals_unrlx))
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.1, right=0.95, bottom=0.15, top=0.9, wspace=0.1, hspace=0.2)
lambda_block = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[:, :], wspace=0.2, hspace=0.05)

ax0 = mpl.subplot(lambda_block[0])
ax1 = mpl.subplot(lambda_block[1])

dos_slg = monolayer_density(energy)
ED=0.786
# Reconstruct the DOS function and plot UNRLX
dos = reconstruct_dos(moments_unrlx, energy, (a_unrlx, b_unrlx), broadening=10e-3 / a_unrlx)
ax0.plot(energy-ED, 2 * dos, c='C0', label='10 meV', linewidth=3)
np.savetxt('dos_unrlx_10meV.dat', np.column_stack((energy, 2 * dos)))
dos = reconstruct_dos(moments_unrlx, energy, (a_unrlx, b_unrlx), broadening=5e-3 / a_unrlx)
ax0.plot(energy-ED, 2 * dos, c='C1', label='5 meV', linewidth=3)
np.savetxt('dos_unrlx_5meV.dat', np.column_stack((energy, 2 * dos)))
dos = reconstruct_dos(moments_unrlx, energy, (a_unrlx, b_unrlx), broadening=2e-3 / a_unrlx)
#ax0.plot(energy-ED, 2 * dos, c='C2', label='2 meV', linewidth=3)
np.savetxt('dos_unrlx_2meV.dat', np.column_stack((energy, 2 * dos)))

dos = reconstruct_dos(moments_unrlx, energy, (a_unrlx, b_unrlx), broadening=1e-3 / a_unrlx)
ax0.plot(energy-ED, 2 * dos, c='C2', label='1 meV', linewidth=3)
np.savetxt('dos_unrlx_1meV.dat', np.column_stack((energy, 2 * dos)))

ax0.plot(energy, dos_slg, c='k', ls='--', label='monolayer', linewidth=3)
np.savetxt('dos_slg.dat', np.column_stack((energy, dos_slg)))
ax0.text(0.15, 0.95, 'Unrlx \nsize {}x{}'.format(l1_unrlx, l2_unrlx), horizontalalignment='center',
         verticalalignment='center', transform=ax0.transAxes)

# Reconstruct the DOS function and plot RLX
dos = reconstruct_dos(moments_rlx, energy, (a_rlx, b_rlx), broadening=10e-3 / a_rlx)
ax1.plot(energy-ED, 2 * dos, c='C3', label='10 meV', linewidth=3)
np.savetxt('dos_rlx_10meV.dat', np.column_stack((energy, 2 * dos)))
dos = reconstruct_dos(moments_rlx, energy, (a_rlx, b_rlx), broadening=5e-3 / a_rlx)
ax1.plot(energy-ED, 2 * dos, c='C4', label='5 meV', linewidth=3)
np.savetxt('dos_rlx_5meV.dat', np.column_stack((energy, 2 * dos)))
dos = reconstruct_dos(moments_rlx, energy, (a_rlx, b_rlx), broadening=2e-3 / a_rlx)
#ax1.plot(energy-ED, 2 * dos, c='C5', label='2 meV', linewidth=3)
np.savetxt('dos_rlx_2meV.dat', np.column_stack((energy, 2 * dos)))

dos = reconstruct_dos(moments_rlx, energy, (a_rlx, b_rlx), broadening=1e-3 / a_rlx)
ax1.plot(energy-ED, 2 * dos, c='C5', label='1 meV', linewidth=3)
np.savetxt('dos_rlx_1meV.dat', np.column_stack((energy, 2 * dos)))

ax1.plot(energy, dos_slg, c='k', ls='--', label='monolayer', linewidth=3)

ax1.text(0.15, 0.95, 'Rlx \nsize {}x{}'.format(l1_rlx, l2_rlx), horizontalalignment='center',
         verticalalignment='center', transform=ax1.transAxes)

mpl.sca(ax0)
mpl.xlabel('E (eV)')
mpl.ylabel(r'${\rm states / (atom \cdot eV)}$')
mpl.sca(ax1)
mpl.xlabel('E (eV)')
ax0.set_xlim([-0.2, 0.2])
ax0.set_ylim([0, 0.04])
ax1.set_xlim([-0.2, 0.2])
ax1.set_ylim([0, 0.04])

mpl.sca(ax1)
mpl.legend()

mpl.savefig('tblg_clean_dos_broadening.png')
mpl.show()
