import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.pyplot as mpl

# ************************************
# keep these definitions for kite website
import seaborn as sns

mpl.rcParams['figure.dpi'] = 640
mpl.rcParams['savefig.dpi'] = 640
sns.set_style("white")
# Kite color scheme
colors = ["dusty purple", "faded green", "windows blue", "amber", "greyish"]
current_palette = sns.xkcd_palette(colors)
sns.set_palette(current_palette)
sns.set_style("ticks")
# sns.set_context("talk", font_scale=1.2)

mpl.rc('text', usetex=True)
mpl.rc('text.latex', unicode=True)
font = {'family': 'serif',
        'weight': 'bold',
        'size': 28}
mpl.rc('font', **font)

# ************************************
# read walltime STRIDE, nCores
nCores = [64, 32, 32, 16, 8, 8, 4]
stride64 = [14.75454545, 9.117977528, 7.995073892, 4.293650794, 2.217213115, 1.887209302, 1]
stride128 = [13.15966387, 8.15625, 7.085972851, 3.857142857, 1.982278481, 1.921472393, 1]
stride256 = [10.58823529, 7.297297297, 6.923076923, 3.941605839, 2.040302267, 1.864211738, 1]

stride64_rel_scaled = [12.94545455, 8, 7.014778325, 3.767195767, 1.945355191, 1.655813953, 0.8773875539]
stride128_rel_scaled = [12.90756303, 8, 6.950226244, 3.783251232, 1.944303797, 1.884662577, 0.9808429119]
stride256_rel_scaled = [10.37908497, 7.153153153, 6.786324786, 3.863746959, 2, 1.827387802, 0.9802469136]

perfect_scale_nc = [16, 8, 8, 4, 2, 2, 1]

np.savetxt('ncores_stride64_wlt_scaled.dat', np.column_stack((nCores, stride64)))
np.savetxt('ncores_stride128_wlt_scaled.dat', np.column_stack((nCores, stride128)))
np.savetxt('ncores_stride256_wlt_scaled.dat', np.column_stack((nCores, stride256)))

np.savetxt('ncores_stride64_wlt_rel_scaled.dat', np.column_stack((nCores, stride64_rel_scaled)))
np.savetxt('ncores_stride128_wlt_rel_scaled.dat', np.column_stack((nCores, stride128_rel_scaled)))
np.savetxt('ncores_stride256_wlt_rel_scaled.dat', np.column_stack((nCores, stride256_rel_scaled)))
# ************************************
# read walltime SIZE Walltime
lxly = [8192 * 8192, 8192 * 16384, 16384 * 8192, 16384 * 16384, 16384 * 32768, 32768 * 16384, 32768 * 32768,
        32768 * 65536, 65536 * 32768, 65536 * 65536]

size64wlt = [1, 2.032258065, 1.838709677, 3.935483871, 7.677419355, 8.290322581, 16.64516129, 33.67741935, 37.77419355,
             74.70967742]
size64mem = [3.103248596, 6.119796753, 6.118953705, 12.13285828, 24.16428757, 24.16258621, 48.18500137, 96.2477684,
             96.24792862, 192.2946587]

size64wlt_rel_scaled = [0.8471391973, 1.721605465, 1.55764304, 3.333902647, 6.503842869, 7.023057216, 14.10076857,
                        28.529462, 32, 63.28949616]

perfect_scale_size = [1, 2, 2, 4, 8, 8, 16, 32, 32, 64]

np.savetxt('lxly_wlt_scaled.dat', np.column_stack((lxly, size64wlt)))
np.savetxt('lxly_mem_scaled.dat', np.column_stack((lxly, size64mem)))

np.savetxt('lxly_wlt_rel_scaled.dat', np.column_stack((lxly, size64wlt_rel_scaled)))
np.savetxt('lxly_mem_rel_scaled.dat', np.column_stack((lxly, size64mem)))

# ************************************
# read walltime nMoments
nM = np.arange(1000, 11000, 1000)
nMoments64 = [1, 2.008695652, 2.947826087, 3.965217391, 5.095652174, 6.086956522, 6.965217391, 8.060869565, 8.965217391,
              10.2]
nMoments64_rel_scaled = [0.9803921569, 1.969309463, 2.890025575, 3.887468031, 4.995737425, 5.967604433, 6.828644501,
                         7.902813299, 8.789428815, 10]
perfect_scale_nm = np.arange(1, 11, 1)
np.savetxt('nM_wlt_scaled.dat', np.column_stack((nM, nMoments64)))

np.savetxt('nM_wlt_rel_scaled.dat', np.column_stack((nM, nMoments64_rel_scaled)))

# ********************************************* #
# PLOTING WITH RESPECT TO THE INITIAL REFERENCE #
# ********************************************* #
# define figure
fig = mpl.figure(figsize=(16, 10))
# fig.suptitle('Benchmarking DOS KITE')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.07, right=0.95, bottom=0.1, top=0.98, wspace=0.2, hspace=0.2)
lambda_block = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[:, :], wspace=0.19, hspace=0.3)

ax0 = mpl.subplot(lambda_block[0])
ax1 = mpl.subplot(lambda_block[1])
ax2 = mpl.subplot(lambda_block[2])
ax3 = mpl.subplot(lambda_block[3])

ax0.scatter(nCores, stride64, s=180, label='KITE STRIDE=64', marker='*', alpha=0.8)
ax0.scatter(nCores, stride128, s=180, label='KITE STRIDE=128', marker='s', alpha=0.8)
ax0.scatter(nCores, stride256, s=180, label='KITE STRIDE=256', marker='d', alpha=0.8)
ax0.plot(nCores, perfect_scale_nc, label='ideal scaling', linestyle='--', linewidth=3, color='grey')
ax0.legend(loc=4, prop={'size': 18})
ax0.text(0.35, 0.9, 'Scaling STRIDE/nCores', horizontalalignment='center',
         verticalalignment='center', transform=ax0.transAxes)
ax0.set_ylabel('Speedup')
ax0.set_xlabel('nCores')

ax1.scatter(lxly, size64wlt, label='KITE', s=180, marker='*', alpha=0.8)
ax1.plot(lxly, perfect_scale_size, label='ideal scaling', linestyle='--', linewidth=3, color='grey')
ax1.legend(loc=4, prop={'size': 18})
ax1.text(0.35, 0.9, 'Scaling SIZE (Runtime)', horizontalalignment='center',
         verticalalignment='center', transform=ax1.transAxes)
ax1.set_ylabel('Relative Eff.')
ax1.set_xlabel('Effective size ($L_x \\times L_y$)')

ax2.scatter(nM, nMoments64, label='KITE', s=180, marker='*', alpha=0.8)
ax2.plot(nM, perfect_scale_nm, linestyle='--', linewidth=3, label='ideal scaling', color='grey')
ax2.legend(loc=4, prop={'size': 18})
ax2.text(0.3, 0.9, 'Scaling nMoments', horizontalalignment='center',
         verticalalignment='center', transform=ax2.transAxes)
ax2.set_ylabel('Relative Eff.')
ax2.set_xlabel('nMoments')

ax3.scatter(lxly, size64mem, s=180, marker='*', alpha=0.8)
ax3.set_ylabel('Memory (GB)')
ax3.set_xlabel('Effective size ($L_x \\times L_y$)')
ax3.text(0.35, 0.9, 'Scaling SIZE (Memory)', horizontalalignment='center',
         verticalalignment='center', transform=ax3.transAxes)

mpl.savefig('benchmarks.svg', format='svg')
mpl.show()

# ********************************************* #
# PLOTING WITH RESPECT TO THE SCALED REFERENCE #
# ********************************************* #

# define figure
fig = mpl.figure(figsize=(16, 10))
# fig.suptitle('Benchmarking DOS KITE')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.07, right=0.95, bottom=0.1, top=0.98, wspace=0.2, hspace=0.2)
lambda_block = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[:, :], wspace=0.19, hspace=0.3)

ax0 = mpl.subplot(lambda_block[0])
ax1 = mpl.subplot(lambda_block[1])
ax2 = mpl.subplot(lambda_block[2])
ax3 = mpl.subplot(lambda_block[3])

ax0.scatter(nCores, stride64_rel_scaled, s=180, label='KITE STRIDE=64', marker='*', alpha=0.8)
ax0.scatter(nCores, stride128_rel_scaled, s=180, label='KITE STRIDE=128', marker='s', alpha=0.8)
ax0.scatter(nCores, stride256_rel_scaled, s=180, label='KITE STRIDE=256', marker='d', alpha=0.8)
ax0.plot(nCores, perfect_scale_nc, label='ideal scaling', linestyle='--', linewidth=3, color='grey')
ax0.legend(loc=4, prop={'size': 18})
ax0.text(0.35, 0.9, 'Scaling STRIDE/nCores', horizontalalignment='center',
         verticalalignment='center', transform=ax0.transAxes)
ax0.set_ylabel('Speedup')
ax0.set_xlabel('nCores')

ax1.scatter(lxly, size64wlt_rel_scaled, s=180, label='KITE', marker='*', alpha=0.8)
ax1.plot(lxly, perfect_scale_size, label='ideal scaling', linestyle='--', linewidth=3, color='grey')
ax1.legend(loc=4, prop={'size': 18})
ax1.text(0.35, 0.9, 'Scaling SIZE (Runtime)', horizontalalignment='center',
         verticalalignment='center', transform=ax1.transAxes)
ax1.set_ylabel('Relative Eff.')
ax1.set_xlabel('Effective size ($L_x \\times L_y$)')

ax2.scatter(nM, nMoments64_rel_scaled, s=180, label='KITE', marker='*', alpha=0.8)
ax2.plot(nM, perfect_scale_nm, linestyle='--', linewidth=3, label='ideal scaling', color='grey')
ax2.legend(loc=4, prop={'size': 18})
ax2.text(0.3, 0.9, 'Scaling nMoments', horizontalalignment='center',
         verticalalignment='center', transform=ax2.transAxes)
ax2.set_ylabel('Relative Eff.')
ax2.set_xlabel('nMoments')

ax3.scatter(lxly, size64mem, s=180, marker='*', alpha=0.8)
ax3.set_ylabel('Memory (GB)')
ax3.set_xlabel('Effective size ($L_x \\times L_y$)')
ax3.text(0.35, 0.9, 'Scaling SIZE (Memory)', horizontalalignment='center',
         verticalalignment='center', transform=ax3.transAxes)

mpl.savefig('benchmarks_rel_scale.svg', format='svg')
mpl.show()
