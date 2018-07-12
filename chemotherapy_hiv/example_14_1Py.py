from forward_backward_sweep import ForwardBackwardSweep
import matplotlib.pyplot as plt
import matplotlib as mpl

s = 10.0
m_1 = 0.02
m_2 = 0.5
m_3 = 4.4
r = 0.03
k = 0.000024
n_weight = 300.0
a = 0.2
t_cell_max = 1500
t_cell_zero = 806.4
t_cell_infected_zero = 0.04
virus_particle_zero = 1.5

fbsm = ForwardBackwardSweep()
fbsm.set_parameters(s, m_1, m_2, m_3, r, k, n_weight, a,
                       t_cell_max, t_cell_zero, t_cell_infected_zero,
                       virus_particle_zero)
t = fbsm.t

x_wc = fbsm.runge_kutta_forward(fbsm.u)
[x, lambda_, u] = fbsm.forward_backward_sweep()

mpl.style.use('ggplot')
# plt.ion()
plt.show()

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0))
ax4 = plt.subplot2grid((2, 2), (1, 1))

'''
ax1.plot(t, x_wc[:, 0], '-',
         ms=3,
         lw=2,
         alpha=1,
         color='darkgreen',
         label='Uncontrolled T cells'
         )
'''
ax1.plot(t, x[:, 0], '--',
         ms=3,
         lw=2,
         alpha=1.0,
         color='orange',
         label='Cotrolled T cells')
ax1.set_ylabel(r'Suceptibles $T_{cells}$')
ax1.set_xlabel(r'Time (days)')
#ax1.legend(loc=0)
#
'''
ax2.semilogy(t, x_wc[:, 1], '-',
             ms=3,
             lw=2,
             alpha=1.0,
             color='darkgreen',
             label=r'Uncontrolled $T_{cells}^{infected}$')

'''
ax2.semilogy(t, x[:, 1], '--',
             ms=3,
             lw=2,
             alpha=1.0,
             color='orange',
             label=r'Controlled $T_{cells}^{infected}$')
ax2.set_ylabel(r'Infected $T_{cells}$')
ax2.set_xlabel(r'Time(days)')
#
'''
ax3.semilogy(t, x_wc[:, 2], '-',
             ms=3,
             lw=2,
             alpha=1.0,
             color='darkgreen',
             label=r'Uncontrolled $Virus$')
'''
ax3.semilogy(t, x[:, 2], '--',
             ms=3,
             lw=2,
             alpha=1.0,
             color='orange',
             label=r'Controlled $Virus')
ax3.set_ylabel(r' $Virus$')
ax3.set_xlabel(r'Time(days)')
ax4.plot(t, u, '--',
         ms=3,
         lw=2,
         alpha=1.0,
         color='orange',
         label='Control')
ax4.set_ylabel('Chemotherapy  $u(t)$')
ax4.set_xlabel('Time (days)')
#
# plt.legend(loc=0)
plt.tight_layout()
plt.show()
