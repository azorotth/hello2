
# coding: utf-8

# Calculation of y parameter profile, through binning to mapping


# load libraries
import sys
sys.path.append('/home/tazkera/pygad/')
get_ipython().magic(u'pylab inline')
pylab.rcParams['figure.figsize'] = (12,8)
import pygad as pg
import pygad.binning


# load data
s = pg.Snap('./snap_p50n288gw_108.hdf5', load_double_prec=True)
print s

# create the pressure field
s.gas['pressure'] = s.gas['temp'].in_units_of("K") * pg.physics.kB.in_units_of('keV/K')*((s.gas['ne']*s.gas['rho'].in_units_of("g/cm**3")*0.76)/pg.physics.m_H.in_units_of("g"))


# create the y-parameter field
s.gas['y_param'] = (s.gas['pressure'].in_units_of("cm**-3 keV")*pg.physics.Th_cs.in_units_of('cm**2'))/(pg.physics.m_e.in_units_of('g')*(pg.physics.c.in_units_of('cm/s')**2))


# create a mask to exclude the star formation >0 particles
s.gas['sfr']
sfr_out = s.gas['sfr'] > 0
mask2=np.array(sfr_out)*1
mask2
s.gas['onlysfr'] = s.gas['y_param']*mask2


# bin y parameter, obtain line of sight integrated grided data
y= pygad.binning.map_qty(s.gas, field=True, qty='y_param_sfrout', extent=[[   562.95245361,  49333.65234375], [   587.37097168,  49577.0078125 ]]
                        )

# plot the binned data into a map
pg.plotting.plot_map(y.in_units_of(pg.units.Unit("h_0**-1"),{'a':1}))
savefig('y_param_sfrout.pdf')


# isolate a single galaxy
dx= 923.72

# bin y parameter only for the single galaxy
y= pygad.binning.map_qty(s.gas, field=True, qty='y_param_sfrout', extent=[[1917.85-dx, 1917.85+dx ], [15639.3-dx, 15639.3+dx]])

# plot the binned map
pg.plotting.plot_map(y.in_units_of(pg.units.Unit("h_0**-1"),{'a':1}), clogscale=True)
savefig('14M_sol.pdf')

# create profile 
prof, r_edges=pygad.binning.profile_from_map(y.in_units_of(pg.units.Unit("h_0**-1"),{'a':1}), Nbins=None,extent=[[1917.85-dx, 1917.85+dx ], [15639.3-dx, 15639.3+dx]] ,reduction='mean')

r_edges=r_edges[:len(prof)]

# plot the profile
xlabel('R')
ylabel('y')
semilogy(r_edges,prof, 'bo')
savefig('y.pdf')


# hello_world
