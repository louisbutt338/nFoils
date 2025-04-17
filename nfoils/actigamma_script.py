import actigamma as ag # type: ignore
import matplotlib.pyplot as plt # type: ignore
import os

folder_path = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/analysis/fispact_gammaspec'
isotope_name = 'Mn56'
activity = 3.82e7

SPECTYPE = "gamma"
db = ag.Decay2012Database()

# get halflife of Co60
print(db.gethalflife(isotope_name))
# get gamma lines of Co60
print(db.getenergies(isotope_name, spectype=SPECTYPE))
print(db.getintensities(isotope_name,spectype=SPECTYPE))

# define an energy grid between 0 and 4 MeV with 5,000 bins
grid = ag.EnergyGrid(bounds=ag.linspace(0, 4e6, 10000))
# bin the lines appropriately using single type aggregator
lc = ag.LineAggregator(db, grid)
inv = ag.UnstablesInventory(data=[
    (db.getzai(isotope_name), activity),
])

hist, bin_edges = lc(inv, spectype=SPECTYPE)

plt.figure(figsize=(12,6))
plt.xlabel('Energy (eV)') 
plt.ylabel('gamma activity')
plt.tick_params(axis='y')
plt.xlim(0,3e6)
#plt.set_xscale("log")
plt.ylim(1e-5,1e15)
plt.yscale("log")
plt.stairs(hist,bin_edges)
plt.savefig(os.path.join(folder_path, '{}_gammaspec.png'.format(isotope_name)), transparent=False, bbox_inches='tight')