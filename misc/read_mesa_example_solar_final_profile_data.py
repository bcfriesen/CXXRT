# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np

# <codecell>

data = np.genfromtxt('/home/brian/final_profile.data', skip_header=5, names=True)

# <codecell>

print(len(data))
print('%15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s' %
      ('temperature',
       'mass_density',
       'radius',
       'h1',
       'h2',
       'he3',
       'he4',
       'li7',
       'be7',
       'b8',
       'c12',
       'c13',
       'n13',
       'n14',
       'n15',
       'o14',
       'o15',
       'o16',
       'o17',
       'o18',
       'f17',
       'f18',
       'f19',
       'ne18',
       'ne19',
       'ne20',
       'mg22',
       'mg24'))
for datum in data:
    print('%15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e %15e' %
          (10.0**datum['logT'],
           10.0**datum['logRho'],
           datum['radius'],
           datum['h1'],
           datum['h2'],
           datum['he3'],
           datum['he4'],
           datum['li7'],
           datum['be7'],
           datum['b8'],
           datum['c12'],
           datum['c13'],
           datum['n13'],
           datum['n14'],
           datum['n15'],
           datum['o14'],
           datum['o15'],
           datum['o16'],
           datum['o17'],
           datum['o18'],
           datum['f17'],
           datum['f18'],
           datum['f19'],
           datum['ne18'],
           datum['ne19'],
           datum['ne20'],
           datum['mg22'],
           datum['mg24']))

# <codecell>


