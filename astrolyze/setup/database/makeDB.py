import os
from astrolyze.setup.paths import prefix

os.system('rm -rf '+prefix.dataBase+'parameter.db')
os.system('python lineParamDB.py')
os.system('python galaxyParamDB.py')
os.system('python mapParamDB.py')
