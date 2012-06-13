# Copyright (C) 2012, Christof Buchbender
# BSD License (License.txt)
import os
import site
from distutils.core import setup

if not os.path.exists(os.path.expanduser('~/.astrolyze')):
    r"""
    Creating a folder where custom setups can be stored.
    So far fized to .astrolyze in the home folder.
    """
    # TODO: Allow custom paths!!! Or NOT, decide!!!
    jskjd
    os.system('mkdir ~/.astrolyze')
    os.system('mkdir ~/.astrolyze/setup')
    os.system('mkdir ~/.astrolyze/database')
    os.system('mkdir ~/.astrolyze/lte')


def extend_path(path):
    r"""
    Creates a .pth file including all custom path for the installation.
    So far only needed for the database setup.
    """
    # TODO: Make sure the file goes to the correct folder. 
    #       So far its maybe only valid for my computer.
    string = os.path.expanduser(path)
    pth_file = open('/usr/lib/python2.7/dist-packages/astrolyze.pth', 'w')
    pth_file.write(string)
    pth_file.close()


def init_data_base():
    os.system('rm -rf ' +
              os.path.expanduser('~/.astrolyze/setup/astrolyze_prefix.py'))
    init_prefix = ('dataBase = "' + 
                    os.path.expanduser('~/.astrolyze/database/') + '"\n')
    file_out = open(os.path.expanduser('~/.astrolyze/setup/astrolyze_prefix'
                    '.py'), 'w')
    file_out.write(init_prefix)
    file_out.close()

def change_permissions():
    r"""
    This asures that the user can modify the files to customize the database
    entries and create the database itself.
    """
    os.system('chmod 777 ' +
              os.path.expanduser('~/.astrolyze/setup/line_parameter.txt'))
    os.system('chmod 777 ' +
              os.path.expanduser('~/.astrolyze/setup/galaxy_parameter.txt'))
    os.system('chmod 777 ' +
              os.path.expanduser('~/.astrolyze/setup/calibration_error.txt'))
    os.system('chmod 777 ' +
              os.path.expanduser('~/.astrolyze/database/'))

init_data_base()
extend_path('~/.astrolyze/setup/')

setup(
    name='astrolyze',
    version='0.1.0',
    author='Christof Buchbender',
    author_email='christof.buchbender@gmail.com',
    packages=['astrolyze', 
              'astrolyze/maps',
              'astrolyze/spectra',
              'astrolyze/sed', 
              'astrolyze/lte', 
              'astrolyze/functions', 
             ],
    scripts=['bin/setup_astrolyze.py'],
    data_files=[(os.path.expanduser('~/.astrolyze/setup/'),
                ['cfg/line_parameter.txt', 'cfg/galaxy_parameter.txt',
                'cfg/calibration_error.txt'])],
    url='http://www.strange-associations.de/astrolyze/',
    license='LICENSE.txt',
    description='Reduction and analysing tools for (mainly) Radioastronomy.',
    long_description=open('README.txt').read(),
    classifiers=['Topic :: Scientific/Engineering :: Astronomy'],
    requires=[
        "numpy",
        "pyfits",
        "matplotlib",
        "scipy",
        "pywcs",
        "copy",
        "random",
        "pysqlite2",
        "socket"
        "docutils"
    ],
)

change_permissions()
