from distutils.core import setup

setup(
    name='AstroLyze',
    version='0.1.0',
    author='Christof Buchbender',
    author_email='christof.buchbender@gmail.com',
    packages=['astrolyze', 'astrolyze/maps', 'astrolyze/spectra','astrolyze/sed',
              'astrolyze/maptools', 'astrolyze/lte', 'astrolyze/functions',
              'astrolyze/setup', 'astrolyze/setup/paths', 'astrolyze/setup/database'],
    scripts=[],
    url='',
    license='LICENSE.txt',
    description='Reduction and analysing tools for (mainly) Radioastronomy.',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.1.1",
        "pyfits >= 0.1.4",
        "matplotlib >= 0.1.0",
        "scipy >= 0.1.0",
        "pywcs >= 0.1.0",
        "copy >= 0.1.0",
        "random >= 0.1.0",
        "pysqlite2 >= 0.1.0",
        "socket >= 0.1.0"
    ],
)

