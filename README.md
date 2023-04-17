# PROMETHEUS
PRObing Mass loss in Exoplanetary Transits with Hydrostatic, Evaporative and User-defined Scenarios.
PROMETHEUS is a radiative transfer tool to compute lightcurves and transmission spectra of an object
transiting its host star, typically an exoplanet. The code calculates the amount of absorption during
the transit for gaseous media in arbitrary geometry. PROMETHEUS supports various density profiles
beyond the canonical hydrostatic (barometric) law for dense atmospheres, such as the outgassed
cloud of an exomoon or a circumplanetary torus. For these tenuous exospheres, line absorption
by various atoms and ions is considered (with line lists from NIST). Additionally, it is possible
to model absorption by molecules based on ExoMOL lookup tables.

## Installation
Note that this code is written in python 3.8.3. Compatibility testing has so far been very limited.
1. Run ```git clone https://github.com/andreagebek/Prometheus.git``` in your terminal in a directory of your choice (e.g. 'exoplanets'). This will
create the 'Prometheus' base folder in the 'exoplanets' folder.
2. Create the following subfolders in the 'exoplanets' folder: setupFiles, output, figures (```mkdir setupFiles output figures```). Now,
these three folders and the 'Prometheus' folder should be on the same level.


## Usage
1. Navigate to the git subfolder (here: '/Users/agebek/exoplanets/Prometheus') and start the setup program
by typing ```python prometheus.py setup``` in your terminal in the git subfolder. This starts the Q&A session
in the terminal to setup a Prometheus calculation.
2. Answer all the questions in the Q&A session. Only answers within the specified intervals are
allowed. This Q&A session determines the density profiles of consideration (e.g. torus), 
the stellar and planetary parameters of the system (e.g. stellar radius), the abundances
of the absorbers in the system, and finally parameters related to the grids. Note that the
runtime of the code is mostly set by the choice of the 5-dimensional grid (wavelength, time,
three spatial dimensions).
3. During setup a name to store the parameter txt file has to be entered (e.g. 'testSimulation').
This file is located at '/Users/agebek/exoplanets/setupFiles/testSimulation.txt'. To run a Prometheus
calculation, navigate to the git subfolder and type ```python prometheus.py testSimulation``` in
the terminal. This will create the file '/Users/agebek/Prometheus/output/testSimulation_output.txt'
which contains the resulting (time-dependent) transit spectrum of the calculation.

## Adding more planets
Adding another exoplanet to the list of available systems is straightforward: Simply add the necessary information
about the planet and its host star to the AvailablePlanets class in the celestialBodies.py script. Also extend
the self.planetList in this class with the newly added object.

## Adding atomic/ionic absorption lines
Adding more absorption lines of atoms and ions is relatively straightforward. The absorption cross sections are
based on data of the NIST line list.
1. Add the necessary information to the AvailableSpecies class in the constants.py file.
2. Append the line list to the LineList.txt file in the Resources folder. The line list can be obtained from
the NIST database as shown in the figure. Carefully copy all the settings shown in this figure (except for 
the Spectrum field at the top where you can select the species).

![NIST example](docs/NISTexample.png?raw=true)

## Adding molecular absorption
Molecular absorption is treated differently than absorption by atoms and ions in prometheus.
The absorption cross section for molecules is absorbed based on a lookup table, where the
cross sections are tabulated based on temperature, pressure, and wavelength. These lookup
tables can be obtained from ExoMol (https://www.exomol.com/). Since these files are quite
heavy, none of them are included in this git folder by default. To add absorption features
from molecules, follow these steps:
1. Open the ExoMol website which contains the line lists for various molecules (https://www.exomol.com/data/molecules/).
2. Select the molecule and isotope of your choice.
3. Navigate to the opacity files, and download the TauREx hdf5 lookup table.
4. Create a folder named 'molecularResources' at the same level as the Prometheus base folder .
5. Copy the file into the molecularResources folder, and rename it to a name that you also use during setup (e.g. 'H2O.h5').
6. During the setup Q&A session, when asked to add absorbers, add a molecular absorber and name it according
to the lookup table name ('H2O' in this example).
