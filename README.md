# PROMETHEUS
PRObing Mass loss in Exoplanetary Transits with Hydrostatic, Evaporative and User-defined Scenarios.
PROMETHEUS is a radiative transfer tool to compute lightcurves and transmission spectra of an object
transiting its host star, typically an exoplanet. The code calculates the amount of absorption during
the transit for gaseous media in arbitrary geometry. PROMETHEUS supports various density profiles
beyond the canonical hydrostatic (barometric) law for dense atmospheres, such as the outgassed
cloud of an exomoon or a circumplanetary torus. For these tenuous exospheres, line absorption
by various atoms and ions is considered (with line lists from NIST).

## Installation
Note that this code is written in python 3.8.3. Compatibility testing has so far been very limited.
1. Run ```git clone https://github.com/andreagebek/Prometheus.git``` in your terminal in a directory of your choice (e.g. 'exoplanets'). This will
create the 'Prometheus' base folder in the 'exoplanets' folder.
2. Create the following subfolders in the 'exoplanets' folder: setupFiles, output, figures (```mkdir setupFiles output figures```). Now,
these three folders and the 'Prometheus' folder should be on the same level.


## Usage
1. Navigate to the git subfolder (here: '/Users/agebek/PrometheusProject/Prometheus') and start the setup program
by typing ```python main.py``` in your terminal in the git subfolder. This starts the Q&A session
in the terminal to setup a Prometheus calculation.
2. Answer all the questions in the Q&A session. Only answers within the specified intervals are
allowed. This Q&A session determines the density profiles of consideration (e.g. torus), 
the stellar and planetary parameters of the system (e.g. stellar radius), the abundances
of the absorbers in the system, and finally parameters related to the grids. Note that the
runtime of the code is mostly set by the choice of the 5-dimensional grid (wavelength, time,
three spatial dimensions).
3. During setup a name to store the parameter txt file has to be entered (e.g. 'testSimulation').
This file is located at '/Users/agebek/Proemtheus/setupFiles/testSimulation.txt'. To run a Prometheus
calculation, navigate to the git subfolder and type ```python main.py testSimulation``` in
the terminal. This will create the file '/Users/agebek/Prometheus/output/testSimulation_lightcurve.txt'
(and, depending on your setup, additional output files) which contains the resulting lightcurve of the calculation.
4. If you want to make use of the built-in plotting scripts, navigate to the plottingScripts subfolder
in the git subfolder (here: '/Users/agebek/Prometheus/git/plottingScripts'). Run the plotting script
of your choice, e.g. 'plotSpectra.py' to plot a transmission spectrum, by typing ```python plotSpectra.py testSimulation```
in the terminal.

## Adding atomic/ionic absorption lines
Adding more absorption lines of atoms and ions is relatively straightforward. The absorption cross sections are
based on data of the NIST line list.
1. Add the necessary information to the speciesInfoDict in the constants.py file.
2. Append the line list to the LineList.txt file in the Resources folder. The line list can be obtained from
the NIST database as shown in the figure. Carefully copy all the settings shown in this figure (except for 
the Spectrum field at the top where you can select the species).

![NIST example](docs/NISTexample.png?raw=true)