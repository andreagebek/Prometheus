# PROMETHEUS
PRObing Mass loss in Exoplanetary Transits with Hydrostatic, Evaporative and User-defined Scenarios.
PROMETHEUS is a radiative transfer tool to compute lightcurves and transmission spectra of an object
transiting its host star, typically an exoplanet. The code calculates the amount of an absorption during
the transit for gaseous media in arbitrary geometry.

## Installation
Note that this code is written in python 3.8.3. Compatibility testing has so far been very limited.
1. Create a folder (e.g. 'Prometheus') at your desired location which will hold the git subfolder containing the code,
as well as subfolders for txt files and figures created when performing a PROMETHEUS calculation.
2. Run ```git clone https://github.com/andreagebek/Prometheus.git``` in your terminal in the 'Prometheus' folder.
3. Create the following subfolders in the Prometheus folder: setupFiles, output, figures (```mkdir setupFiles output figures```).
4. Optional: If you want to include molecular line absorption (currently, only sulfur dioxide is supported), create an additional 
subfolder in the Prometheus folder named molceularResources: ```mkdir molecularResources```. Download the sulfur dioxide line list
in TauRex format from the ExoMOL database under https://www.exomol.com/data/molecules/SO2/32S-16O2/ExoAmes/ and store it as an hdf5 
file ('.h5') in the molecularResources folder.

## Usage
1. Navigate to the git subfolder (here: '/Users/agebek/Prometheus/git') and start the setup program
by typing ```python main.py``` in your terminal in the git subfolder. This starts the Q&A session
in the terminal to setup a Prometheus calculation.
2. During setup a name to store the parameter txt file has to be entered (e.g. 'testSimulation').
This file is located at '/Users/agebek/Proemtheus/setupFiles/testSimulation.txt'. To run a Prometheus
calculation, navigate to the git subfolder and type ```python main.py testSimulation``` in
the terminal. This will create the file '/Users/agebek/Prometheus/output/testSimulation_lightcurve.txt'
(and, depending on your setup, additional output files) which contains the resulting lightcurve of the calculation.
3. If you want to make use of the built-in plotting scripts, navigate to the plottingScripts subfolder
in the git subfolder (here: '/Users/agebek/Prometheus/git/plottingScripts'). Run the plotting script
of your choice, e.g. 'plotSpectra.py' to plot a transmission spectrum, by typing ```python plotSpectra.py testSimulation```
in the terminal.
