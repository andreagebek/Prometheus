# prometheus
Radiative transfer tool to compute transmission spectra of transiting exoplanets and exomoons.

## Installation
Note that this code is written in python 3.8.3. Compatibility testing has so far been very limited.
1. Create a folder (e.g. 'prometheus') at your desired location which will hold the git subfolder containing the code,
as well as all txt files and figures created when performing a prometheus calculation.
2. Run 'git clone https://github.com/andreagebek/prometheus.git' in your terminal in the 'prometheus' folder.

## Usage
1. Navigate to the git subfolder (here: '/Users/agebek/prometheus/git') and start the setup program
by typing 'python setup.py' in your terminal in the git subfolder. This starts the Q&A session
in the terminal to setup a prometheus calculation.
2. During setup a name to store the parameter txt file has to be entered (e.g. 'testSimulation').
This file is located at '/Users/agebek/prometheus/testSimulation.txt'. To run a prometheus
calculation, navigate to the git subfolder and type 'python main.py testSimulation' in
the terminal. This will create the file '/Users/agebek/prometheus/testSimulation_spectrum.txt'
which contains the resulting spectrum of the calculation.
3. If you want to make use of the built-in plotting scripts, navigate to the plottingScripts subfolder
in the git subfolder (here: '/Users/agebek/prometheus/git/plottingScripts'). Run the plotting script
of your choice, e.g. 'spectrum.py' to plot a transmission spectrum, by typing 'python spectrum.py testSimulation'
in the terminal.
