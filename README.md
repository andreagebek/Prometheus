# DISHOOM-ELITE
Radiative transfer tool to compute lightcurves of transiting exoplanets and exomoons.

## Installation
Note that this code is written in python 3.8.3. Compatibility testing has so far been very limited.
1. Create a folder (e.g. 'DishoomElite') at your desired location which will hold the git subfolder containing the code,
as well as subfolders for txt files and figures created when performing a DISHOOM-ELITE calculation.
2. Run ```git clone https://github.com/andreagebek/DishoomElite.git``` in your terminal in the 'DishoomElite' folder.
3. Create the following subfolders in the DishoomElite folder: setupFiles, output, figures (```mkdir setupFiles output figures```).

## Usage
1. Navigate to the git subfolder (here: '/Users/agebek/DishoomElite/git') and start the setup program
by typing ```python main.py``` in your terminal in the git subfolder. This starts the Q&A session
in the terminal to setup a prometheus calculation.
2. During setup a name to store the parameter txt file has to be entered (e.g. 'testSimulation').
This file is located at '/Users/agebek/DishoomElite/setupFiles/testSimulation.txt'. To run a DISHOOM-ELITE
calculation, navigate to the git subfolder and type ```python main.py testSimulation``` in
the terminal. This will create the file '/Users/agebek/DishoomElite/output/testSimulation_lightcurve.txt'
(and, depending on your setup, additional output files) which contains the resulting lightcurve of the calculation.
3. If you want to make use of the built-in plotting scripts, navigate to the plottingScripts subfolder
in the git subfolder (here: '/Users/agebek/DishoomElite/git/plottingScripts'). Run the plotting script
of your choice, e.g. 'plotSpectra.py' to plot a transmission spectrum, by typing ```python plotSpectra.py testSimulation```
in the terminal.
