# pRad

This repository contains all you need for reconstructing a proton radiograph from the simulations performed using the GATE framework in the *gate_pCT* repository.
The output folder should have all the simulated data before reconstruction. Filenames might have changed and there are settings you can change inside the code, so you are required to make changes before running the code and producing an image. There's not much object oriented programming involved here and there's room for a lot of improvement and optimizations. The pRad reconstruction here is threefold, one is preparing the simulation data in the **createTree.cc** to create a single root file with all the variables. This root file is then filtered and performing hull-algorithm in **filterSS.cc** (Single-Sided) or **filterDS.cc** (Double-Sided) to get all the positions and proton data ready for MLP and pRad reconstruction in **unfold_Single.cc** (Single-Sided pRad) or **unfold_Double.cc** (Double-Sided pRad). Below I will list some of the settings to be changed.

**createTree.cc**<br />
First of all, make sure you are using either the ideal trackers or realistic trackers! Change the variable on line 38: *int detectorVersion = 1; //1 for ideal trackers, 2 for Bergen DTC*

Make sure the scanning magnet distance is the same as in the source description. There are some files you need to have prepared before hand, these are: *PlanDescriptionToGate.txt*
and *Water_Proton.dat*. The last one is a list of stopping power values gotten directly from the simulation and used to find the wepl of the proton by measuring the initial and exiting energy of the proton. If you change the physics, you will change these values. But if you don't change the physics in the main.mac then you are fine using the attached file, also, the impact is typically negligble as long as you are using a relevant physicslist.

If you are using the realistic Bergen DTC and want to model the DTC energy reconstruction, you need the two following files (these are already there, but you might have a different path! So check it): *weplErr.dat* for the systematic error, and *weplSigma.dat* for the noise. These are used to create a Gaussian blurring of the wepl based on the capability of the DTC.

REMEMBER to check the paths to the appropriate output files.

This file will prepare all the variables and also remove protons that did not make it out of the phantom, it creates a new root file that contains all the variables from the four tracker planes and also creates the TPS information from the PlanDescriptionToGate file.
