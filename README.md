# pRad

This repository contains all you need for reconstructing a proton radiograph from the simulations performed using the GATE framework in the *gate_pCT* repository.
The output folder should have all the simulated data before reconstruction. Filenames might have changed and there are settings you can change inside the code, so you are required to make changes before running the code and producing an image. There's not much object oriented programming involved here and there's room for a lot of improvement and optimizations. The pRad reconstruction here is threefold, one is preparing the simulation data in the **createTree.cc** to create a single root file with all the variables. This root file is then filtered and performing hull-algorithm in **filterSS.cc** (Single-Sided) or **filterDS.cc** (Double-Sided) to get all the positions and proton data ready for MLP and pRad reconstruction in **unfold_Single.cc** (Single-Sided pRad) or **unfold_Double.cc** (Double-Sided pRad). Below I will list some of the settings to be changed.

REMEMBER to check all the paths to the appropriate output files!

**createTree.cc**<br />
First of all, make sure you are using either the ideal trackers or realistic trackers! Change the variable on line 38: *int detectorVersion = 1; //1 for ideal trackers, 2 for Bergen DTC*

Make sure the scanning magnet distance is the same as in the source description. There are some files you need to have prepared beforehand, these are: *PlanDescriptionToGate.txt*
and *Water_Proton.dat*. The last one is a list of stopping power values gotten directly from the simulation and used to find the wepl of the proton by measuring the initial and exiting energy of the proton. If you change the physics, you will change these values. But if you don't change the physics in the main.mac then you are fine using the attached file, also, the impact is typically negligble as long as you are using a relevant physicslist.

If you are using the realistic Bergen DTC and want to model the DTC energy reconstruction, you need the two following files (these are already there, but you might have a different file path! So check it): *weplErr.dat* for the systematic error, and *weplSigma.dat* for the noise. These are used to create a Gaussian blurring of the wepl based on the capability of the DTC.

This createTree.cc will prepare all the variables and also remove protons that did not make it out of the phantom, it creates a new root file that contains all the variables from the four tracker planes and also creates the TPS information based on the PlanDescriptionToGate file and is used later for MLP.

**filterSS.cc** and **filterDS.cc**<br />
Use filterSS.cc if you are using the single-sided setup, and filterDS.cc if you are using the conventional double-sided setup. The difference is in the extra most likely entrance step and hull algorithm needed for the single-sided setup. The double-sided setup is much simpler as you only need to do the hull-algorithm.<br />
This filtering builds all the filter distribution for use together with the 3-sigma filters during reconstruction, and also forms the new proton positions on the hull/contour of the phantom inside the root file for use in the image reconstruction.

Change the variables for the MLP algorithm if you are using a different pencil beam source and detector resolution. (line 202, the s_pos particularly needs to be changed if you are using a different pencil beam thickness!) The slighly confusing part is likely the hull-algorithm preparations, and I am only showing it for the head phantom here. The hull algorithm requires the entire head phantom (or cylinder pahantom contour) to be imported and placed in the same coordinates as the MC simulation! This way the protons can be projected from the trackers to the phantom contour and this will make the MLP estimations much better! It requires a fair bit of memory however, so be careful with running to many projections at once! If you are using a water box as the phantom, this can be simplified immensily since you know exactly where your edge is (it is 15cm away from the innermost tracker) and don't need to import anything.

You might want to change the reconstruction area for the filtering (line 90-96) if you want bigger or smaller bins/pixel to create the WEPL and angle filters, but make sure the dimensions are the same as your intended pRad! This is beacause you will use the filter bins/pixels to know which protons to remove and keep, which will affect the wepl you build inside each pixel. There is room for improvement here, but this apprach was found to work best for my purposes.

**unfold_Single.cc** and **unfold_Double.cc**<br />
This is where the reconstruction is happening. Define your MLP variables if you change anything (ideal or realistic), and your reconstruction area. This requires the filter distributions from the previous filterSS.cc to run. So this is the last step...
