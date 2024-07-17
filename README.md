# solv
A computer program for the calculation of solvation Gibbs free energy based Based on Fixed Atomic CM5 Charges, Scaled Particle Theory, and the Atom−Atom Potential Method
## How to use
1) Download and unpack the archive solv-1.0.zip on your Linux computer:

a) click on "Source code (zip)" here https://github.com/QuantumChemistryGroup/solv/releases/tag/v1.0

b) click on "solv" here https://github.com/QuantumChemistryGroup/solv/releases/tag/v1.0

c) on your Linux machine ```unzip solv-1.0.zip```

d) put solv into the folder solv-1.0 ```mv solv solv-1.0```

c) ```cd solv-1.0```

d) ```chmod u=rwx *```

2) Prepare the file with the Cartesian coordinates and CM5 charges as described here: https://github.com/QuantumChemistryGroup/orca2cm5/tree/main
3) Run the program: ```/path/to/solv-1.0/solv your_ORCA_job.M51 YOUR_SOLVENT 500 150```

Here: 
```
your_ORCA_job.M51 is your file with Cartesian coordinates and CM5 charges
YOUR_SOLVENT is the solvent name, please have a look at file 'solv.txt' in the folder 'solv-1.0' for available solvents
500 - number of grid points per sphere for calculation of the electrostatic contribution to dGsolv
150 - number of grid points for calculation of the non-electrostatic disp-rep, the more the better, but please see the memory used
```

4) The result is printed on screen (you can use '>' symbol to print the result into file):
```
Solvent: thf                             // solvent name
Rsolv = 2.633 (unscaled)                 // unscaled solvent radius in Angstroem
param: [0.898743, 0.970106, 1.00043]     // parameters to scale the Rsolv, see the main paper
Solvent mass = 72.0575                   // solvent molecular mass in g/mol
Solvent den = 882.172                    // solvent density in kg/m**3
Rsolv(vol) = 3.18772                     // Rsolv based on solvent density (overestimated)
Gtr (1atm->x=1) = 3.37705                // conversion to 1 atm ideal gas -> solvent with solute mole fraction x=1
Warning: no R3 for: Li, assign: X        // warning if no disp-rep available
Warning: no k for: Li, assign: X
Gel = -41.52                             // dGsolv electrostatic contribution in kcal/mol
Gdisp = -26.81                           // dGsolv dispersion contribution in kcal/mol
Grep = 3.38                              // dGsolv repulsion contribution in kcal/mol
Rmol (SPT-S) = 4.83754                   // molecular radius from the molecular surface
Rmol (SPT-V) = 3.76671                   // molecular radius from the molecular volume
Gcav = 22.87 21.47 17.62                 // dGsolv cavitation energy in kcal/mol as calculated by Pierotti-Claverie, SPT-S and SPT-V
Gsolv = -42.07 -43.48 -47.33             // Final dGsolv in kcal/mol as calculated by Pierotti-Claverie, SPT-S and SPT-V (we recommend the latter)
```
> [!IMPORTANT]
> **When using this code (orca2cm5) please cite the following publications:**
> 1) Vyboishchikov, S. F.; Voityuk, A. A. Fast Non-Iterative Calculation of Solvation Energies for Water and Non-Aqueous Solvents. J. Comput. Chem. 2021, 42, 1184−1194.
> 4) Minenkov, Y. Solv: An Alternative Continuum Model Implementation Based on Fixed Atomic Charges, Scaled Particle Theory, and the Atom–Atom Potential Method. J. Chem. Theor. Comput. 2023, 19, 5221 – 5230 (DOI: 10.1021/acs.jctc.3c00410)    
