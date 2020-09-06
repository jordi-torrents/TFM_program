This is a Master final project program simulating the Vicsek model and two new Vicsek variation using levy flights and burst-and-coasts dynamics.

The program can run three different models depending on the input mode:
- **Mode 0**: Classic Vicsek model
- **Mode 1**: Lévy behaviour model
- **Mode 2**: Burst-and-coast model

And different measures and options can be done in the input file. This file includes the following parameters:

 - **Output directory**: directory where to write all output data
 - **Status filename**: filename of the 'status' file (if None, an automatic name is set)
 - **Mode**: desired model (0/1/2)
 - **Noise Scan mode**: take noise values list from the file 'noise_values.dat' (True/False)
 - **Noise value**: desired noise valu if Scan_mode is False
 - **Lévy alpha parameter**: characteristic power-law parameter for the Lévy model
 - **Mean flight**: mean flight parameter for Vicsek Model and Lévy model
 - **burst amplitude**: Delta-v for the burst-and-coast model
 - **N_measures**: number of measurements in every iteration
 - **N_steps**: time-steps between consecutives measurements
 - **N_reset**: time-steps to reset the system (thermalazing)
 - **N_iterations**: number of total iteration of the simulation
 - **N_GNF**: the number of boxes to use in the GNF computation
 - **Polarization**: temporal evolution of the global order (True/False)
 - **GNF**: Giant Nuber Fluctuations measurement (True/False)
 - **Configuration**: Positions of all particles through the simulation to visualize them later (True/False)
 - **Speed**: temporal evolution of the global group speed (True/False)
 - **ignore_last_config**: do not use last scenario from last simulation to start the new one (True/False)
 - **L**: desired system size
 - **Density**: desired system density
 - **Random seed**

Python file presents a fast analysis of data in one output directory. To use the visualization in video use the file visualyzer.py changing the proper directory in its 7th line. 

- **Name**: Jordi Torrents Monegal
- **Email**: jordi.torrentsm@gmail.com
