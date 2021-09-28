# Multipole_Trajectory
Used to calculate the trajectory of molecules through a multipole or multipoles.

Has the following requirements:
- python 3.6+
- pip install numpy
- pip install matplotlib

The search functions utilizes multiple cores to accelerate discovery of parameters, it will limit
the number of cores to the number available on the machine. There is no memory limit and I have not
tested the maximum amount of ram required for operation, however, I imagine 8gb of ram will be more
than sufficient.

Additional_Information.txt outlines what parameters mean and how certain functions were determined.
Search_Information.txt outlines how to setup a search for optimal distances.
