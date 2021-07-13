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

The search function operates through the following method:
1. Randomly produce 50,000 different sets of distances and sample the transmission of molecules
2. Sets are randomly compared against one another until all 50,000 have been compared
  a. When two sets are compared the one with the higher tramission is retained
  b. After all sets are compared only 25,000 sets remain, breed remaining sets
  c. Breed 50,000 new sets, some additional sets will be due to a mutation factor
3. Repeat step 2 until a threshold is met for tranmission of molecules
4. Repeat entire process from beginning using refined distance limits determined for initial run
