#!/bin/csh -fx

#To run:
#./findDQE.csh > outputDQE.txt &

/labdata/allab/michaelc/FindDQE/FindDQE/bin/finddqe.exe << eof
2	#Option #1: Input counts/electron or #2: Input electrons/pixel
49.9	#Electrons/pixel = (PixelSize^2)*(Dose e-/A2)
1	#Multiplying factor
1	#Number of input images 
14jan09a_00002en.mrc	#Micrograph
14jan09a_00002en.residual.mrc	#Residual output file
eof 
