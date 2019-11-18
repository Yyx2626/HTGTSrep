Baseline_Functions.r
--------------------
Contains all the functions needed for the Baseline_Main.r and SHMulation.r



Baseline_Main.r
---------------

The main program that drives BASELINe. You can run the program from the command line using RScript:

e.g.
Rscript Baseline_Main.r [testID] [species] [substitutionModel] [mutabilityModel] [clonal] [fixIndels] [region] [inputFilePath] [outputPath] [ouputID]

* testID: 
	1 = Focused test
	2 = Local test
	The statistical framework used to test for selection. 
	Both the Focused & Local statitics are decribed in Hershberg U et al. (2008) & 
	Uduman M et al. (2011).
	
*	species:
	1 = Human
	2 = Mouse
	The choice of species determines the intrinsic targeting model for somatic 
	hypermutation along with the substitution matrix. This information is used 
	along with the input germline sequence, to calculate the expected pattern of R 
	and S mutations in the absence of selection.

* substitutionMode:
	0 = Uniform substitution
	1 = Model adapated from Smith DS et al. 1996

*	mutabilityModel:
	0 = Uniform mutability
	1 = Model adapated from Shapiro GS et al. 2002

* clonal:
	0 = Sequences are all independent
	1 = Treat each germline group as a single, independent clone.

*	fixIndels:
	0 = Ignore all sequences with insertion/deletions
	1 = Attempt to correct the indels by adding an "N" to the germline (insertion) or the input (deletion), at the position of the indel.

*	region:		
	Definition of the CDR/FWR boundaries (Amino Acid residues)

*	inputFilePath:
	Full path to the input file

* outputPath:
	Full path to output directory

* outputID:
	ID to name output files
	


SHMulation.r
------------

This program simulates somatic hypermutation & selection	
