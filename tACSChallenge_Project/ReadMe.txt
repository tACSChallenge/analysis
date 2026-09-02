tACS Challenge Analysis ReadMe
Written by Anne-Fiona Griesfeller

1. Required Matlab version and toolboxes:

* MATLAB R2020b or later
* Signal Processing Toolbox: 			https://fr.mathworks.com/products/signal.html
* Statistics and Machine Learning Toolbox:	https://fr.mathworks.com/products/statistics.html 
* Circular Statistics Toolbox (CircStat):	https://fr.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

1.1 For optional Bayesian analysis R with: 	BayesFactor package: Morey, R. D., & Rouder, J. N. (2024). BayesFactor: Computation of Bayes factors for common										designs. R package version 0.9.12-4.7.


2.1 How to run:

All analyses are run with Run_All.m in tACSChallenge_Project/scripts/

Within this script:

* Enable or disable H1-H3
* Enable or disable H4-H6, and single hypotheses within them
* Select which labs to include
* Set the number of permutations (standard 1000)
* Force a rebuild that ignores PhaseMatrics and AnalysisTable files that are already in the lab folders. 

No other script needs to be run manually.


!H1-H3 must be run at least once before H4-H6, because H4-H6 read the PhaseMetrics.mat file that H1-H3 writes into each lab folder.

!Behavioral files and Metadata must follow SOP 5.


By default lab_ids is empty, which analyses every lab folder found in data/.
Set it to e.g. lab_ids = [3 45], to restrict the analysis.


The project root is taken from the location of Run_All.m. 
No paths need editing as long as the folder structure is followed (see 2.2 below).


2.2 Folder structure:

tACSChallenge_Project/
|
|-- data/
|   |-- L01/
|   |-- L02/
|   |-- L18/
|   |   |-- sub-L18_S01/
|   |   |   |-- beh/
|   |   |   |-- metadata/
|   |   |-- sub-L18_S02/
|   |   |-- PhaseMetrics.mat    (created automatically by H1-H3)
|   |   |-- AnalysisTable.mat   (created automatically by H4-H6)
|   |-- ...
|
|-- functions/
|   |-- tACSChallenge_ImportData.m         	read one block file
|   |-- tACSChallenge_SortData.m           	trials and tACS phase per participant
|   |-- tACSChallenge_EvalData.m           	phase regression and permutations
|   |-- tACSChallenge_buildAnalysisTable.m 	one row per participant per condition
|   |-- tACSChallenge_runH1H2H3_Pooled.m   	H1 to H3
|   |-- tACSChallenge_runH4H5H6.m          	wrapper for H4 to H6
|   |-- tACSChallenge_runGlobalModels.m    	H4, H5, H6 LMMs
|   |-- tsc_parse_beh_filename.m           	parse lab, subject, block, condition
|   |-- tsc_perf_from_sorted.m             	hits, false alarms, d', criterion, RT
|   |-- tACS_bayesFactorBIC.m              	Bayes factors
|   |-- tACS_condColours.m                 	shared figure palette for conditions
|   |-- tACS_raincloudPanel.m              	raincloud panels (H2, H5, H6)
|   |-- tACS_scatterPanel.m                	scatter panels (H4)
|   |-- circstat/                          	download and add to this folder 
|
|
|-- bayes/       				Bayes factors for H1 and H2, see bayes/ReadMe.txt
|
|-- results/
|   |-- H1H2H3/
|   |-- H4H5H6/
|
|-- scripts/
    |-- Run_All.m



3. Output is written to the results folder, by hypothesis.
 
results/H1H2H3/pooled/
|-- H1/   figures, H1_results.txt
|-- H2/   figures, H2_results.txt
|-- H3/   H3_results.txt
|-- H1H2H3_pooled_metrics.mat
 
results/H4H5H6/
|-- H4/   figures, model report
|-- H5/   figures, model report
|-- H6/   figures, model report
|-- AnalysisTable_allLabs.mat / .csv



4. Bayesian Analyses
 
Bayes factors are preregistered for every hypothesis, but the protocol does
not specify how they are computed.

See bayes_ReadMe.txt in /bayes for details on what was done here. 
