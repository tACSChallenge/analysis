tACS Challenge Analysis README
Written by Anne-Fiona Griesfeller, January 26

Overview

*How to Run the Analyses*

All analyses are run via Run_all.m located in tACSChallenge_Project/scripts/

Within this script you can:

* Enable or disable H1–H3
* Enable or disable H4–H6
* Select which labs to include
* Run the full analysis

No other scripts need to be run manually.

Important: Run_all.m assumes a specific folder structure and separates:

⦁	Data storage from each Lab
⦁	Analysis functions
⦁	Results output
⦁	Main script

Please follow the required structure below to ensure that the analyses run correctly and remain comparable across labs.

*Required Folder Structure*

tACSChallenge_Project/
│
├── data/
│   ├── L01/
│   ├── L02/
│   ├── L18/
│   │   ├── sub-L18_S01/
│   │   │   ├── beh/
│   │   │   └── metadata/
│   │   ├── sub-L18_S02/
│   │   └── PhaseMetrics.mat   (created automatically after running H1-3)
│   └── ...
│
├── functions/
│   ├── tACSChallenge_ImportData.m
│   ├── tACSChallenge_SortData.m
│   ├── tACSChallenge_EvalData.m
│   ├── tACSChallenge_AnalyseData.m
│   ├── tACSChallenge_buildAnalysisTable.m
│   ├── tACSChallenge_runGlobalModels.m
│   ├── tACSChallenge_runH1H2H3_pooled.m
│   ├── tACSChallenge_runH4H5H6.m
│   ├── helper functions (*.m)
│   └── circular_statistics_toolbox/
│
├── results/
│   ├── H1H2H3/
│   └── H4H5H6/
│
└── scripts/
    └── Run_All.m
---


Otherwise important:

The file `PhaseMetrics.mat` is created automatically when running H1–H3
It is stored inside each lab folder
This file is required for Hypothesis 4 and must not be moved manually

---

Where to find results

All outputs are written to the results folder (`results/`).

H1–H3 results will be pooled across all labs, no per-lab result is saved by default. 
Subfolders for each hypothesis will be created upon running H1-3, in which figures and statistics can be found.

results/H1H2H3/
├── H1/
│   ├── figures
│   └── H1_results.txt
├── H2/
│   ├── figures
│   └── H2_results.txt
└── H3/
    └── H3_results.txt
```


H4–H6 Results
Results are saved per hypothesis as well, and output will include mixed effects models and figures.
```
results/H4H5H6/
├── H4/
│   ├── figures
│   └── model reports
├── H5/
│   ├── figures
│   └── model reports
└── H6/
    ├── figures
    └── model reports
```
