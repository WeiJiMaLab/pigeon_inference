# pigeons_inference


# Directory
Every repository should have a directory for files
describing the structure of the repository and a high-level
overview of the contents

```
pigeon_inference
│   README.md 
│
└─── analysis : all analysis files
│   │
│   └───demos : runnable demo of the experiment itself
│   │
│   └───figures : all figures included in the final paper
│   │
│   └───processed_data : results > results_N where N is the model # (see C_specifymodel_v2)
│   │
│   └───src : analysis code
│
│   
└─── experiment : all experiment files
│   │
│   └───raw_data : unmodified data -- DO NOT EDIT
│   │
│   └───src : experiment code goes here
│
└─── resources : tutorials, presentations, papers, posters, etc.

```

# Installation
You will have to install PsychToolBox to run the experiment.

Note: this code was developed using PsychToolbox version 3.0.14
It is only supported by PsychToolbox versions 3.0.14 and some older
versions.
It is not supported by v 3.0.15


# Usage (Experiment)
To run the full experiment, run pigeon_expt.m

# Usage (Analysis)
The analysis files are listed in alphabetical order from A -> F, this is the order in which they should be run.

A_readdata processes the raw data files
B_plotsummarystats produces plots of summary statistics (descriptive)
C_modelfitting is the main script that calls other model fitting files such as:
	- C_modelpredictions_v2 which produces the model predictions according to each model specification
	- C_specifymodel_v2 which has the full list of models and their different parameters, as well as their corresponding model number
	- Model fitting will result in the processed results found in processed_data/results, where results_N corresponds to the result for model number N 
D_callplotdatafits calls lower function D_plotdatawithfits which can be used to plot model fits against raw data
F_ functions are further analyses based on the log likelihood ratios found in the results_N files

## Demo(s)
To run a demo of the experiment, run demo_video.m

# Authors and Acknowledgements
Jennifer Lee and Wei Ji Ma

# License
Contact lab about licensing code
