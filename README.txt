Created by Colin Averill and Bonnie Waring. May 24, 2017

This repository contains all necessary code to replicate model experiments and figures described in Averill and Waring 2017.
This work was done using:
R version 3.4.0 (2017-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.2 LTS

And depends upon the following package versions:
data.table version 1.10.4
rootSolve version 1.6
wesanderson version 0.3.2

The model code is contained in 'model.r'. All parameter values are contained within the file 'parameters.r', though some model experiments explicitly manipulate the input C:N parameter and the v2 (sorptive capacity) parameters. 

All experiments are in the directory 'model_experiments' 
All experiment generated data are in the directory 'experiment_output'
All scripts to generate figures from experimental output are in the directory 'figure_scripts'
All figures are saved in the directory 'figures'

So long as the directory structure is in place, and the user has the files 'parameters.r', 'model.r', and all scripts contained within the directories 'model_experiments' and 'figure_scripts', the user *should* be able to replicate all experiments and figures. However, there may be unforseen complications with version nubmers, dependencies, and other potential software incompatabilities. Who knows! Please contact the authors if you run into any trouble.