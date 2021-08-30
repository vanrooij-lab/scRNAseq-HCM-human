#!/bin/bash

################################################################################

# Note: I got MARGE running in the end, but it yielded empty output files;
# note sure why. In any case, I've switched to using LISA, which seems 
# to be a follow-up of this tool.

# Note2: when installing LISA, i might have accidentally broken
# the MARGE conda env.

################################################################################

# To install MARGE see notes in Word file also,
if [[ "" == "install" ]]; then
  conda create -n MARGE
  conda activate MARGE
  conda install -c bioconda snakemake
  conda install -c conda-forge argparse
  conda install numpy
  # pip install -U scikit-learn
  pip install scikit-learn
  
  conda install -c anaconda pytables
  
  conda install -c bioconda macs2
  conda install scipy=1.5.2
  
  conda install matplotlib
  
  # pip install scikit-learn==0.18
  
  # conda install -c bioconda marge 
  # --> doesn't seem to work on HPC
  # Did seem to work locally
  # which MARGE
  # /opt/anaconda3/envs/MARGE/bin/MARGE
  
  # Followed their install
  # Added stuff to paths
  
  conda install -c bioconda ucsc-bedclip
  conda install -c bioconda ucsc-bigwigtobedgraph
  conda install -c bioconda ucsc-bigWigSummary
  conda install -c bioconda ucsc-bigWigAverageOverBed
  
  # Another library that appears necessary
  conda install -c conda-forge regions
  
  # http://cistrome.org/MARGE/src/Py3_MARGE.zip
  cd /hpc/hub_oudenaarden/mwehrens/bin/MARGE/
  curl -o Py3_MARGE.zip http://cistrome.org/MARGE/src/Py3_MARGE.zip
  python setup.py install --prefix=/hpc/hub_oudenaarden/mwehrens/bin/MARGE/installed

  # To address marge import issue
  export PYTHONPATH=$PYTHONPATH:/hpc/hub_oudenaarden/mwehrens/bin/MARGE/Py3_MARGE
  # Now make marge itself findable
  export PATH=$PATH:/hpc/hub_oudenaarden/mwehrens/bin/MARGE/Py3_MARGE/marge
    # Could be added to bashrc 
    
  marge
  
  # Manual fixes to /hpc/hub_oudenaarden/mwehrens/bin/MARGE/Py3_MARGE/marge/rpregress.py
  # where required..
  
  # #from sklearn.grid_search import GridSearchCV
  # from sklearn.model_selection import GridSearchCV # MW
  # # from sklearn import cross_validation
  # from sklearn import model_selection as cross_validation # MW

  # Two instances:
  # #LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01)
  # LR_l1 = linear_model.LogisticRegression(penalty='l1', tol=0.01, solver='liblinear') # MW
  
  # Also in other files
  
  # (..) to do: also mention those changes (again libs and linear_model)
  
fi

# Also download the ref files @ http://cistrome.org/MARGE/download.html
# I downloaded "hg38 Library (all)"

# To execute MARGE see:
# http://cistrome.org/MARGE/tutorial.html

# Configure the config.json as desired, see example @ http://cistrome.org/MARGE/download.html
# I choose "config.json for GeneList only"

# Note that some of the settings are a bit non-intuitive
# - MARGEdir, when using conda, should point to /opt/anaconda3/pkgs/marge-1.0-py_2/site-packages/marge
# AFTER REINSTALLING --> /Users/m.wehrens/Software_custom/MARGE/marge
# - REFdir should point to the directory with the downloaded refs

################################################################################
# Running locally

# cd /Users/m.wehrens/Data/_2019_02_HCM_SCS/2021_HPC_analysis.3/MARGE/test_run

# conda activate MARGE
# snakemake -n
# snakemake --cores 1

################################################################################
# HPC

# To execute this script @ HPC, use 
# cd /hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/MARGE
# subdir=XXX
# sbatch --job-name=MARGE -c 1 --time=1-00:00:00 --mem=32G --export=ALL,subdir="${subdir}" --output=slurm-${subdir}-%x.%j.out /hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/HCM_SCS_2021_08_MARGE.sh
# sbatch --job-name=MARGE -c 2 --time=1-00:00:00 --mem=32G --export=ALL,subdir="${subdir}" --output=slurm-${subdir}-%x.%j.out /hpc/hub_oudenaarden/mwehrens/scripts/SCS_HCM_analysis/HCM_SCS_2021_08_MARGE.sh

cd /hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/MARGE/${subdir}

# First make conda accessible in this shell (https://github.com/conda/conda/issues/7980)
condapath=$(conda info --base)
source ${condapath}/etc/profile.d/conda.sh
# Activate MARGE environment
conda activate MARGE

# Export pahts
export PYTHONPATH=$PYTHONPATH:/hpc/hub_oudenaarden/mwehrens/bin/MARGE/Py3_MARGE
export PATH=$PATH:/hpc/hub_oudenaarden/mwehrens/bin/MARGE/Py3_MARGE/marge

# Only do this when you're setting up the job
# Will create snakemake thing and config.json
if [[ "init" == "" ]]; then 
  marge init /hpc/hub_oudenaarden/mwehrens/data/HCM_SCS_RHL.3/MARGE/${subdir}
fi

# Actually start the job
snakemake -n
snakemake --cores 20











