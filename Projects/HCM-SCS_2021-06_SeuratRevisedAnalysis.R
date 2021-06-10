

# Revised version of earlier script to analyze Van Rooij single cell data of
# HCM tissue samples from myectomies.
#
# Previously, we used RaceID2 to analyze data (see other Van Rooij repository), 
# but since we now also incorporate published data of healthy hearts, 
# this scripts repeats earlier analyses using the Seurat package.
#
# Original RaceID2 analysis by Joep Eding
# Seurat analysis including healhty hearst by m.wehrens@hubrecht.eu
# 2021-06
# 
# Outline of script:
# First I define the general analysis that is performed as a function;
# this since I want to perform several versions of the analysis. 
# (The Seurat default versions emphasize different features than
# the RaceID2 analysis due to normalization choices.)
#
# I will then
# 1. Create UMAP to compare the different data sets
# (..)












