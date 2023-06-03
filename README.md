# ELISA'
This repository contains a script for data analysis and visualization of ELISA results. The script utilizes various libraries, including pandas, numpy, matplotlib, scipy, pick, and tabulate. The main functionalities of the script are:

    Importing data from CSV files: The script provides functions import1d and importdata to import 1D arrays and tabular data from CSV files, respectively.
    Calculating statistical measures: The script includes functions calcc and calcC to calculate the coefficient of variation (CV) and the median of a dataset, respectively.
    Fitting a logistic 4-parameter curve: The script implements the logistic 4-parameter curve equation and provides functions for fitting the curve to data and evaluating the fitted curve at given x-values.
    Concentration calculation: The script includes a function concentration to calculate the concentration from an OD value using the fitted logistic curve parameters.
    Data preprocessing: The script provides a function cuttable to extract a specific table from a CSV file based on a given keyword.
    Data analysis: The script includes functions for reading samples, calculating averages, performing quality control checks, and creating a results table.
    ELISA unit conversion: The script provides a function ODbyELISA to convert OD values to ELISA units using the 1Standard as a reference.
    Visualization: The script includes functions for creating plots and saving the results.

To use the script, you need to have the following dependencies installed:

    pandas
    numpy
    matplotlib
    scipy
    pick
    tabulate
