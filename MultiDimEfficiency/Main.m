% This script documents and calls each of the sub-scripts used to reproduce
% the figures and analyses from
%
% Efficient Sensory Coding of Multidimensional Stimuli
% Thomas E. Yerxa, Eric Kee, Michael R. DeWeese, Emily A. Cooper
%

% For the code to run, the optimizer must be in your matlab path
addpath(genpath('../Optimizer'));

% Overall summary of 1-D and 2-D populations
Summary1D;
Summary2D;

% Illustration of population parameters in 10-D
PopulationParameters;

% Examples of efficient 2-D populations
ExamplePopsPrep;     % first generate the manifold
ExamplePopsPlotIt;   % then generate plots

% Analysis of correlation between tuning curve characteristics and
% probability for example populations
CorrelationAnalysis;

% Large scale simulations and analysis of correlations
SimsPrep1;      % first generate the manifolds
SimsPrep2;      % the compute tuning curve characteristics
SimsPlotIt;       % then generate plots

% Analysis Fisher Information in example populations
FisherInfo;
