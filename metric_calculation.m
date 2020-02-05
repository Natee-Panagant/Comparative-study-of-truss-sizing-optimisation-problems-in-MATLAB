clear all; close all; clc;

load RefPoint.mat
load RefFront.mat

itest=1;    % test index
ialgo=1;    % algorithm index
irun=1;     % run index
load(['rst_00' num2str(itest) '_00' num2str(ialgo)]);
hv=hypervolume(rst(irun).fpareto{end},RefPoint{itest})
gd=generational_distance(rst(irun).fpareto{end},RefFPareto{itest})
igd=inverted_generational_distance(rst(irun).fpareto{end},RefFPareto{itest})
ste=spacing_to_extent(rst(irun).fpareto{end})