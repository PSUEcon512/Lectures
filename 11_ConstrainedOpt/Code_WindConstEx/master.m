


clc;
clear;
close all

%MASTER file which runs code in order of dependencies..


tic;
%Setup data structures 
setup;

%Estimate the model and find standard errors
Main_mpec;
toc;
