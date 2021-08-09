%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings.m
% Default settings for PEC80 controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear Workspace and Screen
clear;
clc;

%Cycle Time for PEC80 Tasks
Base_cycle_time  = 25e-6;   %[s]
Task_A_mult      = 4;       %[1]
Task_B_mult      = 20;      %[1]
Task_C_mult      = 1;       %[1]
Task_A           = Task_A_mult * Base_cycle_time;   %[s]
Task_B           = Task_B_mult * Task_A;            %[s]
Task_C           = Task_C_mult * Task_B;            %[s]
Ts               = 0.01;   %Background task cycle time in [s]
T                = Task_A; %T Sample in [s]
