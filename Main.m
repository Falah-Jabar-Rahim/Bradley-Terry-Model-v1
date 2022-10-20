% Bradley-Terry-Model
% (c) Falah Jabar (falah.jabar@lx.it.pt)

close all; clear; clc;

%% Example 1 taken from R. A. Bradley, Paired Comparisons: Some Basic Procedures and
 % Examples, Handbook of Statistics, vol. 4, P. R. Krishnaiah and P. K.
% P=[ 0 28 15 23 ; 112 0 46 47 ; 39 17 0 0 ; 34 11 0 0 ];
 %% Example 2 from J. C. Handley, “Comparative Analysis of Bradley-Terry and
                                     % Thurstone-Mosteller Paired Comparison Models for Image Quality
                                     % Assessment,” Proc. of PICS?: Image: Processing, Quality, Capture,
                                     % Systems, Montreal, Canada, Apr. 2001.
P=[0 26 28 22;64 0 46 34; 62 44 0 26; 68 56 64 0];

[Q, Lb, Ub] = BT_model(P);

