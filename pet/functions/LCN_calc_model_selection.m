function [AIC,SC] = LCN_calc_model_selection(N,SS,P)
% format [AIC,SC] = LCN_calc_model_selection(N,SS,P);
% 
% This script will calculate criteria for model selection
% INPUT
%   N = the number of samples
%   P = the number of parameters in the model
%   SS = the sum of squared differences between measurement and model prediction
% OUTPUT
%   AIC = Akaiki criterion
%   SC  = Schwarz criterion
%
% author: Patrick Dupont
% date: July, 2025
% history:
%__________________________________________________________________________
% @(#)LCN_calc_model_selection.m        v0.1      last modified: 2025/07/02

AIC = N.*log(SS/N)+2.*P;
SC  = N.*log(SS/N)+P.*log(N);

end
