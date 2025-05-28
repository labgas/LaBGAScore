function [C_calc] = LCN_calc2_model_2T4k_Vb(params,midscantimes,frameduration,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,calc_option,STEP)
%
% calculates the output value for the model_2T4k_Vb (with parameters params)
% FORMAT C_calc = LCN_calc_model_2T4k_Vb(params, midscantimes, frameduration, calc_option)
%
% params: values for [K1,k2,k3,k4,Vb]
% midscantimes  = list of midscan times of each frame (in min p.i.).
% frameduration = list of frame durations (in min).
% CA_FINETIMES: plasma concentration (kBq/ml) at FINETIMES
% CA_WB_FINETIMES: whole blood concentration (kBq/ml) at FINETIMES
% FINETIMES: vector of times (in minutes) with small step size STEP
% calc_option: parameter determing the way we calculate the output of a
% model.
%   1 - the model is calculated as the integral of the output concentration
%       devided by the frameduration.
%   2 - the model is calculated at the midscantime.
% STEP: step size of finer time sampling (in minutes)
% 
% C_calc: value of the model calculated at time times
%__________________________________________________________________________
%
% author: 	Patrick Dupont
% date: 	November, 2024
% history: 	
%__________________________________________________________________________
% @(#)LCN_calc2_model_2T4k_Vb.m	      v0.1        last modified: 2024/11/22

K1 = params(1);
k2 = params(2);
k3 = params(3);
k4 = params(4);
Vb = params(5);

ksum   = k2 + k3 + k4;
tmp = sqrt(ksum ^ 2 - (4 * k2 * k4));
alfa1 = (ksum - tmp) / 2;
alfa2 = (ksum + tmp) / 2;

factor = K1 / (alfa2 - alfa1);

% calculate the integrals at finetimes
int_finetimes1 = STEP*cumsum(CA_FINETIMES .* exp(alfa1 * FINETIMES));
int_finetimes2 = STEP*cumsum(CA_FINETIMES .* exp(alfa2 * FINETIMES));

C_calc = factor .* ((k3 + k4 - alfa1) * exp(-alfa1 * FINETIMES) .* int_finetimes1 - (k3 + k4 - alfa2) * exp(-alfa2 * FINETIMES) .* int_finetimes2) + Vb*CA_WB_FINETIMES;

if calc_option == 1
   cum_cmod = STEP.*cumsum(C_calc);
   ENDFRAME   = midscantimes + frameduration/2;
   STARTFRAME = midscantimes - frameduration/2;
   cum_cmod_end   = interp1q(FINETIMES,cum_cmod, ENDFRAME);
   cum_cmod_start = interp1q(FINETIMES,cum_cmod, STARTFRAME);
   C_calc = (cum_cmod_end - cum_cmod_start)./(ENDFRAME-STARTFRAME);
elseif calc_option == 2
   C_calc = interp1q(FINETIMES, C_calc, midscantimes);
end

if size(C_calc,2) > 1
   C_calc = C_calc';
end

end