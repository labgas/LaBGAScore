function [C_calc] = LCN_calc_model_2T4k_vasc1k(params,midscantimes,frameduration,CA_FINETIMES,CA_WB_FINETIMES,FINETIMES,calc_option,STEP)
%
% calculates the output value for the LCN_calc_model_2TCM_K1.m (with parameters params)
% FORMAT C_calc = LCN_calc_model_2TCM_K1(params, midscantimes, frameduration, calc_option)
%
% params: values for [K1,k2,k3,k4,Vb,K1v]
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
% 
% C_calc: value of the model calculated at time times
%__________________________________________________________________________
%
% author: 	Patrick Dupont
% date: 	June, 2025
% history: 	
%__________________________________________________________________________
% @(#)LCN_calc_model_2T4k_vasc1k.m	   v0.1       last modified: 2025/06/24

K1  = params(1);
k2  = params(2);
k3  = params(3);
k4  = params(4);
Vb  = params(5);
K1v = params(6);

ksum   = k2+k3+k4;
delta  = sqrt(ksum^2-(4*k2*k4));
theta1 = (ksum+delta)/2;
theta2 = (ksum-delta)/2;
phi1   = K1*(theta1-k3-k4)/delta;
phi2   = K1*(k3+k4-theta2)/delta;

% calculate the integrals at finetimes
int_finetimes1 = STEP*cumsum(CA_FINETIMES.*exp(theta1*FINETIMES));
int_finetimes2 = STEP*cumsum(CA_FINETIMES.*exp(theta2*FINETIMES));
int_finetimes3 = STEP*cumsum(CA_WB_FINETIMES);

C_calc = (1-Vb)*(phi1*exp(-theta1*FINETIMES).*int_finetimes1+phi2*exp(-theta2*FINETIMES).*int_finetimes2)+Vb*(CA_WB_FINETIMES+K1v.*int_finetimes3);

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