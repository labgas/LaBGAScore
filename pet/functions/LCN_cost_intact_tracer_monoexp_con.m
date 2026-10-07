function [res] = LCN_cost_intact_tracer_monoexp_con(params)
%
% costfunction for the calculation of the intact fraction of a tracer
% using a mono exponential function with amplitude 1 (assuming that at time 
% t=0, the intact fraction is 1)
%
% FORMAT: res = LCN_cost_intact_tracer_monoexp_con(params)
%
% params: [alfa] (see LCN_calc_intact_tracer_monoexp_con.m)
% res: residual error
%
% global variables:
%       TIME_METAB: time (in min) when metabolite samples are taken
%       FRACTION_INTACT_TRACER: fraction intact tracer 
%                               (values between 0 and 1, unitless)
%       WEIGHTS_METAB: weights for each sample to be used in costfunction
%                      (values between 0 and 1).
%__________________________________________________________________________
%
% author: Patrick Dupont
% date:   February, 2006	
% history: 	
%__________________________________________________________________________
% @(#)LCN_cost_intact_tracer_monoexp_con.m	1.0  last modified: 20006/02/16

global TIME_METAB
global FRACTION_INTACT_TRACER
global WEIGHTS_METAB

tmp = LCN_calc_intact_tracer_monoexp_con(params,TIME_METAB);
res = norm(WEIGHTS_METAB.*(tmp-FRACTION_INTACT_TRACER));