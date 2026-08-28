function [res] = LCN_cost_intact_tracer_hill(params)
%
% costfunction for the calculation of the intact fraction of a tracer
% using a Hill function (see LCN_calc_intact_tracer_hill.m)
% We don't assume that at time t=0, the intact fraction is 1, but allow 
% a delay before metabolites appear
%
% FORMAT: res = LCN_cost_intact_tracer_hill(params)
%
% params = [a, n] (see LCN_calc_intact_tracer_hill.m)
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
% author: Natalie Nelissen, Patrick Dupont
% date:   February, 2006
% history: 	
%__________________________________________________________________________
% @(#)LCN_cost_intact_tracer_hill.m	    1.0      last modified: 20006/02/16

global TIME_METAB
global FRACTION_INTACT_TRACER
global WEIGHTS_METAB

tmp = LCN_calc_intact_tracer_hill(params,TIME_METAB);
res = norm(WEIGHTS_METAB.*(tmp-FRACTION_INTACT_TRACER));