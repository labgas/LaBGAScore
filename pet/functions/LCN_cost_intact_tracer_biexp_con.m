function [res] = LCN_cost_intact_tracer_biexp_con(params)
%
% costfunction for the calculation of the intact fraction of a tracer
% using a sum of two exponentials 
% We assume that at time t=0, the intact fraction is 1, but allow 
% a fixed delay before metabolites appear
%
% FORMAT: res = LCN_cost_intact_tracer_biexp_con(params,delay)
%
% params: [alfa, a, beta] (see LCN_calc_intact_tracer_biexp_con.m)
% delay: time (in min) since when metabolisation starts
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
% author: Tom Muylle, Natalie Nelissen, Patrick Dupont
% date:   December, 2005	
% history: 	
%__________________________________________________________________________
% @(#)LCN_cost_intact_tracer_biexp_con.m	1.0  last modified: 20005/12/30

global TIME_METAB
global FRACTION_INTACT_TRACER
global WEIGHTS_METAB

tmp = LCN_calc_intact_tracer_biexp_con(params,TIME_METAB);
res = norm(WEIGHTS_METAB.*(tmp-FRACTION_INTACT_TRACER));