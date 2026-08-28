function [intact_fractions] = LCN_calc_intact_tracer_hill(params,times)
%
% the function y(t) = b/(1+(a/t)^n) will be fitted with b set to 1
% b = asymptotic value, assumed 1 at time 0 (100% intact tracer)
% a = t where the function reaches half of its maximum (b/2)
% n determines the steepness of the slope
%
% FORMAT: intact_fractions = LCN_calc_intact_tracer_hill(params,times)
%
% params: [a, n]
% times: array of times in minutes
% intact_fractions: intact fraction of tracer (between 0 and 1)
%__________________________________________________________________________
%
% author: Natalie Nelissen, Patrick Dupont
% date:   February, 2006
% history: 	
%__________________________________________________________________________
% @(#)LCN_calc_intact_tracer_hill.m	     1.0     last modified: 20006/02/16

intact_fractions = 1./(1+(params(1)./times).^params(2));