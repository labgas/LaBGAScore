function [intact_fractions] = LCN_calc_intact_tracer_monoexp(params,times)
%
% calculation of the intact fraction of a tracer using the function
% a*exp(-ln(2)*t/alfa)
%
% FORMAT: intact_fractions = LCN_calc_intact_tracer_monoexp(params,times)
%
% params: [alfa, a]
% times: array of times in minutes
% intact_fractions: intact fraction of tracer (between 0 and 1)
%__________________________________________________________________________
%
% author: Patrick Dupont
% date:   February, 2006	
% history: 	
%__________________________________________________________________________
% @(#)LCN_calc_intact_tracer_monoexp.m	1.0      last modified: 20006/02/16

intact_fractions =  params(2)*exp(-log(2)*times/params(1));