function [intact_fractions] = LCN_calc_intact_tracer_biexp_delay(params,times)
%
% calculation of the intact fraction of a tracer using the function
% a*exp(-ln(2)*(t-t_delay)/alfa)+(1-a)*exp(-ln(2)*(t-t_delay)/beta)
%
% FORMAT: intact_fractions = LCN_calc_intact_tracer_biexp_delay(params,times)
%
% params: [alfa, a, beta, t_delay]
% times: array of times in minutes
% intact_fractions: intact fraction of tracer (between 0 and 1)
%__________________________________________________________________________
%
% author: Tom Muylle, Natalie Nelissen, Patrick Dupont
% date:   December, 2005	
% history: 	
%__________________________________________________________________________
% @(#)LCN_calc_intact_tracer_biexp_delay.m	1.0  last modified: 20005/12/30

intact_fractions =  params(2)*exp(-log(2)*(times-params(4))/params(1)) ...
                + (1-params(2))*exp(-log(2)*(times-params(4))/params(3));