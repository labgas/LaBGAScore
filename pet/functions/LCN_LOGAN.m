function [VD,error] = LCN_LOGAN(MIDSCANTIMES,C_MEASURED,CA_TIMES,CA,logan_start_time,figure_nr)
% LCN_LOGAN.m
%
% calculates the slope of the Logan plot. The slope is the distribution 
% volume VD
%
% FORMAT [VD,error] = LCN_LOGAN(MIDSCANTIMES,C_MEASURED,CA_TIMES,CA,logan_start_time);
%
% INPUT VARIABLES
%-----------------
% MIDSCANTIMES = list of midscan times of each frame (in min p.i.).
% C_MEASURED   = the TAC at MIDSCANTIMES
% CA_TIMES     = list of times when a blood sample was taken (in min p.i.)
% CA           = the metabolite corrected input function in the same units
% as C_MEASURED at CA_TIMES
% logan_start_time = start time in minutes for the Logan analysis
% figure_nr = 0 for no figure or an integer indicating the figure number
%
% OUTPUT VARIABLES
%-----------------
% VD = distribution volume
% error = residual error of the Logan analysis
%__________________________________________________________________________
%
% author: 	Patrick Dupont
% date: 	January, 2015
% history: 	September 2023: we don't plot the first data point since
%                           this is not relevant and the values can be
%                           artificially high.
%__________________________________________________________________________
% @(#)LCN_LOGAN.m                         0.2     last modified: 2023/09/26

STEP = 0.01; % in min to interpolate at finer times to calculate the integrals

% make the inputvectors all the same
if size(MIDSCANTIMES,1) == 1
   MIDSCANTIMES = MIDSCANTIMES';
end
if size(C_MEASURED,1) == 1
   C_MEASURED = C_MEASURED';
end
if size(CA_TIMES,1) == 1
   CA_TIMES = CA_TIMES';
end
if size(CA,1) == 1
   CA = CA';
end

% calculate which frames may be taken into account (time after
% logan_start_time but before end of bloodsampling).
index_logan = find(MIDSCANTIMES > logan_start_time);
index_time  = find(MIDSCANTIMES < max(CA_TIMES));
indexvalues = intersect(index_logan,index_time);

FINETIMES      = [0:STEP:min(max(MIDSCANTIMES),max(CA_TIMES))]';

C_MEASURED_FINETIMES        = interp1q([0; MIDSCANTIMES],[0; C_MEASURED],FINETIMES);
INT_C_MEASURED_FINETIMES    = [0; (FINETIMES(2:end)-FINETIMES(1:end-1))].*cumsum(C_MEASURED_FINETIMES);
INT_C_MEASURED_MIDSCANTIMES = interp1q(FINETIMES,INT_C_MEASURED_FINETIMES,MIDSCANTIMES);

CA_FINETIMES = interp1q([0; CA_TIMES],[0; CA],FINETIMES);
INT_CA_FINETIMES    = [0; (FINETIMES(2:end)-FINETIMES(1:end-1))].*cumsum(CA_FINETIMES);
INT_CA_MIDSCANTIMES = interp1q(FINETIMES,INT_CA_FINETIMES,MIDSCANTIMES);

y_all = INT_C_MEASURED_MIDSCANTIMES./C_MEASURED;
y = y_all(indexvalues);
x_all = INT_CA_MIDSCANTIMES./C_MEASURED;
x = x_all(indexvalues);
       
[Vdfit,error1] = polyfit(x,y,1);
       
if figure_nr > 0
   figure(figure_nr)
   plot(x_all(2:end),y_all(2:end),'or')
   hold on
   plot(x,y,'ob')
   plot([x(1) x(end)],[Vdfit(1)*x(1)+Vdfit(2) Vdfit(1)*x(end)+Vdfit(2)]) 
   title('Logan plot')
   xlabel('(input integral) / tissue')
   ylabel('(tissue integral) / tissue')
   axvalues = axis;
   text(axvalues(1)*1.05,axvalues(3)+0.9*(axvalues(4)-axvalues(3)),['VD = ' num2str(Vdfit(1))]);
end
VD    = Vdfit(1); % slope
error = error1.normr; 
end