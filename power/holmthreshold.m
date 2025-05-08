% This function calculates the threshold for Holm (stepdown Bonferroni)
% multiple testing correction based on the following inputs
%
%   base_alpha      alpha value you want to start from, typically 0.05
%   num_dvs         number of dependent variables/univariate tests

function threshold = holmthreshold(base_alpha,num_dvs)

end_value = 0;

    for n = 1:num_dvs
     end_value = end_value + base_alpha/n;
    end

 threshold = end_value/num_dvs;

end
    