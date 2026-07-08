function [q,aprioriprob] = LaBGAScore_Storey_FDR(p)

% This function implements Storey's positive FDR as in SAS proc multtest
%
% 1. It checks whether the condition of at least one true null hypothesis 
%       to apply positive FDR has been
%       fulfilled, otherwise defaults back to Benjamini-Hochberg FDR
% 2. After applying Storey's fdr, it applies the constraint q >= p to
%       prevent corrected p-values from being lower than raw p-values
%
% INPUT
%
% p : a nx1 vector of raw p-values
%
% OUTPUTS
%
% q: a nx1 vector of FDR-corrected p-values
%
% aprioriprob: estimate of true null hypotheses

% Apply Storey's FDR
[~, q, aprioriprob] = mafdr(p);

if aprioriprob > 0.99
    % Fallback to Benjamini-Hochberg FDR
    q = mafdr(p, 'BHFDR', true);
    
    fprintf('\naprioriprob = %s, falling back to Benjamini-Hochberg FDR\n',num2str(aprioriprob));
   
else
    % Keep Storey FDR & enforce constraint q >= p (as in SAS proc multtest) 
    for j = 1:length(q)
        if q(j) < p(j)
            q(j) = p(j);
        end
    end
    
    fprintf('\naprioriprob = %s, applying Storey FDR\n',num2str(aprioriprob));
    
end

