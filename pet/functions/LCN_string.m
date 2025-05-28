function new_string = LCN_string(stringarray)
% this function will replace _ by \_ which is useful for printing strings
% because _ is leading to a subscript for the next character
%
% This routine is supplied as is. 
%
% Comments or questions can be send to:
% Patrick.Dupont@med.kuleuven.be
%__________________________________________________________________________
%
% author: Patrick Dupont
% date:   December 07, 2016
% history: bug fix
%__________________________________________________________________________
% @(#)LCN_string.m             v0.2               last modified: 2019/06/05

a = strfind(stringarray,'_');
new_string = blanks(length(stringarray)+length(a));
if isempty(a)
   new_string = stringarray;
else
   new_string(1:a(1)-1) = stringarray(1:a(1)-1);
   for i = 1:length(a)-1
       new_string(a(i)+i-1:a(i)+i) = '\_';
       new_string(a(i)+i+1:a(i+1)+(i-1)) = stringarray(a(i)+1:a(i+1)-1);
   end
   new_string(a(end)+(length(a)-1):a(end)+length(a)) = '\_';
   new_string(a(end)+length(a)+1:end) = stringarray(a(end)+1:end);
end
end