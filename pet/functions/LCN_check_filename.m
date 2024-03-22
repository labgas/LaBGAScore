function [filename,go] = LCN_check_filename(dir_data,template_string)
% This script will check if we find uniquely a filename in a directory
% dir_data using the template_string as search string.
% If we find a unique file, this filename will be given as output. The
% variable go is also an output which is 1 if we find a unique filename and
% 0 otherwise.
%
% author: Patrick Dupont
% date:   July, 2023
% history: 
% _________________________________________________________________________
% @(#)LCN_check_filename.m          v0.1          last modified: 2023/07/11       
       
cd(dir_data)
dirlist = dir(template_string);
if size(dirlist,1) == 1
   filename = fullfile(dir_data,dirlist(1).name);
   go = 1;
else
   filename = '';
   go = 0;
end    
end