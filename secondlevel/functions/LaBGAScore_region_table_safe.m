function [rpos, rneg, results_table] = LaBGAScore_region_table_safe(fcn, r, varargin)
% LaBGAScore_region_table_safe  Call a region-table function and always get 3 outputs.
%
%   [rpos, rneg, results_table] = LaBGAScore_region_table_safe(fcn, r, ...)
%
%   Thin wrapper around CanlabCore's @region/table (or the vendored
%   LaBGAScore_region_table), passing everything through unchanged EXCEPT that an
%   empty result yields an empty results_table instead of an error.
%
%   WHY THIS EXISTS
%   @region/table assigns its third output, results_table, only on the path that
%   actually renders a table. Two of its paths return without assigning it:
%
%     - the legacy printer, used when autolabel_regions throws
%     - the "No regions to display" branch, taken when region_table is empty
%
%   A three-output call on either path dies with
%
%       Output argument "results_table" (and maybe others) not assigned
%       during call to "region/table".
%
%   so a contrast with nothing suprathreshold - a perfectly ordinary result -
%   aborts the whole display script. That is what happened to the proj_cfs
%   parcelwise c2a run on 2026-09-09: 571 parcels, nothing displayable at one
%   threshold, script dead.
%
%   INPUTS
%     fcn       function handle, @table or @LaBGAScore_region_table
%     r         region object array
%     varargin  passed straight through to fcn
%
%   OUTPUTS
%     rpos, rneg     as returned by fcn; on the empty path, rpos = r, rneg = []
%     results_table  the table fcn produced, or an EMPTY table if it produced none
%
%   Only the specific "results_table not assigned" failure is absorbed. Any other
%   error is re-thrown, so real problems still surface.
%
%   See also table, LaBGAScore_region_table.

% AUTHOR
% Lukas Van Oudenhove, KU Leuven
%
% ..

results_table = table();

try

    [rpos, rneg, results_table] = fcn(r, varargin{:});

catch ME

    if contains(ME.message, 'results_table')

        % nothing displayable: still return the split regions if we can
        try
            [rpos, rneg] = fcn(r, varargin{:});
        catch
            rpos = r;
            rneg = [];
        end

        results_table = table();
        fprintf('\nno region table produced for this contrast (nothing displayable); continuing\n\n');

    else

        rethrow(ME);

    end

end

end % function
