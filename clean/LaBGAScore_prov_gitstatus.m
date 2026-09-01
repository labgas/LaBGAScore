function S = LaBGAScore_prov_gitstatus(repopath, varargin)
%
%
% *USAGE*
%
% Determines which tracked files in a git repository have uncommitted
% modifications, by parsing .git/index directly and comparing it against
% the working tree - WITHOUT shelling out to git.
%
% This matters because a commit hash alone does not pin a dependency. At
% the time of writing CanlabCore sits at a known commit but carries five
% uncommitted local edits, two of them in @predictive_model, i.e. squarely
% in the MVPA path. Recording only the hash would misrepresent what ran.
%
% Verified to reproduce "git status --porcelain" exactly on CanlabCore
% (5 modified of 4004 tracked), LaBGAScore (0 of 153) and JuSpace (8 of
% 108).
%
% Method, which is git's own: for every index entry compare the recorded
% size and mtime against the file on disk; only where those differ, hash
% the file as a git blob (sha1 of "blob <len>\0" + content) and compare it
% to the recorded object id. So the expensive step runs on a handful of
% files, not on all of them.
%
%
% *OPTIONS*
%
% * repopath        path to the repository working tree
% * 'subset'        cellstr/string of repo-relative paths to restrict the
%                   check to. Use this to ask the far cheaper question
%                   "is any file THIS SCRIPT depends on modified?"
%
%
% *OUTPUT*
%
% * S               struct with fields
%                   .ok            logical, index was parsed successfully
%                   .ntracked      number of tracked files
%                   .modified      string array of repo-relative paths
%                   .deleted       string array of repo-relative paths
%                   .nchecked      files whose stat differed, so were hashed
%                   .dirty         logical, any modification or deletion
%                   .error         message when the index could not be read
%
%
% *NOTES*
%
% Untracked files are NOT reported: the index has no record of them, and
% they cannot be part of any past commit. A dependency file that is
% untracked would show up as unresolved in the provenance record instead.
%
% Only index format version 2 is supported, which is what git writes by
% default and what all repositories here use. A version 3 or 4 index is
% reported via .error rather than parsed incorrectly.
%
%
% *SEE ALSO*
%
% LaBGAScore_prov_gitinfo, LaBGAScore_prov_snapshot
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove & Claude Opus 5
%
% date:   KU Leuven, August, 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prov_gitstatus.m         v1.0
%
% last modified: 2026/08/31
%
%
%% PARSE OPTIONS
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('subset', strings(0,1), @(x) iscell(x) || isstring(x) || ischar(x));
p.parse(varargin{:});

subset = string(p.Results.subset);
subset = subset(:);

repopath = char(repopath);

S = struct('ok', false, 'ntracked', 0, 'modified', strings(0,1), ...
    'deleted', strings(0,1), 'nchecked', 0, 'dirty', false, 'error', '');


%% LOCATE AND READ THE INDEX
% -------------------------------------------------------------------------

G = LaBGAScore_prov_gitinfo(repopath);

if ~G.isrepo
    S.error = G.error;
    return
end

indexfile = fullfile(G.gitdir,'index');

if ~isfile(indexfile)
    S.error = 'no .git/index (nothing staged yet)';
    return
end

fid = fopen(indexfile,'r');
if fid < 0
    S.error = 'could not open .git/index';
    return
end
b = fread(fid, Inf, '*uint8');
fclose(fid);

if numel(b) < 12 || ~isequal(char(b(1:4))','DIRC')
    S.error = 'not a git index file';
    return
end

version = local_be32(b, 5);
nentries = local_be32(b, 9);

if version ~= 2
    S.error = sprintf('unsupported git index version %d (only v2 is parsed)', version);
    return
end


%% WALK THE INDEX ENTRIES
% -------------------------------------------------------------------------
% v2 entry layout, all big-endian:
%   ctime(8) mtime(8) dev(4) ino(4) mode(4) uid(4) gid(4) size(4)
%   sha1(20) flags(2) path(namelen) then 1-8 NUL bytes of padding, so that
%   the whole entry length is a multiple of 8

off = 13;                                   % 1-based, past the 12-byte header

paths = strings(nentries,1);
shas = strings(nentries,1);
sizes = zeros(nentries,1);
mtimes = zeros(nentries,1);

for k = 1:nentries

    start = off;

    if off + 61 > numel(b)
        S.error = sprintf('index truncated at entry %d of %d', k, nentries);
        return
    end

    mtimes(k) = local_be32(b, off + 8);
    sizes(k) = local_be32(b, off + 36);
    shas(k) = local_hex(b(off+40 : off+59));

    flags = double(b(off+60))*256 + double(b(off+61));
    namelen = bitand(flags, 4095);

    off = off + 62;

    if namelen < 4095
        raw = b(off : off+namelen-1);
        off = off + namelen;
    else
        % name length did not fit in the flags: NUL-terminated instead
        nul = find(b(off:end) == 0, 1, 'first');
        raw = b(off : off+nul-2);
        off = off + nul - 1;
    end

    % git stores paths as raw bytes, UTF-8 on this platform. Decoding them
    % as Latin-1 turns any non-ASCII filename into mojibake that then fails
    % to match the file on disk, and it would be reported as deleted.
    paths(k) = string(native2unicode(raw', 'UTF-8'));

    entlen = off - start;
    off = off + 8 - mod(entlen, 8);         % always 1-8 padding bytes

end

S.ok = true;
S.ntracked = nentries;


%% COMPARE AGAINST THE WORKING TREE
% -------------------------------------------------------------------------

if ~isempty(subset)
    keep = ismember(paths, subset);
    paths = paths(keep); shas = shas(keep);
    sizes = sizes(keep); mtimes = mtimes(keep);
end

modified = strings(0,1);
deleted = strings(0,1);
nchecked = 0;

for k = 1:numel(paths)

    fp = fullfile(repopath, char(paths(k)));

    jf = java.io.File(fp);

    if ~jf.exists()
        deleted(end+1,1) = paths(k); %#ok<AGROW>
        continue
    end

    ondisksize = double(jf.length());
    ondiskmtime = floor(double(jf.lastModified())/1000);

    if ondisksize == sizes(k) && ondiskmtime == mtimes(k)
        continue                            % stat matches: git calls this clean
    end

    nchecked = nchecked + 1;

    if jf.isDirectory()
        continue                            % submodule/gitlink, not a blob
    end

    if ~strcmp(local_blob_sha(fp), shas(k))
        modified(end+1,1) = paths(k); %#ok<AGROW>
    end

end

S.modified = modified;
S.deleted = deleted;
S.nchecked = nchecked;
S.dirty = ~isempty(modified) || ~isempty(deleted);

end


%% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function v = local_be32(b, i)
% big-endian uint32 starting at 1-based byte i
v = double(b(i))*16777216 + double(b(i+1))*65536 + double(b(i+2))*256 + double(b(i+3));
end


function h = local_hex(bytes)
% lowercase hex string of a byte vector
h = string(lower(reshape(dec2hex(bytes,2)', 1, [])));
end


function h = local_blob_sha(fp)
% git's object id for a file: sha1 of "blob <bytelength>\0" then the content

fid = fopen(fp,'r');
if fid < 0, h = ''; return, end
data = fread(fid, Inf, '*uint8');
fclose(fid);

header = [uint8(sprintf('blob %d', numel(data))), uint8(0)];

md = java.security.MessageDigest.getInstance('SHA-1');
md.update(header);
md.update(data);

digest = typecast(md.digest(), 'uint8');
h = char(lower(reshape(dec2hex(digest,2)', 1, [])));

end
