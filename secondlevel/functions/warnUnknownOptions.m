function warnUnknownOptions(opts, known, pipelineName)
% warnUnknownOptions  Warn about opts fields the pipeline will silently ignore.
%
%   warnUnknownOptions(opts, known, pipelineName)
%
%   MATLAB struct field names are case-sensitive, so a misspelling such as
%   opts.nrepeats when the pipeline reads opts.nRepeats is not an error -- the
%   field is simply never looked at and the default is used instead. Both
%   calling scripts in this repo shipped exactly that bug for every run, which
%   is what motivated this check.
%
%   Inputs
%     opts          the resolved options struct
%     known         cellstr of every field this pipeline reads
%     pipelineName  char, used in the warning identifier and message
%
%   Call this AFTER the defaults block, so that fields the pipeline itself
%   inserted are already present in both opts and known.
%
%   Near-miss detection: an unknown field whose lowercase form matches a known
%   field's lowercase form is reported as a probable capitalization slip, since
%   that is by far the most common way this goes wrong.
%
%   See also validateCovariates.

extra = setdiff(fieldnames(opts), known);

if isempty(extra)
    return
end

knownLower = lower(known);

for i = 1:numel(extra)
    f = extra{i};
    hit = find(strcmp(lower(f), knownLower), 1);
    if ~isempty(hit)
        warning([pipelineName ':optionCase'], ...
            'opts.%s is ignored -- did you mean opts.%s? Struct field names are case-sensitive, so the default was used instead.', ...
            f, known{hit});
    else
        warning([pipelineName ':optionUnknown'], ...
            'opts.%s is not an option of %s and will be ignored.', f, pipelineName);
    end
end

end
