function values = SearchGammaTable(targets, input, table)
% values = SearchGammaTable(targets, input, table)
%
% Return the [0-1] entry from the passed table that produces
% output closest to the [0-1] target.
%
% The targets are assumed to be a row vector.
% The table is assumed to be a column vector.
% The returned indices and values are are row vectors.
%
% Works by using Matlab's interp1, with the output as the x values and
% the input as the f(x) values.
% 
% I suspect that this is a fast Matlab implementation, but those who want
% to try are welcome to try to do better.  (Remember, though, that this
% routine gains in efficiency the more searches are done at once.
% This is because it contains no dreaded loops.)
%
% 4/2/94		dhb		Added code that checks for special case of zero output.
% 4/4/94		dhb		Fixed code added on 4/2.
% 4/5/94		jms		Fixed code added on 4/2.
% 1/21/95		dhb		Write search as a loop.  Loses time and elegance,
%						but prevents allocation of arrays that may be huge.
% 11/16/06      dhb     Renamed as SearchGammaTable.
%               dhb     Start work on converting to [0-1] universe.  Change
%                       name and interface.
% 11/20/06      dhb     Finish update by calling through MATLAB's interpolation function.
% 9/15/08       dhb     Handle case where there are a bunch of zeros at the beginning of gamma table.
% 5/26/12       dhb     Improve comment.  This was not doing exhaustive search.

% Check dimensions
[m,n] = size(targets);
if (m ~= 1)
    error('Passed targets should be a row vector');
end
[mi,ni] = size(input);
if (ni ~= 1)
    error('Passed input should be a column vector');
end
[mt,nt] = size(table);
if (nt ~= 1)
    error('Passed table should be a column vector');
end
if (mi ~= mt || ni ~= nt)
    error('Input and table must be the same size');
end

% Handle problem that for some monitors, the output is 0 for
% input values up to some threshold.  This causes interp1
% to crash.  We handle this by getting rid of the intermediate
% zeros from the input, if they are there.  This choice means 
% that when we ask for 0 out, we get 0 as the answer.
index = find(table == 0);
index1 = find(table ~= 0);
if ~isempty(index)
    table = table([index(1) ; index1]);
    input = input([index(1) ; index1]);
end

% Invert via linearly interpolation of the passed table
values = interp1(table, input, targets', 'linear')';


