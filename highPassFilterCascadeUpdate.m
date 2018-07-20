function highPassFilterCascade = highPassFilterCascadeUpdate(highPassFilterCascade, input)
%HIGHPASSFILTERCASCADEUPDATE Updates high-pass filter cascade outputs.
%   highPassFilterCascadeUpdate(highPassFilterCascade, input)) returns the
%   updated high-pass filter cascade structure for a specified
%   highPassFilterCascade and input.
%
%   See:
%   https://en.wikipedia.org/wiki/High-pass_filter
%   https://www.electronics-tutorials.ws/filter/filter_2.html

    for filterIndex = 1:highPassFilterCascade.numberOfFilters
        if filterIndex > 1
            input = highPassFilterCascade.outputs(filterIndex - 1);
        end
        highPassFilterCascade.outputs(filterIndex) = highPassFilterCascade.coefficient * (highPassFilterCascade.outputs(filterIndex) + input - highPassFilterCascade.inputs(filterIndex));
        highPassFilterCascade.inputs(filterIndex) = input;
    end

end
