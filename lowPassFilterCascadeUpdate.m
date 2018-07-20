function lowPassFilterCascade = lowPassFilterCascadeUpdate(lowPassFilterCascade, input)
%LOWPASSFILTERCASCADEUPDATE Updates low-pass filter cascade outputs.
%   lowPassFilterCascadeUpdate(lowPassFilterCascade, input)) returns the
%   updated low-pass filter cascade structure for a specified
%   lowPassFilterCascade and input.
%
%   See:
%   https://en.wikipedia.org/wiki/Low-pass_filter
%   https://www.electronics-tutorials.ws/filter/filter_2.html

    for filterIndex = 1:lowPassFilterCascade.numberOfFilters
        if filterIndex > 1
            input = lowPassFilterCascade.outputs(filterIndex - 1);
        end
        lowPassFilterCascade.outputs(filterIndex) = lowPassFilterCascade.outputs(filterIndex) + ((input - lowPassFilterCascade.outputs(filterIndex)) * lowPassFilterCascade.coefficient);
    end

end
