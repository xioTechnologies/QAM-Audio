function lowPassFilterCascade = lowPassFilterCascadeInitialise(numberOfFilters, cornerFrequency, sampleFrequency)
%LOWPASSFILTERCASCADEINITIALISE Initialises low-pass filter cascade.
%   lowPassFilterCascadeInitialise(numberOfFilters, cornerFrequency,
%   sampleFrequency) returns a low-pass filter cascade structure for a
%   specified numberOfFilters, cornerFrequency, and sampleFrequency.
%
%   See:
%   https://en.wikipedia.org/wiki/Low-pass_filter
%   https://www.electronics-tutorials.ws/filter/filter_2.html

    lowPassFilterCascade.numberOfFilters = numberOfFilters;

    samplePeriod = 1 / sampleFrequency;
    cornerFrequency = cornerFrequency / sqrt(2^(1 / numberOfFilters) - 1);
    lowPassFilterCascade.coefficient = (samplePeriod / ((1 / (2 * pi * cornerFrequency)) + samplePeriod));

    lowPassFilterCascade.outputs = zeros(1, numberOfFilters);
end

