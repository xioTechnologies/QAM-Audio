function highPassFilterCascade = highPassFilterCascadeInitialise(numberOfFilters, cornerFrequency, sampleFrequency)
%HIGHPASSFILTERCASCADEINITIALISE Initialises high-pass filter cascade.
%   highPassFilterCascadeInitialise(numberOfFilters, cornerFrequency,
%   sampleFrequency) returns a high-pass filter cascade structure for a
%   specified numberOfFilters, cornerFrequency, and sampleFrequency.
%
%   See:
%   https://en.wikipedia.org/wiki/High-pass_filter
%   https://www.electronics-tutorials.ws/filter/filter_2.html

    highPassFilterCascade.numberOfFilters = numberOfFilters;

    cornerFrequency = cornerFrequency * sqrt(2^(1 / numberOfFilters) - 1);
    highPassFilterCascade.coefficient = 1 / ((2 * pi * cornerFrequency * (1 / sampleFrequency)) + 1);

    highPassFilterCascade.inputs = zeros(1, numberOfFilters);
    highPassFilterCascade.outputs = zeros(1, numberOfFilters);
end

