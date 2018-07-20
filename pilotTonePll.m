function [pilotTone, pilotToneInverseSin] = pilotTonePll(inputSignal, sampleFrequency, pilotToneFrequency)
%PILOTTONEPLL Pilot tone phase lock loop.
%   pilotTonePll(inputSignal, sampleFrequency, pilotToneFrequency) returns
%   phase-locked pilotTone and pilotToneInverseSin for an inputSignal with
%   a specified sampleFrequency and pilotToneFrequency.  The
%   pilotToneInverseSin can be used to generate carrier and symbol clocks
%   with of fixed frequency and phase relative to the pilot tone.
%
%   PI gain tuning method:
%   1) Create worst-case expected error.
%   2) Set Ki to zero and adjust Kp until the error response is a smooth
%   convergence to a steady-state value with low overshoot.  Ignore the
%   fixed tone within the error.
%   3) Increase Ki to remove the steady-state error without inducing
%   excessive oscillations.  Ki will typically be two orders of magnitude
%   greater than Kp.

    numberOfSamples = length(inputSignal);
    samplePeriod = 1 / sampleFrequency;
    time = 0:samplePeriod:((numberOfSamples - 1) * samplePeriod);

    % PI gains
    Kp = 250;
    Ki = 50000;
    %pilotToneFrequency = pilotToneFrequency * 1.1; % deliberate error used to tune PI gains

    % Phase detector low-pass filter
    phaseDetectorFilter = lowPassFilterCascadeInitialise(3, pilotToneFrequency, sampleFrequency);
    phaseDetectorToneFilter = phaseDetectorFilter; % create identical filter

    % Vectors populated during loop
    phaseDetector = zeros(1, numberOfSamples);
    filteredPhaseDetector = zeros(1, numberOfSamples);
    error = zeros(1, numberOfSamples);
    integratedError = zeros(1, numberOfSamples);
    pllPhase = zeros(1, numberOfSamples);
    pllOutput = zeros(1, numberOfSamples);
    phaseDetectorTone = zeros(1, numberOfSamples);
    filteredPhaseDetectorTone = zeros(1, numberOfSamples);
    pilotTone = zeros(1, numberOfSamples);
    pilotToneInverseSin = zeros(1, numberOfSamples);

    % Loop through samples
    for index = 3:numberOfSamples

        % Calculate phase detector
        phaseDetector(index - 1) = inputSignal(index - 1) * pllOutput(index - 1);

        % Low-pass filter phase detector
        phaseDetectorFilter = lowPassFilterCascadeUpdate(phaseDetectorFilter, phaseDetector(index - 1));
        filteredPhaseDetector(index - 1) = phaseDetectorFilter.outputs(end); % use output of last filter in cascade

        % Calculate error
        error(index - 1) = filteredPhaseDetector(index - 1) - filteredPhaseDetectorTone(index - 1);

        % Calculate PI feedback
        integratedError(index - 1) = integratedError(index - 2) + (error(index - 1) * samplePeriod);
        piFeedback = Kp * error(index - 1) + Ki * integratedError(index - 1);

        % Update PPL phase
        pllPhase(index) = pllPhase(index - 1) + piFeedback * (2 * pi * samplePeriod);

        % Calculate PPL output
        pllOutputInverseSin = 2 * pi * pilotToneFrequency * time(index) + pllPhase(index);
        pllOutput(index) = sin(pllOutputInverseSin);

        % Calculate phase error tone
        phaseDetectorTone(index) = -0.5 * sin(2 * pllOutputInverseSin);

        % Low-pass filter phase error tone
        phaseDetectorToneFilter = lowPassFilterCascadeUpdate(phaseDetectorToneFilter, phaseDetectorTone(index));
        filteredPhaseDetectorTone(index) = phaseDetectorToneFilter.outputs(end); % use output of last filter in cascade

        % Calculate pilot tone
        pilotToneInverseSin(index) = pllOutputInverseSin - (pi / 2);
        pilotTone(index) = sin(pilotToneInverseSin(index));
    end

    % Plot
    figure;
    axes1 = subplot(2, 1, 1);
    hold on;
    plot(inputSignal);
    plot(pilotTone);
    title('Pilot tone PLL input and output');
    legend('Input', 'Output - 90^{\circ}');
    axes2 = subplot(2, 1, 2);
    hold on;
    plot(error);
    plot([1, numberOfSamples], [0, 0], 'k');
    title('Pilot tone PLL error');
    linkaxes([axes1, axes2], 'x');

end
