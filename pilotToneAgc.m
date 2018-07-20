function output = pilotToneAgc(input, sampleFrequency, pilotToneFrequency)
%PILOTTONEAGC Automatic gain control using pilot tone.
%   pilotToneAgc(inputSignal, sampleFrequency, pilotToneFrequency) returns
%   the normalised outputSignal for an inputSignal with a specified
%   sampleFrequency and pilotToneFrequency.

    numberOfSamples = length(input);
    samplePeriod = 1 / sampleFrequency;
    pilotToneSamplesPerCycle = ceil(sampleFrequency / pilotToneFrequency);

    % Feedback gain
    feedbackGain = 1000 * samplePeriod;

    % Band-pass filter
    filterOrder = 6;
    outputLowPassFilter = lowPassFilterCascadeInitialise(filterOrder, pilotToneFrequency, sampleFrequency);
    outputHighPassFilter = highPassFilterCascadeInitialise(filterOrder, pilotToneFrequency, sampleFrequency);

    % Determine band-pass filter gain at pilot tone frequency
    referencePilotToneNumberOfSamples = 100 * pilotToneSamplesPerCycle; % 100 cycles
    referencePilotToneTime = 0:samplePeriod:((referencePilotToneNumberOfSamples - 1) * samplePeriod);
    referencePilotTone = sin(2 * pi * pilotToneFrequency * referencePilotToneTime);
    referencePilotToneLowPassFilter = outputLowPassFilter; % create identical filter
    referencePilotToneHighPassFilter = outputHighPassFilter; % create identical filter
    for index = 2:referencePilotToneNumberOfSamples
        referencePilotToneLowPassFilter = lowPassFilterCascadeUpdate(referencePilotToneLowPassFilter, referencePilotTone(index));
        referencePilotTone(index) = referencePilotToneLowPassFilter.outputs(end); % use output of last filter in cascade
        referencePilotToneHighPassFilter = highPassFilterCascadeUpdate(referencePilotToneHighPassFilter, referencePilotTone(index));
        referencePilotTone(index) = referencePilotToneHighPassFilter.outputs(end); % use output of last filter in cascade
    end
    referencePilotToneAmplitude = max(referencePilotTone((end - pilotToneSamplesPerCycle):end)); % amplitude of final cycle

    % Vectors populated during loop
    output = zeros(1, numberOfSamples);
    filteredOutput = zeros(1, numberOfSamples);
    envelope = zeros(1, numberOfSamples);
    error = zeros(1, numberOfSamples);
    gain = ones(1, numberOfSamples);

    % Loop through samples
    for index = 2:numberOfSamples

        % Apply gain
        output(index) = gain(index - 1) * input(index);

        % Low-pass filter output signal
        outputLowPassFilter = lowPassFilterCascadeUpdate(outputLowPassFilter, output(index));
        filteredOutput(index) = outputLowPassFilter.outputs(end); % use output of last filter in cascade
        outputHighPassFilter = highPassFilterCascadeUpdate(outputHighPassFilter, filteredOutput(index));
        filteredOutput(index) = outputHighPassFilter.outputs(end); % use output of last filter in cascade

        if index > pilotToneSamplesPerCycle

            % Envelope follower
            envelope(index) = max(abs(filteredOutput((index - pilotToneSamplesPerCycle):index)));         

            % Calculate error
            error(index) = envelope(index) - referencePilotToneAmplitude;

            % Update gain
            gain(index) = gain(index - 1) - (feedbackGain * error(index));
        end

    end

    % Plot
    figure;
    axes1 = subplot(2, 1, 1);
    hold on;
    plot(filteredOutput);
    plot(envelope, 'r');
    plot(-1 * envelope, 'r');
    title('Pilot tone AGC amplitude');
    legend('Pilot tone', 'Pilot tone envelope');
    axes2 = subplot(2, 1, 2);
    hold on;
    plot(error);
    plot([1, numberOfSamples], [0, 0], 'k');
    title('Pilot tone AGC error');
    linkaxes([axes1, axes2], 'x');

end
