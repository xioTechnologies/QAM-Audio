function constellationPoints = demodulator(qamSignal, sampleFrequency, rrcFilterImpulseResponse, pilotToneInverseSin, pilotToneCarrierRatio, pilotToneSymbolRatio)
%DEMODULATOR Demodulates QAM signal into constellation points.
%   demodulator(qamSignal, sampleFrequency, rrcFilterImpulseResponse,
%   pilotToneInverseSin, pilotToneCarrierRatio, pilotToneSymbolRatio)
%   returns the constellation points for a QAM signal with a specified
%   sampleFrequency, rrcFilterImpulseResponse, pilotToneInverseSin,
%   pilotToneCarrierRatio, and pilotToneSymbolRatio.  Carrier and symbol
%   clock phase will be corrected provided that the pilotToneInverseSin
%   frequency is fixed relative to the carrier and symbol clocks.

    numberOfSamples = length(qamSignal);
    samplePeriod = 1 / sampleFrequency;
    samplesPerRrcFilter = length(rrcFilterImpulseResponse);

    % Phase correction gains
    carrierPhaseCorrectionGain = 1000 * samplePeriod;
    symbolPhaseCorrectionGain = 10000 * samplePeriod;

    % Vectors populated during loop
    carrierPhase = zeros(1, numberOfSamples);
    iSignal = zeros(1, numberOfSamples);
    qSignal = zeros(1, numberOfSamples);
    iBaseband = zeros(1, numberOfSamples);
    qBaseband = zeros(1, numberOfSamples);
    symbolPhase = zeros(1, numberOfSamples);
    symbolClock = zeros(1, numberOfSamples);

    % Loop through samples
    for index = 2:numberOfSamples

        % Create carrier clocks
        iCarrier = cos(pilotToneCarrierRatio * (pilotToneInverseSin(index) + carrierPhase(index - 1)));
        qCarrier = -sin(pilotToneCarrierRatio * (pilotToneInverseSin(index) + carrierPhase(index - 1)));

        % Split I and Q signals
        iSignal(index) = qamSignal(index) * iCarrier;
        qSignal(index) = qamSignal(index) * qCarrier;

        % Process signals through RRC filters
        if index >= samplesPerRrcFilter
            iBaseband(index) = sum(rrcFilterImpulseResponse .* iSignal((index - samplesPerRrcFilter + 1):index));
            qBaseband(index) = sum(rrcFilterImpulseResponse .* qSignal((index - samplesPerRrcFilter + 1):index));
        end

        % Calculate carrier phase detector
        iCarrierPhaseDetector = iBaseband(index) * sign(qBaseband(index));
        qCarrierPhaseDetector = qBaseband(index) * sign(iBaseband(index));
        carrierPhaseDetector = qCarrierPhaseDetector - iCarrierPhaseDetector;

        % Update carrier phase
        carrierPhase(index) = carrierPhase(index - 1) + carrierPhaseCorrectionGain * carrierPhaseDetector * (2 * pi * samplePeriod);

        % Create symbol clock
        symbolClock(index) = sin(pilotToneSymbolRatio * (pilotToneInverseSin(index) + symbolPhase(index - 1)));

        % Calculate symbol phase detector
        iSymbolPhaseDetector = abs(iBaseband(index)) * symbolClock(index);
        qSymbolPhaseDetector = abs(qBaseband(index)) * symbolClock(index);
        symbolPhaseDetector = iSymbolPhaseDetector + qSymbolPhaseDetector;

        % Update symbol phase
        symbolPhase(index) = symbolPhase(index - 1) - symbolPhaseCorrectionGain * symbolPhaseDetector * (2 * pi * samplePeriod);
    end

    % Obtain symbol indexes as each symbol clock zero-crossing from negative to positive
    symbolIndexes = find([diff(sign(symbolClock)) == 2, 0]);

    % Discard symbols obtained prior to carrier and symbol phase correction
    symbolIndexes = symbolIndexes(round(length(symbolIndexes) / 2):end);

    % Create array of constellation points
    constellationPoints = zeros(length(symbolIndexes), 2);
    for index = 1:length(symbolIndexes)
        symbolIndex = symbolIndexes(index);
        constellationPoints(index, :) = [iBaseband(symbolIndex), qBaseband(symbolIndex)];
    end

    % Plot
    figure;
    hold on;
    plot(rad2deg(carrierPhase));
    plot(rad2deg(symbolPhase));
    plot([1, numberOfSamples], [0, 0], 'k');
    ylabel('Phase (degrees)');
    title('Carrier and symbol phase correction');
    legend('Carrier', 'Symbol');

end
