clc
clear
close all;

%% QAM parameters

sampleFrequency = 96000; % Hz
carrierFrequency = 9600; % Hz (must be divisible by sampleFrequency)
pilotToneFrequency = 1000; % Hz
carrierCyclesPerSymbol = 2.0; % increasing this value will increases spectral efficiency
bitsPerSymbol = 2; % must be an even number
symbolsPerRrcFilter = 5; % number of symbols
beta = 1.0; % 0 = best spectral efficiency but requires infinite filter length, 1 = poor spectral efficiency but best ISI for finite filter length
oversampling = 1; % multiplier (1 = no oversampling)
constellationRepetitions = 5;

fprintf('Bit rate:   %0.0f bps\n', bitsPerSymbol * (carrierFrequency / carrierCyclesPerSymbol));
fprintf('Latency:    %0.2f ms\n', 1000 * symbolsPerRrcFilter * carrierCyclesPerSymbol * (1 / carrierFrequency));
fprintf('Jitter:     %0.2f ms\n', 1000 * (1 / (carrierFrequency / carrierCyclesPerSymbol)));

samplesPerSymbol = carrierCyclesPerSymbol * (sampleFrequency / carrierFrequency);

%% RRC filter - https://en.wikipedia.org/wiki/Root-raised-cosine_filter

samplesPerRrcFilterSide = (symbolsPerRrcFilter / 2) * samplesPerSymbol;
samplesPerRrcFilter = 1 + symbolsPerRrcFilter * samplesPerSymbol;

tVector = linspace(-samplesPerRrcFilterSide, samplesPerRrcFilterSide, samplesPerRrcFilter);
T = samplesPerSymbol;

rrcFilterImpulseResponse = zeros(1, length(tVector));
for index = 1:length(tVector)
    t = tVector(index);
    if t == 0
        rrcFilterImpulseResponse(index) = (1 / T) * (1 + beta * ((4 / pi) - 1));
    elseif abs(t) == (T / (4 * beta))
        rrcFilterImpulseResponse(index) = (beta / (T * sqrt(2))) * ((1 + (2 / pi)) * sin(pi / (4 * beta)) + (1 - (2 / pi)) * cos(pi / (4 * beta)));
    else
        numerator = sin(pi * (t / T) * (1 - beta)) + 4 * beta * (t / T) * cos(pi * (t / T) * (1 + beta));
        denominator = pi * (t / T) * (1 - (4 * beta * (t / T))^2);
        rrcFilterImpulseResponse(index) = (1 / T) * (numerator / denominator);
    end
end

rrcFilterImpulseResponse = rrcFilterImpulseResponse * (1 / max(rrcFilterImpulseResponse));

figure;
hold on;
plot(rrcFilterImpulseResponse);
centreIndex = samplesPerRrcFilterSide + 1;
plot([centreIndex centreIndex], [min(rrcFilterImpulseResponse), max(rrcFilterImpulseResponse)], 'k--');
for symbolCount = 1:1:floor(symbolsPerRrcFilter / 2)
    index = centreIndex - (symbolCount * samplesPerSymbol);
    plot([index index], [min(rrcFilterImpulseResponse), max(rrcFilterImpulseResponse)], 'k--');
    index = centreIndex + (symbolCount * samplesPerSymbol);
    plot([index index], [min(rrcFilterImpulseResponse), max(rrcFilterImpulseResponse)], 'k--');
end
plot([0 length(tVector)], [0 0], 'k');
title('RRC filter impulse response');
legend('Impulse response', 'Symbol intervals');

%% Symbol modulation

amplitudes = linspace(-1, 1, sqrt(2^bitsPerSymbol));

iSequence = [];
qSequence = [];
for constellationCount = 1:(1 + constellationRepetitions)
    for iSequenceIndex = 1:length(amplitudes)
        for qSequenceIndex = 1:length(amplitudes)
            iSequence = [iSequence, amplitudes(iSequenceIndex)];
            qSequence = [qSequence, amplitudes(qSequenceIndex)];
        end
    end
end

randomOrderIndexes = randperm(length(iSequence));
iSequence = iSequence(randomOrderIndexes);
qSequence = qSequence(randomOrderIndexes);

numberOfSymbols = length(iSequence);

iSequence = [iSequence, zeros(1, symbolsPerRrcFilter)];
qSequence = [qSequence, zeros(1, symbolsPerRrcFilter)];

iImpulses = zeros(1, length(iSequence) * samplesPerSymbol);
qImpulses = zeros(1, length(iSequence) * samplesPerSymbol);
for sequenceIndex = 1:numberOfSymbols
    symbolIndex = (sequenceIndex - 1) * samplesPerSymbol + 1;
    iImpulses(symbolIndex) = iSequence(sequenceIndex);
    qImpulses(symbolIndex) = qSequence(sequenceIndex);
end

iBaseband = filter(rrcFilterImpulseResponse, 1, iImpulses);
qBaseband = filter(rrcFilterImpulseResponse, 1, qImpulses);

figure;
subplot(2, 1, 1);
hold on;
title('Symbol impulses');
stem(find(iImpulses), iImpulses(iImpulses ~= 0));
stem(find(qImpulses), qImpulses(qImpulses ~= 0));
legend('I', 'Q');
subplot(2, 1, 2);
hold on;
title('TX RRC filter output');
plot(iBaseband);
plot(qBaseband);
txSymbolIndexes = (samplesPerRrcFilterSide + 1) + (0:samplesPerSymbol:((numberOfSymbols - 1) * samplesPerSymbol));
for txSampleIndex = txSymbolIndexes
    plot([txSampleIndex, txSampleIndex], [amplitudes(1), amplitudes(end)], 'k--');
end
for amplitudeIndex = 1:length(amplitudes)
    plot([0 length(iBaseband)], [amplitudes(amplitudeIndex), amplitudes(amplitudeIndex)], 'k--');
end
legend('I', 'Q', 'Symbol');

%% QAM signal

samplePeriod = 1 / sampleFrequency;
time = 0:samplePeriod:((length(iBaseband) * samplePeriod) - samplePeriod);
iCarrier = cos(2 * pi * carrierFrequency * time);
qCarrier = -sin(2 * pi * carrierFrequency * time);

iSignal = iBaseband .* iCarrier;
qSignal = qBaseband .* qCarrier;

qamSignal = iSignal + qSignal;

%% Pilot tone

% pilotTone = sin(2 * pi * pilotToneFrequency * time);
% qamSignal = qamSignal + pilotTone;

%% Audio loopback

% outputGain = 0.05;
% startStopDelay = 0.2; % delay between sound start/stop and record start/stop
%
% audiorecorderObject = audiorecorder(sampleFrequency, 24, 1); % 24-bit, single channel
%
% sound(outputGain * qamSignal, sampleFrequency);
% pause(startStopDelay);
% record(audiorecorderObject);
% pause((length(qamSignal) * samplePeriod) - startStopDelay);
% stop(audiorecorderObject);
%
% qamSignal = getaudiodata(audiorecorderObject);
% time = 0:samplePeriod:((length(qamSignal) * samplePeriod) - samplePeriod);
%
% save('audioLoopback.mat');
% load('audioLoopback.mat'); % comment out all lines above this to use most recent audio

%% Plot QAM signal

figure;
plot(qamSignal);
title('QAM signal');

%% Frequency spectrum - https://uk.mathworks.com/help/matlab/ref/fft.html

% startFrequency = 0; % Hz
% stopFrquency = 24000; % Hz
% resolutionBandwidth = 100; % Hz
%
% numberOfSamples = length(qamSignal);
% fftComplex = fft(qamSignal);
% fftReal = abs(fftComplex / numberOfSamples);
% fftReal = fftReal(1:((numberOfSamples / 2) + 1)); % exclude reflection
% fftReal(2:(end - 1)) = 2 * fftReal(2:(end - 1));
% fftFrequencies = sampleFrequency * (0:(numberOfSamples / 2)) / numberOfSamples;
%
% binEdges = startFrequency:resolutionBandwidth:stopFrquency;
% binCentres = binEdges(2:end) - (resolutionBandwidth / 2);
% numberOfBins = length(binCentres);
% binAverage = zeros(1, numberOfBins);
% for binIndex = 2:numberOfBins
%     selectedIndexes = (fftFrequencies >= binEdges(binIndex - 1)) & (fftFrequencies < binEdges(binIndex));
%     numberOfSelectedIndexes = length(selectedIndexes);
%     binAverage(binIndex) = sum(fftReal(selectedIndexes)) / numberOfSelectedIndexes;
% end
%
% figure;
% hold on;
% plot(binCentres / 1000, binAverage * (1 / max(binAverage)));
% xlim([0, stopFrquency / 1000]);
% title('Frequency spectrum');
% xlabel('Frequency (kHz)');

%% AGC using pilot tone

% qamSignal = pilotToneAgc(qamSignal, sampleFrequency, pilotToneFrequency);

%% Pilot tone PPL

% [~, pilotToneInverseSin] = pilotTonePll(qamSignal, sampleFrequency, pilotToneFrequency);

%% Oversampling

oversampledSampleFrequency = oversampling * sampleFrequency;
oversampledTime = linspace(time(1), time(end), oversampling * length(time));
oversampledQamSignal = interp1(time, qamSignal, oversampledTime, 'linear');
oversampledSamplesPerFilter = oversampling * samplesPerRrcFilter;
oversampledSamplesPerSymbol = oversampling * samplesPerSymbol;
oversampledFilterImpulseResponse = interp1(1:samplesPerRrcFilter, rrcFilterImpulseResponse, linspace(1, samplesPerRrcFilter, oversampling * samplesPerRrcFilter), 'linear');
% oversampledPilotToneInverseSin = interp1(time, pilotToneInverseSin, oversampledTime, 'linear');

%% Demodulation

% pilotToneCarrierRatio = (carrierFrequency / pilotToneFrequency);
% pilotToneSymbolRatio = (carrierFrequency / carrierCyclesPerSymbol / pilotToneFrequency);
% constellationPoints = demodulator(oversampledQamSignal, oversampledSampleFrequency, oversampledFilterImpulseResponse, oversampledPilotToneInverseSin, pilotToneCarrierRatio, pilotToneSymbolRatio);

%% Constellation scale correction

% constellationPoints = constellationScaleCorrection(constellationPoints, bitsPerSymbol);
%
% figure;
% hold on;
% plot(constellationPoints(:, 1), constellationPoints(:, 2), '.');
% axis square;
% grid on;
% title('IQ constellation');
% xlabel('I');
% ylabel('Q');
% axisLimit = sqrt(2^bitsPerSymbol);
% xticks(-axisLimit:2:axisLimit);
% yticks(-axisLimit:2:axisLimit);
% xlim([-axisLimit, axisLimit]);
% ylim([-axisLimit, axisLimit]);
% set(gca, 'Yticklabel', [])
% set(gca, 'Xticklabel', []);

%% Cheat demodulation

recoveredICarrier = cos(2 * pi * carrierFrequency * oversampledTime);
recoveredQCarrier = -sin(2 * pi * carrierFrequency * oversampledTime);
rxSymbolIndexes = oversampledSamplesPerFilter + (0:oversampledSamplesPerSymbol:((numberOfSymbols - 1) * oversampledSamplesPerSymbol));

recoveredISignal = oversampledQamSignal .* recoveredICarrier;
recoveredQSignal = oversampledQamSignal .* recoveredQCarrier;

recoveredIBaseband = filter(oversampledFilterImpulseResponse, 1, recoveredISignal);
recoveredQBaseband = filter(oversampledFilterImpulseResponse, 1, recoveredQSignal);

constellationPoints = zeros(length(rxSymbolIndexes), 2);
for index = 1:length(rxSymbolIndexes)
    symbolIndex = rxSymbolIndexes(index);
    constellationPoints(index, :) = [recoveredIBaseband(symbolIndex), recoveredQBaseband(symbolIndex)];
end

figure;
subplot(2, 1, 1);
hold on;
title('Recovered I and Q signals');
plot(recoveredISignal);
plot(recoveredQSignal);
legend('I', 'Q');
subplot(2, 1, 2);
hold on;
title('RX RRC filter output');
plot(recoveredIBaseband);
plot(recoveredQBaseband);
for rxSampleIndex = rxSymbolIndexes
    plot([rxSampleIndex rxSampleIndex], [min([recoveredIBaseband recoveredQBaseband]), max([recoveredIBaseband recoveredQBaseband])], 'k--');
end
legend('I', 'Q', 'Symbol');

figure;
subplot(2, 1, 1);
hold on;
for rxSymbolIndexesIndex = 2:2:(length(rxSymbolIndexes) - 1)
    rxSymbolStartIndex = rxSymbolIndexes(rxSymbolIndexesIndex) - oversampledSamplesPerSymbol;
    rxSymbolSEndIndex = rxSymbolIndexes(rxSymbolIndexesIndex) + oversampledSamplesPerSymbol;
    plot(recoveredIBaseband(rxSymbolStartIndex:rxSymbolSEndIndex), 'b');
end
title('I eye diagram');
subplot(2, 1, 2);
hold on;
for rxSymbolIndexesIndex = 2:2:(length(rxSymbolIndexes) - 1)
    rxSymbolStartIndex = rxSymbolIndexes(rxSymbolIndexesIndex) - oversampledSamplesPerSymbol;
    rxSymbolSEndIndex = rxSymbolIndexes(rxSymbolIndexesIndex) + oversampledSamplesPerSymbol;
    plot(recoveredQBaseband(rxSymbolStartIndex:rxSymbolSEndIndex), 'b');
end
title('Q eye diagram');

figure;
hold on;
plot(constellationPoints(:, 1), constellationPoints(:, 2), '.');
axis square;
axis equal;
title('IQ constellation');

%% Distribute figures across screen

distributeFigures(false);
