function scaledConstellationPoints = constellationScaleCorrection(constellationPoints, bitsPerSymbol)
%CONSTELLATIONSCALECORRECTION Normalises the constellation scale.
%   scaledConstellationPoints = constellationScaleCorrection(
%   constellationPoints, bitsPerSymbol) returns the constellationPoints as
%   scaledConstellationPoints for a specified bitsPerSymbol.
%   constellationPoints must contain at least one constellation point at
%   the each of the constellation within each window.

    numberOfConstellationPoints = length(constellationPoints);
    maximumConstellationPointValue = sqrt(2^bitsPerSymbol) - 1;

    % Adjustable parameters
    gain = 0.01;
    windowsSize = 500;

    % Vectors populated during loop
    scale = ones(1, numberOfConstellationPoints);
    error = zeros(1, numberOfConstellationPoints);

    % Variables to track maximum symbol value in window
    maximumConstellationPointValueInWindow = 0;
    constellationPointCounter = 0;

    % Loop through symbols
    for index = 2:numberOfConstellationPoints

        % Adjust scale
        constellationPoints(index, :) = scale(index - 1) * constellationPoints(index, :);

        % Find nearest symbol (round to odd grid values)
        nearestSymbol = floor(constellationPoints(index, :));
        if mod(nearestSymbol(1), 2) == 0
            nearestSymbol(1) = nearestSymbol(1) + 1;
        end
        if mod(nearestSymbol(2), 2) == 0
            nearestSymbol(2) = nearestSymbol(2) + 1;
        end

        % Limit maximum symbol value
        if nearestSymbol(1) > maximumConstellationPointValue
            nearestSymbol(1) = maximumConstellationPointValue;
        end
        if nearestSymbol(2) > maximumConstellationPointValue
            nearestSymbol(2) = maximumConstellationPointValue;
        end

        % Adjust scale for each window
        maximumConstellationPointValueInWindow = max([maximumConstellationPointValueInWindow, nearestSymbol(1), nearestSymbol(2)]);
        constellationPointCounter = constellationPointCounter + 1;
        if constellationPointCounter >= windowsSize
            scale(index - 1) = scale(index - 1) * (maximumConstellationPointValue / maximumConstellationPointValueInWindow);
            maximumConstellationPointValueInWindow = 0;
            constellationPointCounter = 0;
        end

        % Calculate error
        error(index) = norm(constellationPoints(index, :)) - norm(nearestSymbol);

        % Update scale
        scale(index) = scale(index - 1) - gain * error(index);

    end

    % Discard symbols obtained prior to scale correction
    scaledConstellationPoints = constellationPoints((2 * windowsSize):end, :);

    % Plot
    figure;
    axes1 = subplot(2, 1, 1);
    hold on;
    plot(scale);
    title('Constellation scale correction');
    axes2 = subplot(2, 1, 2);
    hold on;
    plot(error);
    plot([1, numberOfConstellationPoints], [0, 0], 'k');
    title('Constellation scale error');
    linkaxes([axes1, axes2], 'x');

end
