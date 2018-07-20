function distributeFigures(overlap)
%DISTRIBUTEFIGURES Distributes figures across the screen.
%   distributeFigures(overlap), If overlap is true then the figure windows
%   will overlap to hide the figure toolbars and maximise the figure size.

    % Get figure handles
    figureHandlesUnsorted = get(groot, 'Children');
    numberOfFigures = length(figureHandlesUnsorted);

    % Sort figure handles
    [~, I] = sort([figureHandlesUnsorted.Number]);
    figureHandles(1:numberOfFigures) = figureHandlesUnsorted(I);

    % Calculate number of rows and columns
    numberOfRows = 1;
    numberOfColumns = 1;
    while true

        % Special case for square grid
        squareRoot = sqrt(numberOfFigures);
        if mod(squareRoot, 1) == 0
            numberOfRows = squareRoot;
            numberOfColumns = squareRoot;
            break;
        end

        % Break if current number of rows and columns sufficient
        if numberOfFigures <= (numberOfRows * numberOfColumns)
            break;
        end

        % Increment grid size while maintaining aspect ratio
        if(numberOfColumns <= (numberOfRows + 1))
            numberOfColumns = numberOfColumns + 1;
        else
            numberOfRows = numberOfRows + 1;
        end
    end

    % Get screen size
    screenSize = get( groot, 'ScreenSize');
    screenWidth = screenSize(3);
    screenHeight = screenSize(4);

    % Calculate dimensions of each figure
    figureWidth = screenWidth / numberOfColumns;
    figureHeight = screenHeight / numberOfRows;

    % Set position and size of each figure
    for index = 1:length(figureHandles)

        % Calculate row and number number
        rowNumber = ceil(index / numberOfColumns);
        columnNumber = mod(index, numberOfColumns);
        if columnNumber == 0
            columnNumber = numberOfColumns;
        end

        % Calculate figure position
        x = 1 + ((columnNumber - 1) * figureWidth);
        y = 1 + ((rowNumber - 1) * figureHeight);

        % Set figure position and size
        positionString = 'OuterPosition';
        if overlap == true
             positionString = 'Position';
        end
        set(figureHandles(index), positionString, [x, y, figureWidth, figureHeight]);

        % Restore figure in case it was minimised
        figure(figureHandles(index));
    end

end
