%% The folder containing the insect head images

close all;
dataFolder = uigetdir(pwd, ...
'Dataset Insect Head 360');
if isequal(dataFolder, 0)
    return;
end
%% Find the image taken at 0 degrees

files = dir(fullfile(dataFolder, '*.tif'));
assert(~isempty(files), ...
'No tif images found. Select the folder containing the images.');
names = {files.name};
% Match the angle at the end of the filename.
isZero = ~cellfun('isempty', ...
regexp(names, '_\s*-\s*0\.0\.tif$', 'once'));
assert(sum(isZero) == 1, ...
'Could not identify exactly one 0-degree image.');
zeroFile = files(isZero);
%% Open image

imagePath = fullfile(dataFolder, zeroFile.name);
A = imread(imagePath);
%% Display it

figure;
imagesc(A);
colormap gray;
axis image;
xlabel('Column: left to right');
ylabel('Row: top to bottom');
title('Insect Head Image at 0 Degrees');
%% Mark the center row

hold on;
plot([1 size(A, 2)], [750 750], ...
'r-', 'LineWidth', 2);
plot(500, 750, 'r+', 'MarkerSize', 12, ...
'LineWidth', 2);
%% Print the image dimensions

fprintf('Image height: %d rows\n', size(A, 1));
fprintf('Image width: %d columns\n', size(A, 2));
%% Collect the brightness along the red line

profile = double(A(750, :));
%% Display those numbers as a graph

figure;
plot(1:1000, profile, 'b-', 'LineWidth', 1.5);
xlabel('Position along the red line (column)');
ylabel('Stored image brightness');
title('Brightness Along Row 750');
grid on;
%% Read the angle from each filename

angles = zeros(1, numel(files));
for k = 1:numel(files)
    token = regexp(files(k).name, ...
    '_\s*([-\d\.]+)B?\.tif$', 'tokens', 'once');
    assert(~isempty(token), ...
    'Cannot read angle from filename: %s', files(k).name);
    angles(k) = str2double(token{1});
end
%% Put images in angle order

[angles, order] = sort(angles);
files = files(order);
% 360 degrees repeats the starting viewing direction.
keep = angles >= 0 & angles < 360;
angles = angles(keep);
files = files(keep);
assert(isequal(angles, 0:2:358), ...
'Expected images at 0, 2, 4, ... 358 degrees.');
%% Collect row 750 from every image

rawSinogram = zeros(1000, numel(files));
for k = 1:numel(files)
    imagePath = fullfile(dataFolder, files(k).name);
    currentImage = imread(imagePath);
    rawSinogram(:, k) = double(currentImage(750, :)).';
end
%% Display all the profiles together

figure;
imagesc(angles, 1:1000, rawSinogram);
colormap gray;
colorbar;
xlabel('Viewing angle (degrees)');
ylabel('Position along the selected row (column)');
title('Brightness Sinogram: Row 750');
fprintf('Collected %d viewing angles.\n', numel(angles));
%% Save the raw sinogram and its angles

save('Insecthead_sinogram.mat', ...
'rawSinogram', 'angles');
%% Save the displayed figure

saveas(gcf, 'InsectHead_sinogram.png');
%% Enter the scanner dimensions

sourceToCenter = 105.8;       % Distance in millimeters
centerToDetector = 13.6;      % Distance in millimeters
detectorPixelSize = 0.02;     % Width of one detector pixel, mm
sourceToDetector = sourceToCenter + centerToDetector;
%% Each detector position a distance from the axis

axisColumn = 500;
numberOfColumns = size(rawSinogram, 1);
detectorPosition = ...
((1:numberOfColumns)' - axisColumn) * detectorPixelSize;
% Map those positions to an equivalent line through the rock center.
u = detectorPosition * sourceToCenter / sourceToDetector;
%% Calculate the fan-beam correction

fanWeight = sourceToCenter ./ ...
sqrt(sourceToCenter^2 + u.^2);
%% Checking the Correction

figure;
plot(u, fanWeight, 'b-', 'LineWidth', 2);
grid on;
xlabel('Distance from the Center on the Equivalent Detector (mm)');
ylabel('Fan-beam correction weight');
title('Correction for the spreading X-ray rays');
%% Working copy of the measurements

projectionData = double(rawSinogram);
%% Applying the fan-beam correction at every angle

weightedSinogram = zeros(size(projectionData));
for k = 1:numel(angles)
    weightedSinogram(:, k) = ...
    projectionData(:, k) .* fanWeight;
end
%% Compare one profile before and after correction

figure;
plot(u, projectionData(:, 1), ...
'k-', 'LineWidth', 1.5);
hold on;
plot(u, weightedSinogram(:, 1), ...
'r--', 'LineWidth', 1.5);
xlabel('Equivalent detector position (mm)');
ylabel('Projection value');
title('First Viewing Angle: Before and After Correction');
legend('Before Correction', 'After Correction');
grid on;
%% Finding the spacing between detector samples

du = u(2) - u(1);     % Spacing in millimeters
%% Choosing a padded length for filtering

numberOfSamples = size(weightedSinogram, 1);
nFFT = 2^nextpow2(2 * numberOfSamples);
%% Spatial-frequency coordinates

frequency = (0:nFFT-1)' / (nFFT * du);
% FFT stores negative frequencies in the second half.
negativePart = frequency > 1 / (2 * du);
frequency(negativePart) = ...
frequency(negativePart) - 1 / du;
%% Ramp filter

rampFilter = abs(frequency);
%% Display the ramp filter 

[sortedFrequency, order] = sort(frequency);
figure;
plot(sortedFrequency, rampFilter(order), ...
'b-', 'LineWidth', 2);
xlabel('Spatial frequency (cycles per mm)');
ylabel('Filter weight');
title('Ramp Filter');
grid on;
%% Space for the filtered profiles

filteredSinogram = zeros(size(weightedSinogram));
%% Filter one viewing angle at a time

for k = 1:numel(angles)
    % Get the profile for this angle.
    currentProfile = weightedSinogram(:, k);
    % Express the profile in terms of spatial frequencies.
    profileFFT = fft(currentProfile, nFFT);
    % Apply the ramp weights we created.
    filteredFFT = profileFFT .* rampFilter;
    % Convert back to values along the detector.
    filteredProfile = real(ifft(filteredFFT));
    % Keep the original samples; discard the extra padding.
    filteredSinogram(:, k) = ...
    filteredProfile(1:numberOfSamples);
end
%% Compare the first profile before and after filtering

figure;
subplot(2, 1, 1);
plot(u, weightedSinogram(:, 1), 'k-');
xlabel('Equivalent detector position (mm)');
ylabel('Before filtering');
title('First viewing angle');
grid on;
subplot(2, 1, 2);
plot(u, filteredSinogram(:, 1), 'b-');
xlabel('Equivalent detector position (mm)');
ylabel('After filtering');
grid on;
%% Reconstruction grid

imageSize = 350;         % 350 rows and 350 columns
outputPixelSize = 0.04;  % Millimeters per output pixel
%% Give each pixel a physical position

coordinates = ...
((1:imageSize) - (imageSize + 1)/2) * outputPixelSize;
[X, Y] = meshgrid(coordinates, coordinates);
%% Select the first viewing angle

viewNumber = 1;
beta = angles(viewNumber) * pi / 180;  % Degrees to radians
R = sourceToCenter;
%% Locate each pixel relative to the source direction

radialPosition = X*cos(beta) + Y*sin(beta);
sidewaysPosition = -X*sin(beta) + Y*cos(beta);
%% Find the detector ray passing through each pixel

detectorLookup = ...
R * sidewaysPosition ./ (R - radialPosition);
%% Read the filtered measurement for each ray

rayValues = interp1( ...
u, filteredSinogram(:, viewNumber), ...
detectorLookup, 'linear', 0);
%% Apply the fan-beam backprojection weight

backprojectionWeight = ...
(R ./ (R - radialPosition)).^2;
singleView = rayValues .* backprojectionWeight;
%% Start a fresh reconstruction

reconstruction = zeros(imageSize, imageSize);
R = sourceToCenter;
numberOfViews = numel(angles);
rotationDirection = -1;
%% Backprojection every viewing angle

for k = 1:numberOfViews
    beta = rotationDirection * angles(k) * pi / 180;
    % Locate pixels relative to this viewing direction.
    radialPosition = X*cos(beta) + Y*sin(beta);
    sidewaysPosition = -X*sin(beta) + Y*cos(beta);
    % Find the detector ray through each pixel.
    detectorLookup = ...
    R * sidewaysPosition ./ (R - radialPosition);
    % Look up the filtered measurements.
    rayValues = interp1( ...
    u, filteredSinogram(:, k), ...
    detectorLookup, 'linear', 0);
    % Apply the backprojection geometry weight.
    backprojectionWeight = ...
    (R ./ (R - radialPosition)).^2;
    % Add this angle's contribution to the running image.
    reconstruction = reconstruction + ...
    rayValues .* backprojectionWeight;
end
%% Scale the sum for the full 360-degree scan

reconstruction = reconstruction * pi / numberOfViews;
%% Display the full reconstruction

figure;
imagesc(coordinates, coordinates, reconstruction);
axis image;
axis xy;
colormap gray;
colorbar;
xlabel('x (mm)');
ylabel('y (mm)');
title('Insect Head Slice: Trial FBP using Stored tif Values');
%% Save 

save('InsectHead_first_reconstruction.mat', ...
'reconstruction', 'coordinates', 'angles', ...
'sourceToCenter', 'sourceToDetector', ...
'detectorPixelSize', 'axisColumn', ...
'outputPixelSize', 'rotationDirection');
saveas(gcf, 'InsectHead_first_reconstruction.png');