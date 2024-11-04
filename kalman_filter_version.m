clear all;
clc;

% Video dosyasını yükle
videoFile = 'C:\Users\pc\Desktop\tapa_video1.mp4';  % Video dosyasının yolu
videoObj = VideoReader(videoFile);

% Enter camera height
camHeight = 100;

% Enter object initial height
objHeight = 150;

% Enter total distance between camera and wall
distanceTotal = 400;

% Enter distance between camera and wall
distanceCamWall = 50;

% Pixel to meter conversion factor
pixeltometerfac = 0.00222;

% Calculate height distance between camera and object
HeightdistancesCamObj = abs(objHeight - camHeight);

% Calculate distance factor
distanceCamWall = abs(distanceTotal - distanceCamWall);
desiredHeight = (HeightdistancesCamObj * distanceTotal) / distanceCamWall;
distanceFactor = desiredHeight / HeightdistancesCamObj;

% Read frames
frames = {};
frameIndex = 1;

while hasFrame(videoObj)
    frames{frameIndex} = readFrame(videoObj);
    frameIndex = frameIndex + 1;
end

numFrames = length(frames); % Total number of frames

% Allocate for positions
positions = zeros(numFrames, 3); % Each frame for (frame number, x, y)

% Create a figure for visualization
figure;

for k = 1:numFrames
    frame = frames{k};
    
    % Apply gamma correction to brighten the image
    gammaValue = 1.2; % Adjust gamma value as needed
    frame = imadjust(frame, [], [], gammaValue);
    
    % Apply adaptive histogram equalization to each color channel
    R = adapthisteq(frame(:,:,1), 'ClipLimit', 0.02, 'Distribution', 'rayleigh');
    G = adapthisteq(frame(:,:,2), 'ClipLimit', 0.02, 'Distribution', 'rayleigh');
    B = adapthisteq(frame(:,:,3), 'ClipLimit', 0.02, 'Distribution', 'rayleigh');
    frame = cat(3, R, G, B);
    
    % Convert to HSV color space
    hsvFrame = rgb2hsv(frame);
    
    % Create a mask for yellow color range (adjustable values)
    mask = (hsvFrame(:,:,1) >= 0.1 & hsvFrame(:,:,1) <= 0.2) & ... % Yellow hue range
           (hsvFrame(:,:,2) >= 0.4 & hsvFrame(:,:,2) <= 1.0) & ...
           (hsvFrame(:,:,3) >= 0.5 & hsvFrame(:,:,3) <= 1.0);
       
    % Fill the mask
    mask = imfill(mask, 'holes');
    
    % Circle detection (for small and medium circles)
    [centers, radii] = imfindcircles(mask, [10, 20], 'ObjectPolarity', 'bright', 'Sensitivity', 0.93, 'EdgeThreshold', 0.1);
    
    % Plot the frame and detected circles
    imshow(frame);
    hold on;
    if ~isempty(centers)
        viscircles(centers, radii, 'EdgeColor', 'r'); % Plot circles in red
        positions(k, 2:3) = centers(1, :);  % Get the center of the first detected circle
    end
    positions(k, 1) = k;  % Store frame number
    hold off;
    drawnow;
end

% Check valid positions
validPositions = any(positions(:, 2:3), 2);  % Non-zero positions
positions = positions(validPositions, :);

% Check if there are sufficient valid positions
if size(positions, 1) < 2
    error('Yeterli geçerli pozisyon bulunamadı. Çember tespiti başarısız olabilir.');
end

% Enter video fps
videofps = 30;

% Frame rate (fps)
fps = videoObj.FrameRate;

% Enter camera fps
camfps = 1000;

% Kalman filter parameters
dt = 1/fps; % Time step
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; % State transition model
H = [1 0 0 0; 0 1 0 0]; % Measurement model
Q = 0.001 * eye(4); % Process noise covariance (increased for more smoothing)
R = 5 * eye(2); % Measurement noise covariance (increased for more smoothing)
P = eye(4); % Initial estimate error covariance

% Initialize state
x = [positions(1, 2); positions(1, 3); 0; 0]; % Initial state [x; y; vx; vy]

% Allocate for filtered positions
filteredPositions = zeros(size(positions));

% Kalman filter loop
for i = 1:size(positions, 1)
    % Prediction
    x = A * x;
    P = A * P * A' + Q;
    
    % Update
    z = positions(i, 2:3)'; % Measurement
    y = z - H * x; % Measurement residual
    S = H * P * H' + R; % Residual covariance
    K = P * H' / S; % Kalman gain
    x = x + K * y; % Updated state estimate
    P = (eye(size(K, 1)) - K * H) * P; % Updated estimate error covariance
    
    % Store filtered positions
    filteredPositions(i, :) = [positions(i, 1), x(1), x(2)];
end

% Calculate distances with frame number using filtered positions
distances = sqrt(diff(filteredPositions(:,2)).^2 + diff(filteredPositions(:,3)).^2);

% Create a matrix for distances including frame number
distancesWithFrame = [filteredPositions(2:end, 1), distances];

% Calculate time intervals between frames
timeIntervals = diff(filteredPositions(:,1)) / fps; % Time interval in seconds

% fpsfactor = camfps / videofps;
fpsfactor = camfps / videofps;

% Calculate velocities
velocities = (distances ./ timeIntervals) * distanceFactor * pixeltometerfac * fpsfactor;  % Velocity (meters/second)

% Apply a more aggressive low-pass filter to smooth the velocity data
smoothVelocities = sgolayfilt(velocities, 3, 21); % Savitzky-Golay filter

% Plot velocities
figure;
plot(distancesWithFrame(:, 1), smoothVelocities, '-o');  % Plot velocities against frame numbers
xlabel('Frame Number');
ylabel('Velocity (meters/second)');
title('Object Velocity');

% Calculate average velocity (optional)
averageVelocity = mean(smoothVelocities);
disp(['Ortalama Hız: ', num2str(averageVelocity), ' meters/second']);
