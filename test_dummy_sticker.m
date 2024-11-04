clear all;
clc;

% Load Video File
videoFile = 'C:\Users\pc\Desktop\video.mp4';  % video file path
videoObj = VideoReader(videoFile);

% Pixel to meter conversion factor
pixeltometerfac =1/1172.225; % better results sets 0.001 to 0.0011617775196
camfps = 3285; % Enter camera fps

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
           (hsvFrame(:,:,2) >= 0.3 & hsvFrame(:,:,2) <= 1.0) & ...
           (hsvFrame(:,:,3) >= 0.3 & hsvFrame(:,:,3) <= 1.0);
       
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

positions(:,4) = positions(:,1)./camfps;
size_positions = size(positions);

ortalama=12;

for i=1:fix(size_positions(1)/ortalama)-1
    positions_new(i,:) = mean(positions((i-1)*ortalama+1:(i-1)*ortalama+ortalama,:));
end

positions_new(i+1,:)=mean(positions((i)*ortalama+1:(i)*ortalama+mod(size_positions(1),ortalama),:));

 % Check if there are sufficient valid positions
 if size(positions, 1) < 2
    error('Valid positions is not enough.Circle detection is failed.');
 end

% Calculate velocities
distances = sqrt(diff(positions_new(:,2)).^2 + diff(positions_new(:,3)).^2);
timeIntervals = diff(positions_new(:,4)); % Time interval in seconds
velocities = (distances .* pixeltometerfac ./ timeIntervals);

% Plot velocity-time graph
figure;
plot(positions_new(2:end,4), velocities, '-o');
title('Velocity-Time Graph');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

% Plot position-time graph (using y-position)
figure;
plot(positions_new(:,4), positions_new(:,3) * pixeltometerfac, '-o');
title('Position-Time Graph');
xlabel('Time (s)');
ylabel('Position (m)');
