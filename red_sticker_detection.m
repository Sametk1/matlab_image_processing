clear all;
clc;

% Video dosyasını yükle
videoFile = 'C:\Users\pc\Desktop\tapa_video1.mp4';  % İşlenecek video dosyasının yolu
videoObj = VideoReader(videoFile);  % VideoReader nesnesi ile videoyu okuma

% Pixel to meter conversion factor
pixeltometerfac = 0.0010891664246; % Piksel-metre dönüşüm faktörü. Bu, piksellerin fiziksel uzunluğa dönüştürülmesini sağlar.

% Videodaki kareleri okuma işlemi için bir dizi oluşturuyoruz
frames = {};  % Video karelerini saklayacak hücre dizisi
frameIndex = 1;  % Kare numarasını tutacak sayaç

% Videodaki tüm kareleri okuma döngüsü
while hasFrame(videoObj)
    frames{frameIndex} = readFrame(videoObj);  % Kareleri sırayla okuyor ve dizide saklıyoruz
    frameIndex = frameIndex + 1;  % Bir sonraki kare için sayaç artışı
end

numFrames = length(frames); % Toplam kare sayısını hesaplama

% Pozisyonları saklamak için bir matris oluşturuyoruz
positions = zeros(numFrames, 3); % Her kare için (kare numarası, x, y) pozisyonları saklanacak

% Görselleştirme için bir figure (grafik penceresi) oluşturuyoruz
figure;

% Her kareyi işleyip pozisyonları hesaplayan döngü
for k = 1:numFrames
    frame = frames{k};  % Mevcut kareyi alıyoruz
    
    % Gamma düzeltmesi uygulayarak görüntüyü parlaklaştırıyoruz
    gammaValue = 1.2; % Gamma değeri ayarı, parlaklık düzeyini belirler
    frame = imadjust(frame, [], [], gammaValue);  % Gamma düzeltmesi uygulama
    
    % Her renk kanalına adaptif histogram eşitleme uyguluyoruz
    R = adapthisteq(frame(:,:,1), 'ClipLimit', 0.02, 'Distribution', 'rayleigh');  % Kırmızı kanal
    G = adapthisteq(frame(:,:,2), 'ClipLimit', 0.02, 'Distribution', 'rayleigh');  % Yeşil kanal
    B = adapthisteq(frame(:,:,3), 'ClipLimit', 0.02, 'Distribution', 'rayleigh');  % Mavi kanal
    frame = cat(3, R, G, B);  % Renk kanallarını tekrar birleştiriyoruz
    
    % Görüntüyü HSV renk uzayına dönüştürüyoruz
    hsvFrame = rgb2hsv(frame);
    
    % Kırmızı renk aralığı için bir maske oluşturuyoruz (Bu aralıklar ayarlanabilir)
    mask = (hsvFrame(:,:,1) >= 0.9 | hsvFrame(:,:,1) <= 0.1) & ... % Kırmızı ton aralığı
           (hsvFrame(:,:,2) >= 0.4 & hsvFrame(:,:,2) <= 1.0) & ... % Doygunluk aralığı
           (hsvFrame(:,:,3) >= 0.5 & hsvFrame(:,:,3) <= 1.0);       % Parlaklık aralığı
       
    % Maskeyi dolduruyoruz (Delikleri kapatıyoruz)
    mask = imfill(mask, 'holes');
    
    % Küçük ve orta boyutlu çemberler için çember tespiti yapıyoruz
    [centers, radii] = imfindcircles(mask, [10, 20], 'ObjectPolarity', 'bright', 'Sensitivity', 0.93, 'EdgeThreshold', 0.1);
    
    % Mevcut kareyi ve tespit edilen çemberleri çizdiriyoruz
    imshow(frame);  % Kareyi görüntüleme
    hold on;
    if ~isempty(centers)
        viscircles(centers, radii, 'EdgeColor', 'r'); % Tespit edilen çemberleri kırmızı renkte çizdiriyoruz
        positions(k, 2:3) = centers(1, :);  % İlk tespit edilen çemberin merkezini pozisyonlara ekliyoruz
    end
    
    positions(k, 1) = k;  % Kare numarasını pozisyonlara ekliyoruz
    hold off;
    drawnow;  % Grafik penceresini güncelliyoruz
end

% Kamera fps değerini giriyoruz (Videonun karesi başına düşen saniye sayısı)
camfps = 2577;

% Geçerli pozisyonları kontrol ediyoruz (Geçerli pozisyonlar sıfırdan farklı olanlardır)
validPositions = any(positions(:, 2:3), 2);  % Geçerli pozisyonları belirleme
positions = positions(validPositions, :);  % Geçerli pozisyonları filtreleme

% Zaman ekseninde pozisyonlar
positions(:,4) = positions(:,1)./camfps;  % Kare numarasını zamana dönüştürme
size_positions = size(positions);  % Pozisyonların boyutunu alıyoruz

ortalama = 12;  % Ortalama alınacak kare sayısı (Gürültüyü azaltmak için)

% Ortalama pozisyonları hesaplıyoruz
for i = 1:fix(size_positions(1)/ortalama)-1
    positions_new(i,:) = mean(positions((i-1)*ortalama+1:(i-1)*ortalama+ortalama,:));  % Ortalama alma
end

% Son grup için ortalama pozisyon hesaplama
positions_new(i+1,:) = mean(positions((i)*ortalama+1:(i)*ortalama+mod(size_positions(1),ortalama),:));

% Yeterli geçerli pozisyon olup olmadığını kontrol ediyoruz
if size(positions, 1) < 2
    error('Yeterli geçerli pozisyon bulunamadı. Çember tespiti başarısız olabilir.');
end

% Hızları hesaplıyoruz
distances = sqrt(diff(positions_new(:,2)).^2 + diff(positions_new(:,3)).^2);  % İki ardışık pozisyon arasındaki mesafeyi hesaplama
timeIntervals = diff(positions_new(:,4)); % Zaman aralıklarını hesaplama (saniye cinsinden)
velocities = (distances .* pixeltometerfac ./ timeIntervals);  % Hızı hesaplama (mesafe/zaman)

% Hız-zaman grafiğini çiziyoruz
figure;
plot(positions_new(2:end,4), velocities, '-o');  % Zaman eksenine karşı hızları çizdirme
title('Velocity-Time Graph');  % Grafik başlığı
xlabel('Time (s)');  % X ekseni etiketi (Zaman)
ylabel('Velocity (m/s)');  % Y ekseni etiketi (Hız)

% Pozisyon-zaman grafiğini çiziyoruz (y pozisyonunu kullanarak)
figure;
plot(positions_new(:,4), positions_new(:,3) * pixeltometerfac, '-o');  % Zaman eksenine karşı pozisyonları çizdirme
title('Position-Time Graph');  % Grafik başlığı
xlabel('Time (s)');  % X ekseni etiketi (Zaman)
ylabel('Position (m)');  % Y ekseni etiketi (Pozisyon)
