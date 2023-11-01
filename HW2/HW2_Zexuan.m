%% HW2 Adv topics in MRI
clear; close all; clc; warning off;
%% Problem 1 EPG
T1 = [200:100:1500]; % ms
T2 = [50:30:300]; % ms
T1_temp = repmat(T1,length(T2),1); 
T1T2 = [T1_temp(:) repmat(T2',length(T1),1)];
% chosen T1: 200 500 800 1100 1400
% chosen T2: 50 110 170 230 290
% indices for chosen T1 and T2 combinations in T1T2: 
chosen_T1T2 = [1 30 59 88 117];

TE = 5; % ms
alpha_degree = [180 120 60]; % refocusing pulse (degree)
alpha = alpha_degree/180*pi; % refocusing pulse (radian)
alpha_1_degree = [180 120 60]+(90-[180 120 60]/2); % first refocusing pulse (degree)
alpha_1 = ([180 120 60]+(90-[180 120 60]/2))*pi/180; % first refocusing pulse (radian)
rf_pulse_direct_x = pi/2; % rf pulse direction - x
rf_pulse_direct_y = 0; % rf pulse direction - y
echo_length = 64; 

S = NaN(length(T1T2),echo_length,length(alpha)); % 126 (T1T2 combo) x 64 (# of echos) x 3 (# of alphas)
for l = 1:length(alpha)
    figure % plot signal
    count = 0; % figure subplot count
    for k = 1:length(T1T2)
        Q = [0 0 1]'; % initial FpFmZ
        Q = epg_rf(Q,90*pi/180,rf_pulse_direct_x); % apply 90-degree rf
        Q = epg_grelax(Q,T1T2(k,1),T1T2(k,2),TE/2); % relaxation until the 1st refocusing pulse
        
        for j = 1:echo_length % 64 echos in one TR
            
            % apply refocusing pulse
            if j == 1
                Q = epg_rf(Q,alpha_1(l),rf_pulse_direct_y); % apply special rf pulse on 1st refocusing pulse
            else
                Q = epg_rf(Q,alpha(l),rf_pulse_direct_y); % apply alpha pulse
            end

            % relaxation & record echo
            Q = epg_grelax(Q,T1T2(k,1),T1T2(k,2),TE/2); % echo after 2.5 ms
            S(k,j,l) = Q(1,1);

            % relaxation until next refocusing pulse
            if j < 64 % no need to continue after the last echo
                Q = epg_grelax(Q,T1T2(k,1),T1T2(k,2),TE/2); % relax for another 2.5 ms
            end
        end
        
        if ismember(k,chosen_T1T2)
            count = count + 1;
            subplot(length(chosen_T1T2),1,count)
            plot(abs(S(k,:,l)))
            ylim([0 1])
            title(['alpha=',num2str(alpha_1_degree(l)),'^o, T1=',num2str(T1T2(chosen_T1T2(count),1)),'ms, T2=',num2str(T1T2(chosen_T1T2(count),2)),'ms'])
            xlabel('Echo index')
            ylabel('M_x_y')
        end
    end
end

% Plot contour plots
chosen_echo = [6 16 32 48];
sig = abs(S(:,chosen_echo,:));
[X,Y] = meshgrid(T1,T2); % X - T1, Y - T2
for j = 1:length(alpha) % 3 alpha
    figure
    for k = 1:length(chosen_echo) % 4 chosen echo
        subplot(2,2,k)
        contourf(X,Y,reshape(sig(:,k,j),size(X)))
        xlabel('T1 (ms)'); ylabel('T2 (ms)'); colorbar
        title(['Echo#=',num2str(chosen_echo(k)),', alpha=',num2str(alpha_degree(j)),'^o'])
    end
end

%% Problem 2a: Single spin echo
load brain_maps.mat
figure
% T1w
TR = 800; % ms
TE = 5; % ms
img_T1w = NaN(size(T1map));
for j = 1:numel(img_T1w)
    img_T1w(j) = sesignal(T1map(j),T2map(j),TE,TR,0);
    img_T1w(j) = img_T1w(j)*M0map(j);
end
subplot(1,3,1); imshow(abs(img_T1w));
title('T1-weighted image')

% T2w
TR = 8000; % ms
TE = 50; % ms
img_T2w = NaN(size(T2map));
for j = 1:numel(img_T2w)
    img_T2w(j) = sesignal(T1map(j),T2map(j),TE,TR,0);
    img_T2w(j) = img_T2w(j)*M0map(j);
end
subplot(1,3,2); imshow(abs(img_T2w));
title('T2-weighted image')

% PDw
TR = 8000; % ms
TE = 5; % ms
img_PDw = NaN(size(T2map));
for j = 1:numel(img_T2w)
    img_PDw(j) = sesignal(T1map(j),T2map(j),TE,TR,0);
    img_PDw(j) = img_PDw(j)*M0map(j);
end
subplot(1,3,3); imshow(abs(img_PDw));
title('Proton density weighted image')

%% Problem 2b: Fast spin echo
TR = 3000; % ms
TE = 5; % ms
ETL = 32;
ESP = 5; % ms

% 2bi Simulate 5 TRs of FSE sequence
nTRs = 5; % five TRs
T1T2 = [1000 50; 1000 100; 2000 50; 2000 100];
figure
for j = 1:height(T1T2)
    Msignal = fsesignal(T1T2(j,1),T1T2(j,2),TE,TR,0,ETL,nTRs);
    plot(abs(Msignal((end-ETL+1):end))) % Plot the last TR
    xlim([1 ETL])
    xlabel('Echo index')
    ylabel('Signal')
    hold on
end
legend('T1=1000ms,T2=50ms','T1=1000ms,T2=100ms','T1=2000ms,T2=50ms','T1=2000ms,T2=100ms')

% 2bii Simulate FSE to create an image
[X,Y] = meshgrid([1:256],[1:256]);

TE_eff = 80; % ms
TEcho = linspace(ESP,ESP*ETL,ETL); % Varying TE
all_kspace = NaN([size(T1map),length(TEcho)]); % preallocate kspace for all TE
all_img = NaN([size(T1map),length(TEcho)]); % preallocate images for all TE
for j = 1:width(T1map)
    for k = 1:height(T1map)
        Msignal = fsesignal(T1map(k,j),T2map(k,j),ESP,TR,0,ETL,1);
        Msignal(isnan(Msignal)) = 0;
        all_img(k,j,:) = Msignal;
    end
end
all_img = repmat(M0map,[1 1 length(TEcho)]).*all_img;
for j = 1:length(TEcho)
    all_kspace(:,:,j) = fftshift(fft2(all_img(:,:,j)));
end

kspace = NaN(size(T1map));
kspace_demo = NaN(size(T1map));
kcenter = find(TEcho==TE_eff);
line_start = 1:256/ETL:256;
line_end = line_start+256/ETL-1;
line_ind = kcenter;
echo = kcenter;
for j = 1:16 % going from center to the edges
    lower = kcenter-j;
    if lower<1
        lower = ETL+kcenter-j; % wrap around if reach edge
    end
    upper = kcenter+j;
    if upper>ETL
        upper = kcenter+j-ETL; % wrap around if reach edge
    end
    line_ind = [line_ind lower upper]; % start from center
    echo = [echo upper upper];
end
for j = 1:length(TEcho)
    kspace(line_start(line_ind(j)):line_end(line_ind(j)),:) = all_kspace(line_start(line_ind(j)):line_end(line_ind(j)),:,echo(j));
    kspace_demo(line_start(line_ind(j)):line_end(line_ind(j)),:) = line_ind(j);
end
img_eff = ifft2(kspace);
figure
subplot(1,2,1)
imshow(abs(img_eff))
title('TE_e_f_f = 80ms')
subplot(1,2,2)
contourf(X,Y,kspace_demo,ETL); axis equal; axis tight; axis off; colorbar;
title('kspace')

% 2biii TEeff=40 and 120ms
TE_eff = 40; % ms
kcenter = find(TEcho==TE_eff);
echo = kcenter;
kspace = NaN(size(T1map));
kspace_demo = NaN(size(T1map));
for j = 1:16
    upper = kcenter+j;
    if upper>ETL
        upper = kcenter+j-ETL; % wrap around if reach edge
    end
    echo = [echo upper upper];
end
for j = 1:length(TEcho)
    kspace(line_start(line_ind(j)):line_end(line_ind(j)),:) = all_kspace(line_start(line_ind(j)):line_end(line_ind(j)),:,echo(j));
    kspace_demo(line_start(line_ind(j)):line_end(line_ind(j)),:) = echo(j);
end
img_eff = ifft2(kspace);
figure
subplot(2,2,1)
imshow(abs(img_eff))
title('TE_e_f_f = 40ms')
subplot(2,2,2)
contourf(X,Y,kspace_demo,ETL); axis equal; axis tight; axis off; colorbar;
title('kspace')

TE_eff = 120; % ms
kcenter = find(TEcho==TE_eff);
echo = kcenter;
kspace = NaN(size(T1map));
kspace_demo = NaN(size(T1map));
for j = 1:16
    upper = kcenter+j;
    if upper>ETL
        upper = kcenter+j-ETL; % wrap around if reach edge
    end
    echo = [echo upper upper];
end
for j = 1:length(TEcho)
    kspace(line_start(line_ind(j)):line_end(line_ind(j)),:) = all_kspace(line_start(line_ind(j)):line_end(line_ind(j)),:,echo(j));
    kspace_demo(line_start(line_ind(j)):line_end(line_ind(j)),:) = echo(j);
end
img_eff = ifft2(kspace);
subplot(2,2,3)
imshow(abs(img_eff))
title('TE_e_f_f = 120ms')
subplot(2,2,4)
contourf(X,Y,kspace_demo,ETL); axis equal; axis tight; axis off; colorbar;
title('kspace')

% 2biv TEeff=80ms and ETL=16, 32, 64, 128
tic
ETL_iv = [16 32 64 128];
figure
for l = 1:length(ETL_iv)
    TE_eff = 80; % ms
    ETL = ETL_iv(l);
    TEcho = linspace(ESP,ESP*ETL,ETL); % Varying TE
    all_kspace = NaN([size(T1map),length(TEcho)]); % preallocate kspace for all TE
    all_img = NaN([size(T1map),length(TEcho)]); % preallocate images for all TE
    for j = 1:width(T1map)
        for k = 1:height(T1map)
            Msignal = fsesignal(T1map(k,j),T2map(k,j),ESP,TR,0,ETL,1);
            Msignal(isnan(Msignal)) = 0;
            all_img(k,j,:) = Msignal;
        end
    end
    all_img = repmat(M0map,[1 1 length(TEcho)]).*all_img;
    for j = 1:length(TEcho)
        all_kspace(:,:,j) = fftshift(fft2(all_img(:,:,j)));
    end
    
    kspace = NaN(size(T1map));
    kspace_demo = NaN(size(T1map));
    kcenter = find(TEcho==TE_eff);
    line_start = 1:256/ETL:256;
    line_end = line_start+256/ETL-1;
    line_ind = kcenter;
    echo = kcenter;
    for j = 1:ETL/2
        lower = kcenter-j;
        if lower<1
            lower = ETL+kcenter-j; % wrap around if reach edge
        end
        upper = kcenter+j;
        if upper>ETL
            upper = kcenter+j-ETL; % wrap around if reach edge
        end
        line_ind = [line_ind lower upper];
        echo = [echo upper upper];
    end
    for j = 1:length(TEcho)
        kspace(line_start(line_ind(j)):line_end(line_ind(j)),:) = all_kspace(line_start(line_ind(j)):line_end(line_ind(j)),:,echo(j));
        kspace_demo(line_start(line_ind(j)):line_end(line_ind(j)),:) = echo(j);
    end
    img_eff = ifft2(kspace);
    
    subplot(2,2,l)
    imshow(abs(img_eff))
    title(['TE_e_f_f = 80ms,ETL=',num2str(ETL)])
    toc
end