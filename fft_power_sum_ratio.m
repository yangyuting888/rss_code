clear;
clc;
close all;
Fn=200;
set_mw=0.1;
T=8;
num_samples = Fn * T;
time_axis = linspace(0, T, num_samples);
[~,~,~,rssis]=load_data("D:\code\data\1201\static.dat");
csi_trace=read_bf_file("D:\code\data\1201\static.dat");
agc=zeros(1,size(csi_trace,2));
for ii=1:size(agc,2)
    agc(:,ii)=csi_trace{ii}.agc;
end
rssi_mag=0;
rssi_mag = rssi_mag + dbinv(rssis(:,1).');
rss = db(rssi_mag, 'pow') - 44 - agc;
RSS0=rss(Fn+1:9*Fn);
RSS0_denoised=zeros(size(RSS0));
RSS0_mw=zeros(size(RSS0));
RSS0_denoised_mw=zeros(size(RSS0));

for ii=1:size(RSS0_mw,2)
    RSS0_mw(:,ii)= 10^(RSS0(:,ii)/ 10);
end
for ii=1:size(RSS0_denoised_mw,2)
    RSS0_denoised_mw(:,ii)=(RSS0_mw(:,ii)-min(RSS0_mw(RSS0_mw~=0)))/min(RSS0_mw(RSS0_mw~=0)); 
end
RSS0_denoised_mw(RSS0_denoised_mw==0)=min(RSS0_denoised_mw(RSS0_denoised_mw~=0))-set_mw;
for ii=1:size(RSS0_denoised,2)
    RSS0_denoised(:,ii)=10*log10(RSS0_denoised_mw(:,ii));
end

r=9;
r_sigma=4;
Rayleightemp=ones(1,r*2-1);
for i=1:r*2-1
    Rayleightemp(i) = (i-1 )/ (r_sigma^2) * exp(-(i-1)^2 / (2 * r_sigma^2));
end
Rayleightemp = Rayleightemp / sum(Rayleightemp);
[maxr,max_position]=max(Rayleightemp);
rssi0_smooth=zeros(size(RSS0_denoised));
for ii=1:size(rssi0_smooth,2)
    if ii<max_position
       rssi0_smooth(:,ii) = [zeros(1, max_position-ii),RSS0_denoised(:,1:ii+2*r-1-max_position)]*Rayleightemp';
    elseif ii+2*r-1-max_position>size(rssi0_smooth,2)
        rssi0_smooth(:,ii) = [RSS0_denoised(:,ii-max_position+1:size(rssi0_smooth,2)),zeros(1,ii+2*r-1-max_position-size(rssi0_smooth,2))]*Rayleightemp';
    else
        rssi0_smooth(:,ii) = RSS0_denoised(:,ii-max_position+1 : ii+2*r-1-max_position)*Rayleightemp';
    end
end


%高斯滤波
% sigma=4;
% r=9;
% pi=3.14;
% GaussTemp = ones(1,r*2-1);
% for i=1 : r*2-1
%     GaussTemp(i) = exp(-(i-r)^2/(2*sigma^2))/(sigma*sqrt(2*pi));
% end
% rssi0_smooth=zeros(size(RSS0_denoised));
% for ii=1:size(rssi0_smooth,2)
%         if ii<r
%             rssi0_smooth(:,ii) = [zeros(1, r-ii),RSS0_denoised(:,1:ii+r-1)]*GaussTemp';
%         elseif ii+r-1>size(rssi0_smooth,2)
%             rssi0_smooth(:,ii)=[RSS0_denoised(:,ii-r+1:size(rssi0_smooth,2)),zeros(1, ii+r-1-size(rssi0_smooth,2))]*GaussTemp';
%         else
%             rssi0_smooth(:,ii)=RSS0_denoised(:,ii-r+1 : ii+r-1)*GaussTemp';
%         end
%     
% end

Fs = 200; 
windowSize = 400; 
RSS0length = length(rssi0_smooth); 
N = windowSize;
energy_below_30Hz = zeros(1,RSS0length-windowSize);
for i = 1:RSS0length-windowSize
    Y = fft(rssi0_smooth(i:i+windowSize-1)); 
    f = (0:N-1)*(Fs/N); 
    P2 = abs(Y/N); 
    P1 = P2(1:N/2+1); 
    P1(2:end-1) = 2*P1(2:end-1); 
    f1 = f(1:N/2+1); 
    f1 = f1(2:end); 
    P1 = P1(2:end);
    energy_total = sum(P1); 
    indices = f1 < 10; 
    energy_below_10Hz = sum(P1(indices));
    percentage_below_10Hz(1,i) = (energy_below_10Hz / energy_total);
end
[threshold, max_idx] = max(percentage_below_10Hz);

%tx-rx:1(corridor)
RSS1=get_rss("D:\code\rssi\rss_data\corridor\txrx1\1.dat");
RSS2=get_rss("D:\code\rssi\rss_data\corridor\txrx1\2.dat");
RSS3=get_rss("D:\code\rssi\rss_data\corridor\txrx1\3.dat");
RSS4=get_rss("D:\code\rssi\rss_data\corridor\txrx1\4.dat");
RSS5=get_rss("D:\code\rssi\rss_data\corridor\txrx1\5.dat");
RSS6=get_rss("D:\code\rssi\rss_data\corridor\txrx1\6.dat");
RSS7=get_rss("D:\code\rssi\rss_data\corridor\txrx1\7.dat");
RSS8=get_rss("D:\code\rssi\rss_data\corridor\txrx1\8.dat");
RSS9=get_rss("D:\code\rssi\rss_data\corridor\txrx1\9.dat");
RSS10=get_rss("D:\code\rssi\rss_data\corridor\txrx1\10.dat");
%tx-rx:2(corridor)
% RSS1=get_rss("D:\code\rssi\rss_data\corridor\txrx2\1.dat");
% RSS2=get_rss("D:\code\rssi\rss_data\corridor\txrx2\2.dat");
% RSS3=get_rss("D:\code\rssi\rss_data\corridor\txrx2\3.dat");
% RSS4=get_rss("D:\code\rssi\rss_data\corridor\txrx2\4.dat");
% RSS5=get_rss("D:\code\rssi\rss_data\corridor\txrx2\5.dat");
% RSS6=get_rss("D:\code\rssi\rss_data\corridor\txrx2\6.dat");
% RSS7=get_rss("D:\code\rssi\rss_data\corridor\txrx2\7.dat");
% RSS8=get_rss("D:\code\rssi\rss_data\corridor\txrx2\8.dat");
% RSS9=get_rss("D:\code\rssi\rss_data\corridor\txrx2\9.dat");
% RSS10=get_rss("D:\code\rssi\rss_data\corridor\txrx2\10.dat");
%tx-rx:3(corridor)
% RSS1=get_rss("D:\code\rssi\rss_data\corridor\txrx3\1.dat");
% RSS2=get_rss("D:\code\rssi\rss_data\corridor\txrx3\2.dat");
% RSS3=get_rss("D:\code\rssi\rss_data\corridor\txrx3\3.dat");
% RSS4=get_rss("D:\code\rssi\rss_data\corridor\txrx3\4.dat");
% RSS5=get_rss("D:\code\rssi\rss_data\corridor\txrx3\5.dat");
% RSS6=get_rss("D:\code\rssi\rss_data\corridor\txrx3\6.dat");
% RSS7=get_rss("D:\code\rssi\rss_data\corridor\txrx3\7.dat");
% RSS8=get_rss("D:\code\rssi\rss_data\corridor\txrx3\8.dat");
% RSS9=get_rss("D:\code\rssi\rss_data\corridor\txrx2\9.dat");
% RSS10=get_rss("D:\code\rssi\rss_data\corridor\txrx3\10.dat");

%tx-rx:1(meeting)
% RSS1=get_rss("D:\code\rssi\rss_data\meeting\txrx1\1.dat");
% RSS2=get_rss("D:\code\rssi\rss_data\meeting\txrx1\2.dat");
% RSS3=get_rss("D:\code\rssi\rss_data\meeting\txrx1\3.dat");
% RSS4=get_rss("D:\code\rssi\rss_data\meeting\txrx1\4.dat");
% RSS5=get_rss("D:\code\rssi\rss_data\meeting\txrx1\5.dat");
% RSS6=get_rss("D:\code\rssi\rss_data\meeting\txrx1\6.dat");
% RSS7=get_rss("D:\code\rssi\rss_data\meeting\txrx1\7.dat");
% RSS8=get_rss("D:\code\rssi\rss_data\meeting\txrx1\8.dat");
%tx-rx:2(meeting)
% RSS1=get_rss("D:\code\rssi\rss_data\meeting\txrx2\1.dat");
% RSS2=get_rss("D:\code\rssi\rss_data\meeting\txrx2\2.dat");
% RSS3=get_rss("D:\code\rssi\rss_data\meeting\txrx2\3.dat");
% RSS4=get_rss("D:\code\rssi\rss_data\meeting\txrx2\4.dat");
% RSS5=get_rss("D:\code\rssi\rss_data\meeting\txrx2\5.dat");
% RSS6=get_rss("D:\code\rssi\rss_data\meeting\txrx2\6.dat");
% RSS7=get_rss("D:\code\rssi\rss_data\meeting\txrx2\7.dat");
% RSS8=get_rss("D:\code\rssi\rss_data\meeting\txrx2\8.dat");
%tx-rx:3(meeting)
% RSS1=get_rss("D:\code\rssi\rss_data\meeting\txrx3\1.dat");
% RSS2=get_rss("D:\code\rssi\rss_data\meeting\txrx3\2.dat");
% RSS3=get_rss("D:\code\rssi\rss_data\meeting\txrx3\3.dat");
% RSS4=get_rss("D:\code\rssi\rss_data\meeting\txrx3\4.dat");
% RSS5=get_rss("D:\code\rssi\rss_data\meeting\txrx3\5.dat");
% RSS6=get_rss("D:\code\rssi\rss_data\meeting\txrx3\6.dat");
% RSS7=get_rss("D:\code\rssi\rss_data\meeting\txrx3\7.dat");
% RSS8=get_rss("D:\code\rssi\rss_data\meeting\txrx3\8.dat");




%未处理数据
%RSS_all={RSS1,RSS2,RSS3,RSS4,RSS5,RSS6,RSS7,RSS8,RSS9,RSS10};%1m-10m(cor)
RSS_all={RSS1,RSS2,RSS3,RSS4,RSS5,RSS6,RSS7,RSS8};%1m-8m(meeting)
RSS_denoised_mw=cell(size(RSS_all));
RSS_denoised_all=cell(size(RSS_all));
RSS_denoised_all_mw=cell(size(RSS_all));


for ii=1:size(RSS_denoised_mw,2)
    for jj=1:size(RSS_all{ii},2)
        RSS_denoised_mw{ii}(:,jj)=10^(RSS_all{ii}(:,jj)/ 10);
    end
end

for ii=1:size(RSS_denoised_all,2)
    for jj=1:size(RSS_denoised_mw{ii},2)
        RSS_denoised_all_mw{ii}(:,jj)=(RSS_denoised_mw{ii}(:,jj)-min(RSS_denoised_mw{ii}(RSS_denoised_mw{ii}~=0)))/min(RSS_denoised_mw{ii}(RSS_denoised_mw{ii}~=0));
    end
end
for ii=1:size(RSS_denoised_all,2)
    RSS_denoised_all_mw{ii}(RSS_denoised_all_mw{ii}==0)=min(RSS_denoised_all_mw{ii}(RSS_denoised_all_mw{ii}~=0))-set_mw;
end
for ii=1:size(RSS_denoised_all,2)
    for jj=1:size(RSS_denoised_all_mw{ii},2)
        RSS_denoised_all{ii}(:,jj)=10*log10(RSS_denoised_all_mw{ii}(:,jj));
    end
end
r=9;
r_sigma=4;
Rayleightemp=ones(1,r*2-1);
for i=1:r*2-1
    Rayleightemp(i) = (i-1 )/ (r_sigma^2) * exp(-(i-1)^2 / (2 * r_sigma^2));
end
Rayleightemp = Rayleightemp / sum(Rayleightemp);
[maxr,max_position]=max(Rayleightemp);

rssi_smooth=cell(size(RSS_denoised_all));
for ii=1:size(rssi_smooth,2)
    for jj=1:size(RSS_denoised_all{ii},2)
        if jj<max_position
            rssi_smooth{ii}(:,jj) = [zeros(1, max_position-jj),RSS_denoised_all{ii}(:,1:jj+2*r-1-max_position)]*Rayleightemp';
        elseif jj+2*r-1-max_position>size(RSS_denoised_all{ii},2)
            rssi_smooth{ii}(:,jj)=[RSS_denoised_all{ii}(:,jj-max_position+1:size(RSS_denoised_all{ii},2)),zeros(1, jj+2*r-1-max_position-size(RSS_denoised_all{ii},2))]*Rayleightemp';
        else
            rssi_smooth{ii}(:,jj)=RSS_denoised_all{ii}(:,jj-max_position+1 : jj+2*r-1-max_position)*Rayleightemp';
        end
        
    end
end


% %高斯滤波
% sigma=4;
% r=9;
% pi=3.14;
% GaussTemp = ones(1,r*2-1);
% for i=1 : r*2-1
%     GaussTemp(i) = exp(-(i-r)^2/(2*sigma^2))/(sigma*sqrt(2*pi));
% end
% rssi_smooth=cell(size(RSS_denoised_all));
% for ii=1:size(rssi_smooth,2)
%     for jj=1:size(RSS_denoised_all{ii},2)
%         if jj<r
%             rssi_smooth{ii}(:,jj) = [zeros(1, r-jj),RSS_denoised_all{ii}(:,1:jj+r-1)]*GaussTemp';
%         elseif jj+r-1>size(RSS_denoised_all{ii},2)
%             rssi_smooth{ii}(:,jj)=[RSS_denoised_all{ii}(:,jj-r+1:size(RSS_denoised_all{ii},2)),zeros(1, jj+r-1-size(RSS_denoised_all{ii},2))]*GaussTemp';
%         else
%             rssi_smooth{ii}(:,jj)=RSS_denoised_all{ii}(:,jj-r+1 : jj+r-1)*GaussTemp';
%         end
%     end
% end
%hold on;
%plot(time_axis,rssi_smooth{9},'b');



RSS_intrusion = rssi_smooth{1};
RSS_intrusion_length = length(RSS_intrusion); 
N = windowSize;
energy_below_20Hz_intrusion = zeros(1,RSS_intrusion_length-windowSize);
RSS_detection=zeros(1,size(rssi_smooth,2));
for ii=1:size(RSS_detection,2)
    for i = 1:RSS_intrusion_length-windowSize
        Y = fft(rssi_smooth{ii}(i:i+windowSize-1));
        f = (0:N-1)*(Fs/N);
        P2 = abs(Y/N);
        P1 = P2(1:N/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f1 = f(1:N/2+1);
        f1 = f1(2:end);
        P1 = P1(2:end);
        energy_total = sum(P1);
        indices = f1 < 10;
        energy_below_10Hz = sum(P1(indices));
        percentage_below_10Hz_intrusion(1,i) = (energy_below_10Hz / energy_total);

    end
    num_greater_than_threshold = sum(percentage_below_10Hz_intrusion > threshold);
    detection_percent = num_greater_than_threshold/length(percentage_below_10Hz_intrusion(~isnan(percentage_below_10Hz_intrusion)));
    RSS_detection(1,ii)=detection_percent*100;

end

figure;
hold on;
plot(RSS_detection,'b-o','LineWidth', 2);
xticklabels({'1','2','3','4','5','6','7','8','9','10'});  
set(gcf, 'Color', 'white');
set(gca, 'Color', 'white');
xlabel('Distance(m)');
ylabel('Detection rate(%)');
box(gca,'on');
