clear;
clc;
close all;
Fn=200;
set_mw=0.1;
T = 1/Fn;

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
zzz=zeros(5,10); 


RSS1_ALL=files('D:\code\data\FFT_DATA\1\*.dat');
RSS2_ALL=files('D:\code\data\FFT_DATA\2\*.dat');
RSS3_ALL=files('D:\code\data\FFT_DATA\3\*.dat');
RSS4_ALL=files('D:\code\data\FFT_DATA\4\*.dat');
RSS5_ALL=files('D:\code\data\FFT_DATA\5\*.dat');
RSS6_ALL=files('D:\code\data\FFT_DATA\6\*.dat');
RSS7_ALL=files('D:\code\data\FFT_DATA\7\*.dat');
RSS8_ALL=files('D:\code\data\FFT_DATA\8\*.dat');
RSS9_ALL=files('D:\code\data\FFT_DATA\9\*.dat');
RSS10_ALL=files('D:\code\data\FFT_DATA\10\*.dat');


% RSS1_ALL=files('D:\code\data\1227\volunteer\zsy\1\*.dat');
% RSS2_ALL=files('D:\code\data\1227\volunteer\zsy\2\*.dat');
% RSS3_ALL=files('D:\code\data\1227\volunteer\zsy\3\*.dat');
% RSS4_ALL=files('D:\code\data\1227\volunteer\zsy\4\*.dat');
% RSS5_ALL=files('D:\code\data\1227\volunteer\zsy\5\*.dat');
% RSS6_ALL=files('D:\code\data\1227\volunteer\zsy\6\*.dat');
% RSS7_ALL=files('D:\code\data\1227\volunteer\zsy\7\*.dat');
% RSS8_ALL=files('D:\code\data\1227\volunteer\zsy\8\*.dat');
% RSS9_ALL=files('D:\code\data\1227\volunteer\zsy\9\*.dat');
% RSS10_ALL=files('D:\code\data\1227\volunteer\zsy\10\*.dat');


%静态
RSS0_denoised=zeros(size(RSS0));
RSS0_mw=zeros(size(RSS0));
RSS0_denoised_mw=zeros(size(RSS0_mw));
%静态
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

%瑞利滤波
r=9;
r_sigma=4;
Rayleightemp=ones(1,r*2-1);
for i=1:r*2-1
    Rayleightemp(i) = (i-1 )/ (r_sigma^2) * exp(-(i-1)^2 / (2 * r_sigma^2));
end
plot(Rayleightemp)
Rayleightemp = Rayleightemp / sum(Rayleightemp);
[maxr,max_position]=max(Rayleightemp);
%静态
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
plot(rssi0_smooth)


% 设定采样频率
Fs = 200; % 采样频率为200Hz
% 对数据做FFT
windowSize = length(rssi0_smooth); %窗口大小
RSS0length = length(rssi0_smooth); %静态数据总长度
N = length(rssi0_smooth);
Y = fft(rssi0_smooth); % 对数据进行FFT
% 计算频率轴
f = (0:N-1)*(Fs/N); % 频率分布
% 只取一半频率（因为频谱是对称的）
P2 = abs(Y/N); % 归一化幅度
P1 = P2(1:N/2+1); % 只取前半部分
P1(2:end-1) = 2*P1(2:end-1); % 处理双边谱
f1 = f(1:N/2+1); % 相应的频率范围
f1 = f1(2:end); % 移除第一个频率点（0 Hz）
P1 = P1(2:end); % 相应地移除第一个幅度值
% 计算能量
energy_total = sum(P1); % 总能量

indices = f1 < 10; 
energy_below_10Hz = sum(P1(indices)); 

threshold= (energy_below_10Hz / energy_total);



%动态
RSSall={RSS1_ALL,RSS2_ALL,RSS3_ALL,RSS4_ALL,RSS5_ALL,RSS6_ALL,RSS7_ALL,RSS8_ALL,RSS9_ALL,RSS10_ALL};
accuracy_all_f=zeros(1,size(RSSall,2));
accuracy_all_f_n=zeros(1,size(RSSall,2));
for zz=1:size(RSSall,2)
    RSS_all=RSSall{zz};
    RSS_denoised_all=cell(size(RSS_all));
    RSS_denoised_mw=cell(size(RSS_all));
    RSS_denoised_all_mw=cell(size(RSS_denoised_mw));
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
    %动态
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
    sum_ratio=zeros(1,size(rssi_smooth,2));
    for ii=1:size(rssi_smooth,2)
        RSS_intrusion =rssi_smooth{ii};
        RSS_intrusion_length = length(RSS_intrusion); 
        energy_below_10Hz_intrusion = zeros(1,RSS_intrusion_length);
       
            Y = fft(RSS_intrusion); 
            f = (0:RSS_intrusion_length-1)*(Fs/RSS_intrusion_length); 
            P2 = abs(Y/RSS_intrusion_length);
            P1 = P2(1:RSS_intrusion_length/2+1); 
            P1(2:end-1) = 2*P1(2:end-1); 
            f1 = f(1:RSS_intrusion_length/2+1); 
            f1 = f1(2:end); 
            P1 = P1(2:end); 
            energy_total = sum(P1);
            indices = f1 < 10; 
            energy_below_10Hz = sum(P1(indices)); 
            percentage_below_10Hz_intrusion= (energy_below_10Hz / energy_total);
        
        sum_ratio(:,ii)=percentage_below_10Hz_intrusion;
    end
    zzz(:,zz)=sum_ratio';

end
test=zzz;
boxplot(test(:,1:end-1));
therods=ones(1,9)*threshold;
hold on
plot(therods,'r--')
xlabel('Distance(m)');
ylabel('Power sum ratio');
set(findobj(gca,'type','line'),'linewidth',1.5)
set(gcf, 'Color', 'white');
set(gca, 'Color', 'white');








