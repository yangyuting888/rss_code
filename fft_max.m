clear;
clc;
close all;
Fn=200;
set_mw=0.1;
[~,~,~,rssis]=load_data("D:\code\rssi\rss_data\fft_max\static.dat");
csi_trace=read_bf_file("D:\code\rssi\rss_data\fft_max\static.dat");
agc=zeros(1,size(csi_trace,2));
for ii=1:size(agc,2)
    agc(:,ii)=csi_trace{ii}.agc;
end
rssi_mag=0;
rssi_mag = rssi_mag + dbinv(rssis(:,1).');
rss = db(rssi_mag, 'pow') - 44 - agc;
RSS0=rss(Fn+1:9*Fn);

RSS1=get_rss("D:\code\rssi\rss_data\fft_max\1.dat");
RSS2=get_rss("D:\code\rssi\rss_data\fft_max\2.dat");
RSS3=get_rss("D:\code\rssi\rss_data\fft_max\3.dat");
RSS4=get_rss("D:\code\rssi\rss_data\fft_max\4.dat");
RSS5=get_rss("D:\code\rssi\rss_data\fft_max\5.dat");
RSS6=get_rss("D:\code\rssi\rss_data\fft_max\6.dat");
RSS7=get_rss("D:\code\rssi\rss_data\fft_max\7.dat");
RSS8=get_rss("D:\code\rssi\rss_data\fft_max\8.dat");
RSS9=get_rss("D:\code\rssi\rss_data\fft_max\9.dat");
RSS10=get_rss("D:\code\rssi\rss_data\fft_max\10.dat");
boud=[10,9,8,7,6,5,4,3,2,1];
RSS_all={RSS10,RSS9,RSS8,RSS7,RSS6,RSS5,RSS4,RSS3,RSS2,RSS1};%10m-1m
T = 1/Fn;
RSS_amplitude=zeros(1,size(RSS_all,2));
RSSamplitude=cell(size(RSS_all));
for ii=1:size(RSS_amplitude,2)
    L = length(RSS_all{ii}); 
    t = (0:L-1)*T;
    RSS_fft=fft(RSS_all{ii});
    P2 = abs(RSS_fft/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fn*(0:(L/2))/L;
    P1(1)=0;
    RSSamplitude{ii}=P1;
    RSS_amplitude(:,ii)=max(P1);
end
L0 = length(RSS0); 
t0 = (0:L0-1)*T; 
RSS0_fft=fft(RSS0);
P20 = abs(RSS0_fft/L0);
P10 = P20(1:L0/2+1);
P10(2:end-1) = 2*P10(2:end-1);
f0 = Fn*(0:(L0/2))/L0;
P10(1)=0;
[RSS0_amplitude,index0]=max(P10);
position=find(RSS_amplitude>RSS0_amplitude,1,'first');
dmax=boud(:,position);
fprintf("能够检测到人员入侵的最远距离为：%.2f米",dmax);
therods=ones(size(RSS_amplitude(:,1:end)))*RSS0_amplitude;
figure;
hold on;
plot(flip(RSS_amplitude),'b-o','LineWidth', 2);
plot(therods,'r--','LineWidth', 2);
xticklabels({'1','2','3','4','5','6','7','8','9','10'});  
set(gcf, 'Color', 'white');
set(gca, 'Color', 'white');
xlabel('Distance(m)');
ylabel('FFT amplitude');
box(gca,'on');



























