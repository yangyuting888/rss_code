clear;
clc;
close all;
Fn=200;
windowsize=200;
set_mw=0.1;
[~,~,~,rssis]=load_data("D:\code\rssi\rss_data\cwt\CWT_DATA\static.dat");
csi_trace=read_bf_file("D:\code\rssi\rss_data\cwt\CWT_DATA\static.dat");
agc=zeros(1,size(csi_trace,2));
for ii=1:size(agc,2)
    agc(:,ii)=csi_trace{ii}.agc;
end
rssi_mag=0;
rssi_mag = rssi_mag + dbinv(rssis(:,1).');
rss = db(rssi_mag, 'pow') - 44 - agc;
RSS0=rss(Fn+1:9*Fn);
RSS1_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\1\*.dat');
RSS2_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\2\*.dat');
RSS3_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\3\*.dat');
RSS4_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\4\*.dat');
RSS5_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\5\*.dat');
RSS6_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\6\*.dat');
RSS7_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\7\*.dat');
RSS8_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\8\*.dat');
RSS9_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\9\*.dat');
RSS10_ALL=files('D:\code\rssi\rss_data\cwt\CWT_DATA\10\*.dat');



zzz=zeros(5,10);


RSS0_denoised=zeros(size(RSS0));
RSS0_mw=zeros(size(RSS0));
RSS0_denoised_mw=zeros(size(RSS0_mw));

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
plot(Rayleightemp)
Rayleightemp = Rayleightemp / sum(Rayleightemp);%归一化
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


wname = 'db4'; 
level = 5; 

[rssi0_C,rssi0_L] = wavedec(rssi0_smooth,level,wname);
c_thresh_rss0={0,0,0};
rssi0_A3 = appcoef (rssi0_C,rssi0_L,wname,level);
for i = 1:level
    cD_rss0 = detcoef(rssi0_C,rssi0_L,i);
    thr_rss0 = thselect(cD_rss0,'rigrsure');
    cD_rss0 = wthresh(cD_rss0,'s',thr_rss0);
    c_thresh_rss0{i}=cD_rss0;
end
cl_rss0=[rssi0_A3 c_thresh_rss0{level} c_thresh_rss0{level-1} c_thresh_rss0{level-2} c_thresh_rss0{level-3} c_thresh_rss0{level-4} ];
rssi0_dwt = waverec(cl_rss0,rssi0_L,wname);


wavename='db4';
totalscal=100;
wcf=centfrq(wavename);
cparam=2*wcf*totalscal;
scal=cparam./(1:totalscal);
freq=scal2frq(scal,wavename,1/Fn);


coefs_rss0_n=cwt(RSS0,scal,wavename);
sp_rss0_n = max(abs(coefs_rss0_n),[],2);
[sp_max_rss0_n,sp_position_rss0_n]=max(sp_rss0_n);
freq_interest_rss0_n=freq(sp_position_rss0_n);
N_rss0_n = ceil(log2((Fn/2)/ freq_interest_rss0_n));
[rss0_c_n,rss0_l_n] = wavedec(RSS0,N_rss0_n,wname);
detail_coefficient_rss0_n = detcoef(rss0_c_n,rss0_l_n,N_rss0_n);
var_detail0_n=var(detail_coefficient_rss0_n);


coefs_rss0=cwt(rssi0_dwt,scal,wavename);
sp_rss0 = max(abs(coefs_rss0),[],2);
[sp_max_rss0,sp_position_rss0]=max(sp_rss0);
freq_interest_rss0=freq(sp_position_rss0);
N_rss0 = ceil(log2((Fn/2)/ freq_interest_rss0));
[rss0_c,rss0_l] = wavedec(rssi0_dwt,N_rss0,wname);
detail_coefficient_rss0 = detcoef(rss0_c,rss0_l,N_rss0);
var_detail0=var(detail_coefficient_rss0);



RSSall={RSS1_ALL,RSS2_ALL,RSS3_ALL,RSS4_ALL,RSS5_ALL,RSS6_ALL,RSS7_ALL,RSS8_ALL,RSS9_ALL,RSS10_ALL};
accuracy_all_c=zeros(1,size(RSSall,2));
accuracy_all_c_n=zeros(1,size(RSSall,2));
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
    rssi_dwt=cell(size(rssi_smooth));
    for ii=1:size(rssi_dwt,2)
        [C,L] = wavedec(rssi_smooth{ii},level,wname); % 小波分解
        c_thresh_rss={0,0,0};
        A3 = appcoef (C,L,wname,level);
        for jj=1:level
            cD_rss = detcoef(C,L,jj);
            thr_rss = thselect(cD_rss,'rigrsure');
            cD_rss = wthresh(cD_rss,'s',thr_rss);
            c_thresh_rss{jj}=cD_rss;
        end
        cl_rss=[A3 c_thresh_rss{level} c_thresh_rss{level-1} c_thresh_rss{level-2} c_thresh_rss{level-3} c_thresh_rss{level-4} ];
        rssi_dwt{ii}=waverec(cl_rss,L,wname);
    end


    var_detail_n=zeros(1,size(RSS_all,2));
    for ii=1:size(var_detail_n,2)
        coefs_rss=cwt(RSS_all{ii},scal,wavename);%得到小波系数
        sp_rss = max(abs(coefs_rss),[],2); % 沿着时间方向（第二维）求每一行（每一尺度）的最大值，得到尺度投影向量sp
        [sp_max_rss,sp_position_rss]=max(sp_rss);
        freq_interest_rss=freq(sp_position_rss);
        N_rss = ceil(log2((Fn/2)/ freq_interest_rss));
        [rss_c,rss_l] = wavedec(RSS_all{ii},N_rss,wname);
        detail_coefficient_rss = detcoef(rss_c,rss_l,N_rss);
        var_detail_n(:,ii)=var(detail_coefficient_rss);
    end
    var_detail=zeros(1,size(rssi_dwt,2));
    for ii=1:size(var_detail,2)
        coefs_rss=cwt(rssi_dwt{ii},scal,wavename);
        sp_rss = max(abs(coefs_rss),[],2); 
        [sp_max_rss,sp_position_rss]=max(sp_rss);
        freq_interest_rss=freq(sp_position_rss);
        N_rss = ceil(log2((Fn/2)/ freq_interest_rss));
        [rss_c,rss_l] = wavedec(rssi_dwt{ii},N_rss,wname);
        detail_coefficient_rss = detcoef(rss_c,rss_l,N_rss);
        var_detail(:,ii)=var(detail_coefficient_rss);
    end

   %zzz(:,zz)=var_detail';
   zzz(:,zz)=var_detail_n';

   

end

boxplot(zzz(:,1:end-1));
hold on
%therods=ones(1,9)*var_detail0;
therods=ones(1,9)*var_detail0_n';
plot(therods,'r--')
xlabel('Distance(m)');
ylabel('Detail coefficient var');
xticklabels({'1','2','3','4','5','6','7','8','9'}); 
set(gcf, 'Color', 'white');
set(gca, 'Color', 'white');





