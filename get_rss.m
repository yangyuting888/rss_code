function RSS=get_rss(file_name)
Fn=200;
[~,~,~,rssis]=load_data(file_name);
csi_trace=read_bf_file(file_name);
agc=zeros(1,size(csi_trace,1));
for ii=1:size(agc,2)
    agc(:,ii)=csi_trace{ii}.agc;
end
 rssi_mag=0;
 %rssi_mag = rssi_mag + dbinv(rssis(:,1).');
 rssi_mag = rssi_mag + dbinv(rssis(:,3).');
 rss = db(rssi_mag, 'pow') - 44 - agc;
 RSS=rss(Fn+1:9*Fn);
end

