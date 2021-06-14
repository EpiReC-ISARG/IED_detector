function [out, discharges, d_decim, envelope, background, envelope_pdf]=spike_detector_hilbert_v23(d,fs,settings)

% Mutichannel spike detector using statistical description of bacground
% activity to finding of spike suspected section. Algorithm allows resampling 
% of signal (DEFAULT 200 Hz), detect fast discharges (spike, poly-spike) in defined
% frequency band. Algorithm contains main hum filtering (50 Hz DEFAULT).
% Parallel computing, matobject read is supported. Detailed description is in
% paper: Janca, Radek, et al. "Detection of Interictal Epileptiform Discharges 
% Using Signal Envelope Distribution Modelling: Application to Epileptic and 
% Non-Epileptic Intracranial Recordings." Brain topography 28.1 (2015): 172-183.
% DOI: 10.1007/s10548-014-0379-1
%
% _________________________________________________________________________
% recomended calling:
% [...]=spike_detector_hilbert_v20(d,fs);
% [...]=spike_detector_hilbert_v20('c:/datafile.mat "sig_var"',fs); % "name of signal matrix"
% [...]=spike_detector_hilbert_v20('c:/datafile.mat "sig_var" "fs_var"'); % "name of signal matrix" and "name of fs value" 
% or
% [...]=spike_detector_hilbert_v20(d,fs,settings);
% _________________________________________________________________________
%
%
% inputs: -----------------------------------------------------------------
% d ... signal (time x channel)
% fs ... sampling frequency (Hz)
% settings ... string options (exmple: settings='-bl 10 -bh 60 -k1 3.0')
%   -fl low frequency of filtering ('-bl 10' DEFAULT)
%   -fh high frequency of filtering ('-bh 60' DEFAULT)
%   -ft filter type: 1-Chebyshev II, 2-Butterworth, 3-FIR ('-ft 1' DEFAULT)
%   -k1 threshold value for obvious spike decision ('-k1 3.65' DEFAULT)
%   -k2 defines ambiguous spike treshold. Ambiguous 
%       spike is accepted, when simultaneous obvious detection is in other 
%       channel k1 >= k2 (k1 in DEFAULT)
%   -k3 decrease the threshold value (0 in DEFAULT) k1*(mode+median)-k3*(mean-mode);
%   -w winsize - size of segment in sample around spike for background
%                   definition (5*fs DEFAULT)
%   -n noverlap - overlap of segment in sample or absolute (novelap=0.9 ~ 90%)
%                   (4*fs DEFAULT)
%   -buf buffer ... value  of length subsection of signal for batch
%               analyse in seconds. High value with combination of parallel computing, 
%               many analyzed channels is memory-consuming. Minimaly 10x 
%               winsize time is recomended. ('-buf 300' DEFAULT)
%   -h main_hum_freq ... main hum frequency in Hz ('-h 50' DEFAULT)
%   -b low hrequency boundary of beta activity detection beta-25 Hz ('-b 15' recomended)
%       Inf - beta-activity detector off ('-b Inf' DEFAULT)
%   -bw length of tested segment on beta-activity in second
%   -br autoregresive model order of beta-activity detector ('-br 12' DEAFAULT)
%   -dt time between spikes in events
%   -pt polyspike union time: spike in same channel nearest then time will
%        be united
%   -dec allowes signal resampling (frequency in Hz) before processing ('-dec 200' DEFAULT)
%        '-dec 0' disables decimation, preserves original fs
%   -ti time indexing of spike: '-ti 1' maximum of envelope (inflex point)
%                               '-ti 2' maximum of absolute amplitude
%
%         
%
% outputs:-----------------------------------------------------------------
% out ... structure describing each detected spike. Usefull for export to
%         GDF files using biosig toolbox (biosig.sourceforge.net)
%           out.pos ... position (second)
%           out.dur ... duration (second) - fixed 1/fs 
%           out.chan ... channel
%           out.con ... type (1-obvious 0.5-ambiguous)
%           out.weigth ... statistical significance "CDF"
%           out.pdf ... statistical significance "PDF"
% envelope ... instant hilbert envelope of filtered signal "d"
% bacground ... threshold curves with same size as signal "d"
%               background(:,:,1) ... for high threshold (k1)
%               background(:,:,2) ... for low threshol (k2)
% discharges ... structure of multichannel event describes occurence of
%               spikes suddenly 
%           discharges.MV ... matrix of type (n x ch) 
%                             (1-obvious spike, 0.5-ambiguous)
%           discharges.MA ... matrix of max. amplitude of envelope above 
%                             backround (n x ch)
%           discharges.MP ... start position of multichannel event 
%           discharges.MD ... duration of event
%           discharges.MW ... statistical significance "CDF"
%           discharges.MPDF ... probability of occurence
%
% required toolboxes:
%   Signal Processing Toolbox, Statistic Toolbox, 
%
% recomended toolboxes:
%   Parallel Computing Toolbox
%
% MAKING BY ISARG - Radek Janca, 19.9.2013
% REVISION 9.4.2018


% settings defaults:------------------------------------------------------
global bandwidth; bandwidth=[10 60]; % [(-fl) (-fh)]
global k1; k1=3.65; % (-k1)
global k2; k2=k1; % (-k2)
global k3; k3=0; % (-k3)
global buffering; buffering=300; % (-buf)
global main_hum_freq; main_hum_freq=50; % (-h)
global beta; beta=Inf; % (-b)
global beta_win; beta_win=20; % (-bw)
global beta_AR; beta_AR=12; % (-br)
global f_type; f_type=1; % 1-cheby2, 2-but, 3-fir (-ft)
global discharge_tol; discharge_tol=0.005; % (-dt)
global polyspike_union_time; polyspike_union_time=0.12; % (-pt)
global decimation; decimation=200; % (-dec)
global ti_switch; ti_switch=1;

if nargin>2; 
    poz=strfind(settings,'  '); % multi spaces removing
    settings(poz)=[];
    SET=textscan(settings,'%s %f','delimiter',' '); % paremeters reading
else    
    SET{1,1}=[]; 
end

% setting reading
for i=1:size(SET{1,1},1)
    prefix=SET{1,1}{i,:};
    prefix(1)=[];
    
    switch prefix
        case 'fl'
            bandwidth(1)=SET{1,2}(i);
        case 'fh'
            bandwidth(2)=SET{1,2}(i);
        case 'k1'
            k1=SET{1,2}(i);
        case 'k2'
            k2=SET{1,2}(i);
        case 'k3'
            k3=SET{1,2}(i);
        case 'w'
            winsize=SET{1,2}(i);
        case 'n'
            noverlap=SET{1,2}(i);
        case 'buf'
            buffering=SET{1,2}(i);
        case 'h'
            main_hum_freq=SET{1,2}(i);
        case 'b'
            beta=SET{1,2}(i);
        case 'bw'
            beta_win=SET{1,2}(i);
        case 'br'
            beta_AR=SET{1,2}(i);
        case 'ft'
            f_type=SET{1,2}(i);
        case 'dt'
            discharge_tol=SET{1,2}(i);
        case 'pt'
            polyspike_union_time=SET{1,2}(i);
        case 'dec'
            decimation=SET{1,2}(i);
        case 'ti'
            ti_switch=SET{1,2}(i);
    end
end


% matobject read
global useMatObj; useMatObj=false;
if ischar(d)
    useMatObj=true;
    strpoz=strfind(d,'"');
    if isempty(strpoz);
        dataPath=d; dataVarName='d'; fsVarName='fs';
    else
        dataPath=d(1:strpoz(1)-1);
        dataVarName=d(strpoz(1)+1:strpoz(2)-1);
        
        if length(strpoz)>2
            fsVarName=d(strpoz(3)+1:strpoz(4)-1);
        end
    end
    dataObj=matfile(dataPath);
    
    fs=eval(['dataObj.' fsVarName]);
    
    dataSize=size(dataObj,'d');
else
    dataSize=size(d);
end


if ~exist('winsize'); winsize=5*fs; end % (-w)
if ~exist('noverlap');noverlap=4*fs; end % (-n)

if buffering/(winsize/fs)<10; buffering=(winsize/fs)*10; disp(['buffer size increased to ' num2str(buffering) ' sec.']); end
if buffering>(dataSize(1)/fs); buffering=dataSize(1)/fs; disp(['buffer size decreased to ' num2str(buffering) ' sec.']); end

if decimation==0
   decimation=fs; 
end


if bandwidth(2)>decimation
   error('filter -fh frequency is bigger than fs/2') 
end

% signal buffering ---------------------------------------------------
N_seg=floor(dataSize(1)/(buffering*fs));
if N_seg<1; N_seg=1; end
T_seg=round(dataSize(1)/N_seg/fs);

% indexs of segments with two-side overlap
index_start=1:T_seg*fs:dataSize(1);
if length(index_start)>1
    index_start(2:end)=index_start(2:end)-(3*winsize);
    
    index_stop=index_start+T_seg*fs+2*(3*winsize)-1;
    index_stop(1)=index_stop(1)-(3*winsize);
    index_stop(end)=dataSize(1);
    
    if index_stop(end)-index_start(end)<T_seg*fs
        index_start(end)=[];
        index_stop(end-1)=[];
    end
else
    index_stop=dataSize(1);
end



% detection calling ---------------------------------------------------
out.pos=[]; % spike position
out.dur=[]; % spike duration - fix value 5 ms
out.chan=[]; % spike in channel
out.con=[]; % spike condition
out.weight=[]; % spike weight "CDF"
out.pdf=[]; % spike probability "PDF"

discharges.MV=[]; % spike type 1-obvious, 0.5- ambiguous
discharges.MA=[]; % max. amplitude of envelope above backround
discharges.MP=[]; % event start position
discharges.MD=[]; % event duration 
discharges.MW=[]; % statistical weight
discharges.MPDF=[]; % statistical weight
discharges.MRAW=[]; % amplitude of signal

fprintf(1,'progress:   0 %%')

idx_itt=0;
for i=1:length(index_stop)
    
    
    % subsection signal spike detection ----------------
    if useMatObj
        [sub_d,sub_envelope,sub_background,sub_discharges,sub_out,sub_envelope_pdf,RFactor]=spike_detector(eval(['dataObj.' dataVarName '(index_start(i):index_stop(i),:)']),fs,winsize,noverlap);
    else
        [sub_d,sub_envelope,sub_background,sub_discharges,sub_out,sub_envelope_pdf,RFactor]=spike_detector(d(index_start(i):index_stop(i),:),fs,winsize,noverlap);
    end
    
    
    
    % progress stat
    percent=num2str(100*i/length(index_stop),'%3.0f');
    percent=[repmat(' ',1,3-length(percent)) , percent,' %%'];
    fprintf(1,['\b\b\b\b\b' percent])
    
    
    % conecting of subsection ----------------
    wR=round(winsize/RFactor);
    ist=ceil(index_start(i)/RFactor);
    isp=ceil(index_stop(i)/RFactor);
    if nargout>2; d_decim(ist+(i>1)*(3*wR):isp-(i<length(index_stop))*(3*wR),:)=sub_d(1+(i>1)*(3*wR):end-(i<length(index_stop))*(3*wR),:); end
    if nargout>3; envelope(ist+(i>1)*(3*wR):isp-(i<length(index_stop))*(3*wR),:)=sub_envelope(1+(i>1)*(3*wR):end-(i<length(index_stop))*(3*wR),:); end
    if nargout>4; background(ist+(i>1)*(3*wR):isp-(i<length(index_stop))*(3*wR),:,:)=sub_background(1+(i>1)*(3*wR):end-(i<length(index_stop))*(3*wR),:,:); end
    if nargout>5; envelope_pdf(ist+(i>1)*(3*wR):isp-(i<length(index_stop))*(3*wR),:)=sub_envelope_pdf(1+(i>1)*(3*wR):end-(i<length(index_stop))*(3*wR),:); end
    
    % removing of two side overlap detections
    
    if ~isempty(sub_out.pos)
        if length(index_stop)>1
           
            idx_evt=sub_out.pos<((i>1)*(3*winsize)/fs) | sub_out.pos>((index_stop(i)-index_start(i))-(i<length(index_stop))*(3*winsize))/fs;
            idx_disch=min(sub_discharges.MP,[],2)<((i>1)*(3*winsize)/fs) | min(sub_discharges.MP,[],2)>((index_stop(i)-index_start(i))-(i<length(index_stop))*(3*winsize))/fs;
            
            sub_out.pos(idx_evt)=[];
            sub_out.dur(idx_evt)=[];
            sub_out.chan(idx_evt)=[];
            sub_out.con(idx_evt)=[];
            sub_out.weight(idx_evt)=[];
            sub_out.pdf(idx_evt)=[];
            
            sub_discharges.MV(idx_disch,:)=[];
            sub_discharges.MA(idx_disch,:)=[];
            sub_discharges.MRAW(idx_disch,:)=[];

            sub_discharges.MP(idx_disch,:)=[];
            sub_discharges.MD(idx_disch,:)=[];
            sub_discharges.MW(idx_disch,:)=[];
            sub_discharges.MPDF(idx_disch,:)=[];
            
        end
    end
    
    out.pos=[out.pos; sub_out.pos+index_start(i)/fs-1/fs];
    out.dur=[out.dur; sub_out.dur];
    out.chan=[out.chan; sub_out.chan];
    out.con=[out.con; sub_out.con];
    out.weight=[out.weight; sub_out.weight];
    out.pdf=[out.pdf; sub_out.pdf];
    
    discharges.MV=[discharges.MV; sub_discharges.MV];
    discharges.MA=[discharges.MA; sub_discharges.MA];
    discharges.MP=[discharges.MP; sub_discharges.MP+index_start(i)/fs-1/fs];
    discharges.MD=[discharges.MD; sub_discharges.MD];
    discharges.MW=[discharges.MW; sub_discharges.MW];
    discharges.MPDF=[discharges.MPDF; sub_discharges.MPDF];
    discharges.MRAW=[discharges.MRAW; sub_discharges.MRAW];
    
    idx_itt=idx_itt+size(sub_d,1);
end

% removing of detection on start and end of signal (filtering artifact etc.)
out=structfun(@(x) x(out.pos>2 & out.pos<(size(d,1)/fs-2)),out,'UniformOutput',0);
% switch ti_switch
%     case 2
        discharges=structfun(@(x) x(nanmin(discharges.MP,[],2)>2 & nanmin(discharges.MP,[],2)<(size(d,1)/fs-2),:),discharges,'UniformOutput',0);
%     case 1
%         discharges=rmfield(discharges,'MRAW');
%         discharges=structfun(@(x) x(nanmin(discharges.MP,[],2)>2 & nanmin(discharges.MP,[],2)<(size(d,1)/fs-2),:),discharges,'UniformOutput',0);
%         discharges.MRAW=[];
% end

fprintf(1,'\n')
end



function [d_decim,envelope,background,discharges,out,envelope_pdf,RFactor]=spike_detector(d,fs,winsize,noverlap)
global decimation;
global k1; 
global k2;
global k3;
global main_hum_freq;
global beta; 
global beta_win; 
global beta_AR;
global polyspike_union_time;
global discharge_tol;
global ti_switch;

% signal downsampling ------------------------------------------------------

RFactor=fs/decimation;
if RFactor>1 || decimation~=fs
    winsize=winsize/fs;
    noverlap=noverlap/fs;
    
    decItterations=ceil(log10(RFactor));
    
    for i=1:decItterations
        if i==decItterations
            fs_out=decimation;
        else
            fs_out=round((fs/RFactor^(1/decItterations)));
        end
        sub_d_res=zeros(ceil(size(d,1)/(fs/fs_out)), size(d,2));
        
        parfor ch=1:size(d,2)
            sub_d_res(:,ch)=resample(d(:,ch),fs_out,fs,100);
        end
        d=sub_d_res; 
        fs=fs_out;
    end
    
    winsize=round(winsize*fs);
    noverlap=round(noverlap*fs);
end


%--------------------------------------------------------------------------
% segmentation index
%--------------------------------------------------------------------------
if noverlap<1
    index=1:round(winsize*(1-noverlap)):size(d,1)-winsize+1; % indexy zaèátkù oken
else
    index=1:winsize-noverlap:size(d,1)-winsize+1;
end



% 50/60 Hz filterring -----------------------------------------------------

d=filt50Hz(d,fs,main_hum_freq);
[bb,aa]=butter(2,2*1/fs,'high');
d_decim=filtfilt(bb,aa,d); % return decimated signal without main hums

% mu/beta activity filterring (part 1) ------------------------------------
if beta<fs/2 && beta_win>0;M_beta=beta_detect(d,fs,beta,beta_win,beta_AR);end

% badwidth filterring  ----------------------------------------------------
d=filtering(d,fs);


%--------------------------------------------------------------------------
% ONE CHANNEL DETECTING
%--------------------------------------------------------------------------
p1=k1;
p2=k2;
p3=k3;
p4=polyspike_union_time;
p5=ti_switch;

parfor ch=1:size(d,2)
    if sum(d(:,ch)==0)==length(d(:,ch))
        envelope(:,ch) = zeros(size(d(:,ch)));
        markers_high(:,ch) = zeros(size(d(:,ch)));
        markers_low(:,ch) = zeros(size(d(:,ch)));
        if k1==k2
            background(:,ch,:) = zeros(size(d(:,ch),1),1,1);
        else
            background(:,ch,:) = zeros(size(d(:,ch),1),1,2);
        end
        envelope_cdf(:,ch) = zeros(size(d(:,ch)));
        envelope_pdf(:,ch) = zeros(size(d(:,ch)));
    else
        [envelope(:,ch),markers_high(:,ch),markers_low(:,ch),background(:,ch,:),envelope_cdf(:,ch),envelope_pdf(:,ch)]=one_channel_detect(d(:,ch),fs,index,winsize,p1,p2,p3,p4,p5,d_decim(:,ch));
    end
end
    



% first and last second is not analyzed (filter time response etc.) -------
markers_high([1:fs, end-fs+1:end],:)=false;
markers_low([1:fs, end-fs+1:end],:)=false;

% last stage defined by discharge_tol (added in revision 25.10.2017)
markers_high(end-ceil(discharge_tol*fs+1):end,:)=false;
markers_low(end-ceil(discharge_tol*fs+1):end,:)=false;


% mu/beta activity filterring (part 2) ------------------------------------
if beta<fs/2 && beta_win>0
    markers_high(M_beta)=false;
    markers_low(M_beta)=false;
end

ovious_M=sum(markers_high,2)>0;

% obvious spike events output
out.pos=[];
out.dur=[];
out.chan=[];
out.con=[];
out.weight=[];
out.pdf=[];


t_dur=0.005;
for ch=1:size(markers_high,2)
    if sum(markers_high(:,ch))>0
        idx=find(markers_high(:,ch)); 
        
        out.pos=[out.pos; idx(:)/fs];
        out.dur=[out.dur; t_dur*ones(length(idx),1)];
        out.chan=[out.chan; ch*ones(length(idx),1)];
        out.con=[out.con; ones(length(idx),1)];
        out.weight=[out.weight; envelope_cdf(idx(:),ch)];
        out.pdf=[out.pdf; envelope_pdf(idx(:),ch)];
    end
end


if ~(k2==k1)
    % ambiguous spike events output
    for ch=1:size(markers_low,2)
        if sum(markers_low(:,ch))>0
            idx=find(markers_low(:,ch));  
            
            idx(logical(markers_high(idx(:),ch)))=[];
            for i=1:length(idx)
                if sum(ovious_M(round(idx(i)-0.01*fs:idx(i)-0.01*fs)))>0
                    out.pos=[out.pos; idx(i)/fs];
                    out.dur=[out.dur; t_dur];
                    out.chan=[out.chan; ch];
                    out.con=[out.con; 0.5];
                    out.weight=[out.weight; envelope_cdf(idx(i),ch)];
                    out.pdf=[out.pdf; envelope_pdf(idx(i),ch)];
                end
            end
        end
    end
end



%--------------------------------------------------------------------------
% making M stack pointer of events
%--------------------------------------------------------------------------
M=zeros(size(d));
for k=1:size(out.pos,1)
    M(round(out.pos(k)*fs:out.pos(k)*fs+discharge_tol*fs),out.chan(k))=out.con(k);
end

%--------------------------------------------------------------------------
% definition of multichannel events vectors
%--------------------------------------------------------------------------

point(:,1)=find(diff([0;(sum(M,2))>0])>0);
point(:,2)=find(diff([(sum(M,2))>0;0])<0);

discharges.MV=[]; % spike type 1-obvious, 0.5- ambiguous
discharges.MA=[]; % max. amplitude of envelope above backround
discharges.MP=[]; % event start position
discharges.MD=[]; % event duration 
discharges.MW=[]; % statistical weight
discharges.MPDF=[]; % statistical weight
discharges.MRAW=[]; % max. amplitude of signal

for k=1:size(point,1)
    seg=M(point(k,1):point(k,2),:);
    mv=max(seg,[],1);
    
    

    switch ti_switch
        case 1
            seg=envelope(point(k,1):point(k,2),:)-(background(point(k,1):point(k,2),:,1)/k1);
            ma=max(abs(seg),[],1);
            
            seg=d_decim(point(k,1):point(k,2),:);
            [mraw,poz]=max(abs(seg),[],1);
            poz=poz+(0:size(seg,2)-1)*size(seg,1);
            mraw=mraw.*sign(seg(poz));
%             mraw=[];
        case 2
            seg=envelope(point(k,1)-round(fs*10e-3):point(k,2)+round(fs*10e-3),:)-(background(point(k,1)-round(fs*10e-3):point(k,2)+round(fs*10e-3),:,1)/k1);
            ma=max(abs(seg),[],1);
%             
            seg=d_decim(point(k,1):point(k,2),:);
            [mraw,poz]=max(abs(seg),[],1);
            poz=poz+(0:size(seg,2)-1)*size(seg,1);
            mraw=mraw.*sign(seg(poz));
    end
        
    seg=envelope_cdf(point(k,1):point(k,2),:);
    mw=max(seg,[],1);
    
    seg=envelope_pdf(point(k,1):point(k,2),:);
    mpdf=max(seg.*(M(point(k,1):point(k,2),:)>0),[],1);
    
     % ------ precise position in discharges.MP -------
    [row,col]=find(M(point(k,1):point(k,2),:)>0);
    mp=nan(1,size(d,2));
    mp(col)=row+point(k,1)-1;
    
    
    m1s=size(discharges.MRAW,1);
    m2s=size(mv,1);
    if isempty(mv); continue; end
    if isempty(discharges.MV)
        discharges.MRAW=mraw;
        discharges.MV=mv;
        discharges.MA=ma;
        discharges.MW=mw;
        discharges.MPDF=mpdf;
        discharges.MD=repmat(point(k,2)-point(k,1),1,size(d,2))/fs;
        discharges.MP=mp;
    end
    
    discharges.MRAW(m1s+1:m1s+m2s,:,:)=mraw;
    discharges.MV(m1s+1:m1s+m2s,:,:)=mv;
    discharges.MA(m1s+1:m1s+m2s,:,:)=ma;
    discharges.MW(m1s+1:m1s+m2s,:,:)=mw;
    discharges.MPDF(m1s+1:m1s+m2s,:,:)=mpdf;
    discharges.MD(m1s+1:m1s+m2s,:,:)=repmat(point(k,2)-point(k,1),1,size(d,2))/fs;
    discharges.MP(m1s+1:m1s+m2s,:,:)=mp/fs;
   
end

end


function [envelope,markers_high,markers_low,prah_int,envelope_cdf,envelope_pdf]=one_channel_detect(d,fs,index,winsize,k1,k2,k3,polyspike_union_time,ti_switch,d_decim)


envelope=abs(hilbert(d)); % Hilbert's envelope (intense envelope)

for k=1:length(index) % for each segment
    
    segment=envelope(index(k):index(k)+winsize-1);
    segment(segment<=0)=[];
    
    % estimation of segment's distribution using MLE 
    phat(k,1)=mean(log(segment)); 
    phat(k,2)=std(log(segment)); 
end

r=size(envelope,1)/length(index);
n_average=winsize/fs;


if round(n_average*fs/r)>1
    phat=filtfilt(ones(round(n_average*fs/r),1)/(round(n_average*fs/r)),1,phat);
end

% interpolation of thresholds value to threshold curve (like backround)
phat_int=[];
if size(phat,1)>1
    phat_int(:,1) = interp1(index+round(winsize/2),phat(:,1),(index(1):index(end))+round(winsize/2),'spline');
    phat_int(:,2) = interp1(index+round(winsize/2),phat(:,2),(index(1):index(end))+round(winsize/2),'spline');
    
    phat_int=[ones(floor(winsize/2),size(phat,2)).*repmat(phat_int(1,:),floor(winsize/2),1); phat_int ; ones(size(envelope,1)-(length(phat_int)+floor(winsize/2)),size(phat,2)).*repmat(phat_int(end,:),size(envelope,1)-(length(phat_int)+floor(winsize/2)),1)];   
else
    phat_int=phat.*ones(size(d,1),2);
end

lognormal_mode= exp(phat_int(:,1)-phat_int(:,2).^2);
lognormal_median=exp(phat_int(:,1));
lognormal_mean=exp(phat_int(:,1)+(phat_int(:,2).^2)/2);

prah_int(:,1)=k1*(lognormal_mode+lognormal_median)-k3*(lognormal_mean-lognormal_mode);
if ~(k2==k1)
    prah_int(:,2)=k2*(lognormal_mode+lognormal_median)-k3*(lognormal_mean-lognormal_mode);
end



envelope_cdf=0.5+0.5*erf((log(envelope)-phat_int(:,1))./sqrt(2*phat_int(:,2).^2)); % CDF of lognormal distribution
envelope_pdf=exp(-0.5 .* ((log(envelope) - phat_int(:,1))./phat_int(:,2)).^2) ./ (envelope  .* phat_int(:,2) .* sqrt(2*pi)); % PDF of lognormal distribution

% -------------------------------------------------------------------------
% detection of obvious and ambiguous spike
% -------------------------------------------------------------------------
markers_high=local_maxima_detection(envelope,prah_int(:,1),fs,polyspike_union_time,ti_switch,d_decim);
markers_high=detection_union(markers_high,envelope,polyspike_union_time*fs);

if ~(k2==k1)
    markers_low=local_maxima_detection(envelope,prah_int(:,2),fs,polyspike_union_time,ti_switch,d_decim);
    markers_low=detection_union(markers_low,envelope,polyspike_union_time*fs);
else
    markers_low=markers_high;
end

end






function marker1=local_maxima_detection(envelope,prah_int,fs,polyspike_union_time,ti_switch,d_decim)

marker1=zeros(size(envelope));
marker1(envelope(:)>prah_int(:))=1; % crossing of high threshold

point=[];
point(:,1)=find(diff([0;marker1])>0); % strat crossing
point(:,2)=find(diff([marker1;0])<0); % end crossing


switch ti_switch
    case 2
        envelope=abs(d_decim);
end

marker1=false(size(envelope));
for k=1:size(point,1)
    
    % detection of local maxima in section which crossed threshold curve
    if point(k,2)-point(k,1)>2
        seg=envelope(point(k,1):point(k,2));
        seg_s=diff(seg);
        seg_s=sign(seg_s); 
        seg_s=find(diff([0;seg_s])<0); % positions of local maxima in the section
        
        marker1(point(k,1)+seg_s-1)=true;
    elseif point(k,2)-point(k,1)<=2
        seg=envelope(point(k,1):point(k,2));
        [~,s_max]=max(seg); % positions of local maxima in the section
        marker1(point(k,1)+s_max-1)=true;
    end
end


% union of section, where local maxima are close together <(1/f_low + 0.02 sec.)~ 120 ms
pointer=find(marker1==true); % index of local maxima
state_previous=false;
for k=1:length(pointer)
    if ceil(pointer(k)+polyspike_union_time*fs)>size(marker1,1)
        seg=marker1(pointer(k)+1:end);
    else
        seg=marker1(pointer(k)+1:ceil(pointer(k)+polyspike_union_time*fs));
    end
    
    if state_previous
        if sum(seg)>0
            state_previous=true;
        else
            state_previous=false;
            marker1(start:pointer(k))=true;
        end
        
    else
        
        if sum(seg)>0
            state_previous=true;
            start=pointer(k);
        end
    end
end

% finding of the highes maxima of the section with local maxima
point=[];
point(:,1)=find(diff([0;marker1])>0); % start 
point(:,2)=find(diff([marker1;0])<0); % end 

% local maxima with gradient in souroundings
for k=1:size(point,1)
    if point(k,2)-point(k,1)>1
        lokal_max=pointer(pointer>=point(k,1) & pointer<=point(k,2)); % index of local maxima
        
        marker1(point(k,1):point(k,2))=false;
        
        lokal_max_val=envelope(lokal_max); % envelope magnitude in local maxima
        lokal_max_poz=(diff(sign(diff([0;lokal_max_val;0]))<0)>0);
        
        marker1(lokal_max(lokal_max_poz))=true;
    end
end
end



function marker2=detection_union(marker1,envelope,union_samples)

marker1=marker1(:);

union_samples=ceil(union_samples);

if mod(union_samples,2)==0; union_samples=union_samples+1; end
MASK=ones(union_samples,1);
marker1=convn(marker1,MASK,'same')>0; % dilatation
marker1=~logical(convn(~marker1,MASK,'same')); % erosion

marker2=false(size(marker1));
point=[];
point(:,1)=find(diff([0;marker1])>0); % start 
point(:,2)=find(diff([marker1;0])<0); % end 

for i=1:size(point,1)
    [~,maxp]=max(envelope(point(i,1):point(i,2)));
    marker2(point(i,1)+maxp-1)=true;
end

end





function d=filt50Hz(d,fs,hum_fs)
global bandwidth;

if nargin<3
    hum_fs=50;
end

if min(size(d))==1
   d=d(:); 
end

f0 = hum_fs:hum_fs:fs/2; % Hz
f0(f0>1.1*bandwidth(2))=[];


R = 1; r = 0.985;
for i=1:length(f0)
    b = [1 -2*R*cos(2*pi*f0(i)/fs) R*R];
    a = [1 -2*r*cos(2*pi*f0(i)/fs) r*r];

    d=filtfilt(b,a,d);
end
end



function d=filtering(d,fs)

global decimation
global bandwidth
global f_type

if decimation~=200 && f_type==1;
    f_type=2;
    warning('filter type changed to IIR-Butterworth')
end


% bandpass filtering
switch f_type
    case 1
        % IIR-cheby2
        % low pass
        Wp = 2*bandwidth(2)/fs; Ws = 2*bandwidth(2)/fs+ 0.1;
        Rp = 6; Rs = 60;
        [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);
        [bl,al] = cheby2(n,Rs,Ws);
        
        % high pass
        Wp = 2*bandwidth(1)/fs; Ws = 2*bandwidth(1)/fs- 0.05;
        Rp = 6; Rs = 60;
        [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);
        [bh,ah] = cheby2(n,Rs,Ws,'high');
        
        % pp
        
    case 2
        % IIR butterworth
        Wp = 2*bandwidth(2)/fs; Ws = 2*bandwidth(2)/fs + 0.1; Ws(Ws>1)=1;
        Rp = 6; Rs = 60;
        [but_order,Ws] = buttord(Wp,Ws,Rp,Rs);
        [bl,al]=butter(but_order,Ws);
        
        Wp = 2*bandwidth(1)/fs; Ws = 2*bandwidth(1)/ fs-0.05; Ws(Ws<0)=0.1;
        Rp = 6; Rs = 60;
        [but_order,Ws] = buttord(Wp,Ws,Rp,Rs);
        [bh,ah]=butter(but_order,Ws,'high');
        
%         [bl,al]=butter(but_order,2*bandwidth(2)/fs);
%         [bh,ah]=butter(but_order,2*bandwidth(1)/fs,'high');
        

        
    case 3
        % FIR
        % low pass
        bl = fir1(fs/2,2*bandwidth(2)/fs); al=1;
        
        % high pass
        bh = fir1(fs/2,2*bandwidth(1)/fs,'high'); ah=1;
end

% figure(); freqz(bl,al,10*fs,fs); figure(); freqz(bh,ah,10*fs,fs)

% filter testing
hl=freqz(bl,al,10*fs,fs);
hh=freqz(bh,ah,10*fs,fs);
if max(abs(hl))>1.001 || max(abs(hh))>1.001;
    error('filters is probably unstable !!!')
end

d=filtfilt(bh,ah,d);
if bandwidth(2)==fs/2; return; end
d=filtfilt(bl,al,d);

end




function M=beta_detect(d,fs,beta,winsize,beta_AR)


M=[];
winsize=winsize*fs;
noverlap=round(0.5*winsize);

index=1:winsize-noverlap:size(d,1)-winsize+1;
if isempty(index)
   index=1;
   winsize=size(d,2);
end

[bb,aa]=butter(4,2*30/fs);
parfor ch=1:size(d,2)
    MM=[];
%     winsize=5*fs;
%     noverlap=2.5*fs;
    
    for i=1:length(index)
        seg=d(index(i):index(i)+winsize-1,ch);
        seg=filtfilt(bb,aa,seg);
        
        a=lpc(seg-mean(seg),beta_AR);
        [h,f]=freqz(1,a,[],fs);
        h=abs(h);
        poz=f(diff([0; sign(diff([0; h]))])<0);
        
        MM(i)=sum(poz<25 & poz>beta)>0;
    end
    MM=[MM MM(end)];
    M(:,ch)=interp1([index size(d,1)],MM,1:size(d,1),'nearest');
end
M=logical(M);
end
