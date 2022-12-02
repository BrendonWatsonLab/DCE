%% shuffle by local jitter
function [CCS,CCSs,CCSn,FFm]=CCFunc(varargin)
name=bz_BasenameFromBasepath(pwd);
spikeGroupso='all';
spikeso = bz_GetSpikes('spikeGroups',spikeGroupso);
UIDstarto=spikeso.UID(1);
UIDendo=length(spikeso.UID);
if exist([name '.UPDOWNIntervals.mat'],'file')
    ud=load([name '.UPDOWNIntervals.mat']);
    UPIntst=ud.UPInts;
    DNIntst=ud.DNInts;
else
    ud=load([name '.SlowWaves.events.mat']);
    UPIntst=ud.SlowWaves.ints.UP;
    DNIntst=ud.SlowWaves.ints.DOWN;
end

p = inputParser;
addParameter(p,'JW',0.5); %local jitter window
addParameter(p,'binSize',0.01);
addParameter(p,'epoch',[0 inf]);
addParameter(p,'spikeGroups',spikeGroupso);
addParameter(p,'slice',0.125);
addParameter(p,'spikes',spikeso);
addParameter(p,'UIDstart',UIDstarto);
addParameter(p,'UIDend',UIDendo);
addParameter(p,'shuffles',1);
addParameter(p,'UPInts',UPIntst);
addParameter(p,'DNInts',DNIntst);
parse(p,varargin{:})
JW=p.Results.JW;
shuffles=p.Results.shuffles;
binSize = p.Results.binSize;
epoch=p.Results.epoch;
spikeGroups=p.Results.spikeGroups;
spikes=p.Results.spikes;
slice=p.Results.slice;
UPInts=p.Results.UPInts;
DNInts=p.Results.DNInts;
UIDstart=p.Results.UIDstart;
UIDend=p.Results.UIDend;

load([name '.SleepState.states.mat']);

if ~strcmp(spikeGroups,'all')
spikes = bz_GetSpikes('spikeGroups',spikeGroups);
end

spiketimes=[];
spikes_shankID=[];
spikes_cluID=[];
spikes_UID=[];
spikes_all=[];
spikeIDs=[];
spike_pointer=1;
for i=1:length(spikes.times)
    
    spikes_UID(spike_pointer:spike_pointer+length(spikes.times{1,i})-1)=spikes.UID(i);
    spikes_shankID(spike_pointer:spike_pointer+length(spikes.times{1,i})-1)=spikes.shankID(i);
    spikes_cluID(spike_pointer:spike_pointer+length(spikes.times{1,i})-1)=spikes.cluID(i);
    spiketimes=cat(1,spiketimes,spikes.times{1,i});
    spike_pointer=spike_pointer+length(spikes.times{1,i});
end
spikeIDs=[spikes_shankID;spikes_cluID;spikes_UID]';
spikes_all=cat(2,spiketimes,spikeIDs);
spikes_all=sortrows(spikes_all,1);



duration = 3;  %sec
conv_w = .010/binSize;  % 10ms window
alpha = 0.001;

SF={};
maxgroup=UIDend;
CCS=zeros(ceil((duration/binSize)/2)*2+1,maxgroup,maxgroup);
for i=1:size(epoch,1)
                SF{i}=[spikes_all(spikes_all(:,1)>epoch(i,1) & spikes_all(:,1)<epoch(i,2),:)];
    if size(SF{i},1)>1
        CCS=CCS+CCG64(SF{i}(:,1),SF{i}(:,4),'binSize',binSize,'duration',duration,'maxgroup',maxgroup);
    end
end

%%%%%%%%%%%calculate the firing rate, consider epoch

FFm=[];
for m=UIDstart:UIDend
    ct=0;
    for i=1:size(epoch,1)
        ct=ct+length(find(spikes.times{m-UIDstart+1}>epoch(i,1) & spikes.times{m-UIDstart+1}<epoch(i,2)));
    end
    tt=sum(epoch(:,2)-epoch(:,1));
    if tt==Inf
        tt=spiketimes(end)-epoch(1,1);
    end
    FFm(m)=ct/tt;
end
%%%%%%%%%%%%%%
for i=UIDstart:UIDend
    for j=i+1:UIDend
        CCSs(:,i,j)=tri_filter(CCS(:,i,j),JW/binSize);
CCSn(:,i,j)=log10(CCS(:,i,j)./CCSs(:,i,j));
CCSs(:,j,i)=CCSs(:,i,j);
CCSn(:,j,i)=CCSn(:,i,j);
    end
end
CCSn(isnan(CCSn)| isinf(CCSn))=0;