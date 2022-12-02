% shuffle by local jitter (gaussian smooth), CCG with only UP 
function [CCS,CCSs,CCSn,FFm]=CCFuncsUP1dc(varargin)

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



maxgroup=UIDend;
CCS=zeros(ceil((duration/binSize)/2)*2+1,maxgroup,maxgroup);
CCSn=zeros(ceil((duration/binSize)/2)*2+1,maxgroup,maxgroup);
CCSs=zeros(ceil((duration/binSize)/2)*2+1,maxgroup,maxgroup);
blocksize=1000;
parforcell = mat2cell([1:size(UPInts,1)]',diff([0:blocksize:size(UPInts,1)-1,size(UPInts,1)]));
spikes_allst=spikes_all(InIntervals(spikes_all(:,1),UPInts),:);
for n1=UIDstart:UIDend
    for n2=n1+1:UIDend
        spikes_alls=spikes_allst((spikes_allst(:,4)==n1 | spikes_allst(:,4)==n2),:);
        spikes_alls((spikes_alls(:,4)==n1),4)=1;
        spikes_alls((spikes_alls(:,4)==n2),4)=2;
        asCCS=zeros(ceil((duration/binSize)/2)*2+1,2,2,size(parforcell,1));
        parfor i=1:size(parforcell,1)
            for j=1:size(parforcell{i},1)
                k=(i-1)*blocksize+j;
                asSF=[spikes_alls((spikes_alls(:,1)>UPInts(k,1) & spikes_alls(:,1)<UPInts(k,2)),:)];
                if (size(asSF,1)>1)
                    asCCS(:,:,:,i)=asCCS(:,:,:,i)+CCG64(asSF(:,1),asSF(:,4),'binSize',binSize,'duration',duration,'maxgroup',2);
                end
            end
        end
       
        temp=sum(asCCS,4);
        CCS(:,n1,n2)=temp(:,1,2);
        CCS(:,n2,n1)=temp(:,1,2);
        CCSs(:,n1,n2)=tri_filter(temp(:,1,2),JW/binSize);
        CCSs(:,n2,n1)=CCSs(:,n1,n2);
    end
end

for i=UIDstart:UIDend
    for j=UIDstart:UIDend
CCSn(:,i,j)=log10(CCS(:,i,j)./CCSs(:,i,j));
    end
end
CCSn(isnan(CCSn)| isinf(CCSn))=0;

%%Firing rate in upstate
FFm=[];
for m=UIDstart:UIDend
    ct=0;
    for i=1:size(UPInts,1)
        ct=ct+length(find(spikes.times{m-UIDstart+1}>UPInts(i,1) & spikes.times{m-UIDstart+1}<UPInts(i,2)));
    end
    tt=sum(UPInts(:,2)-UPInts(:,1));
    if tt==Inf
        tt=spiketimes(end)-UPInts(1,1);
    end
    FFm(m)=ct/tt;
end
