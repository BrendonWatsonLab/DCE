cp=struct;
%%
dt = 1/1000; % sec
tend = 500000; % (sec) total spike time
nTrials = 1; % number of trials
binSize=0.001; % sec
duration=3; % sec


%% creating up down states interval from real data
name=bz_BasenameFromBasepath(pwd);
load([name '.UPDOWNIntervals.mat']);
UPIntt=diff(UPInts,1,2);
DNIntt=diff(DNInts,1,2);
InttLength=min(length(UPIntt),length(DNIntt));
Dntemp=reshape([DNIntt(1:InttLength)';UPIntt(1:InttLength)'],[],1);
Dntemp1=[0;cumsum(Dntemp)];
DNIntss=round(1000*reshape(Dntemp1(1:end-1),2,[])');
UPIntss=round(1000*reshape(cumsum(Dntemp),2,[])');
DNIntsss=DNIntss/1000;
UPIntsss=UPIntss/1000;
tend=round(UPIntss(end)/1000);


dcem=zeros(UIDend-UIDstart+1,UIDend-UIDstart+1);
load([name '.FFm_UP.mat'])
for n1=UIDstartc:UIDendc
    parfor n2=n1+1:UIDendc
        %%   up and down states in all SWS epoch
        sp1=spikes.times{n1};
        sp2=spikes.times{n2};
        spupints1={};
        spupints2={};
        for i=1:size(UPInts,1)
            spupints1{i}=sp1(find(sp1>UPInts(i,1) & sp1<UPInts(i,2)))-UPInts(i,1);
            spupints2{i}=sp2(find(sp2>UPInts(i,1) & sp2<UPInts(i,2)))-UPInts(i,1);
        end
        %% percentage of the upstate instead
        block=50;
        hspt1s=zeros(size(UPInts,1),block);hspt2s=zeros(size(UPInts,1),block);
        for it=1:size(UPInts,1)
            [hspt1s(it,:),cspt1]=hist([0;spupints1{it};UPInts(it,2)-UPInts(it,1)],block);
            [hspt2s(it,:),cspt2]=hist([0;spupints2{it};UPInts(it,2)-UPInts(it,1)],block);
            hspt1s(it,1)=hspt1s(it,1)-1;hspt1s(it,end)=hspt1s(it,end)-1;
            hspt2s(it,1)=hspt2s(it,1)-1;hspt2s(it,end)=hspt2s(it,end)-1;
            hspt1s(it,:)=hspt1s(it,:)/sum(hspt1s(it,:));
            hspt2s(it,:)=hspt2s(it,:)/sum(hspt2s(it,:));
        end
        hspt1=mean(hspt1s,1,'omitnan');
        hspt2=mean(hspt2s,1,'omitnan');

        %%

        FFm=FFm_UP;
        spikesn=spikes;
        dcen=[];
        for it=1:1
            fr1=firingRateGen(DNIntss,hspt1,FFm(n1),tend);
            fr2=firingRateGen(DNIntss,hspt2,FFm(n2),tend);

            [spikeMat1, tVec] = poissonSpikeGen(fr1, tend, nTrials);
            [spikeMat2, tVec] = poissonSpikeGen(fr2, tend, nTrials);

            xaxf=[-duration/2:binSize:duration/2];
            sp1=(find(spikeMat1)'-1)/1000;
            sp2=(find(spikeMat2)'-1)/1000;
            spikesn.times=[];
            spikesn.times{1}=sp1;
            spikesn.times{2}=sp2;
            [CCS,CCSs,CCSn,FFmn]=CCFunc('spikes',spikesn,'binSize',binSize1,'JW',JW,'shuffles',shuffles,'UIDstart',1,'UIDend',2);
            [dce,tce,dcetableform]=ttestdce(CCS,CCSn,'subset',[1:2],'sig',1e-2,'UIDstart',1,'binSize',binSize1);
            dcen(it)=dce(1,2);
        end
        dcem(n1,n2)=mean(dcen);

    end
end
poolobj = gcp ('nocreate');
%%
it=1;
for n1=UIDstartc:UIDendc
    for n2=n1+1:UIDendc
        cp(it).dce_SWS=dce_SWS(n1,n2);
        cp(it).dce=dcem(n1,n2);
        cp(it).n1=n1;
        cp(it).n2=n2;
        cp(it).name=name;
        cp(it).dce_WAKE=dce_WAKE(n1,n2);
        cp(it).dce_REM=dce_REM(n1,n2);
        it=it+1;
    end
end
%%
n1=31;
n2=43;
sp1=spikes.times{n1};
sp2=spikes.times{n2};
spupints1={};
spupints2={};
for i=1:size(UPInts,1)
    spupints1{i}=sp1(find(sp1>UPInts(i,1) & sp1<UPInts(i,2)))-UPInts(i,1);
    spupints2{i}=sp2(find(sp2>UPInts(i,1) & sp2<UPInts(i,2)))-UPInts(i,1);
end
block=100;
hspt1s=zeros(size(UPInts,1),block);hspt2s=zeros(size(UPInts,1),block);
for it=1:size(UPInts,1)
    [hspt1s(it,:),cspt1]=hist([0;spupints1{it};UPInts(it,2)-UPInts(it,1)],block);
    [hspt2s(it,:),cspt2]=hist([0;spupints2{it};UPInts(it,2)-UPInts(it,1)],block);
end
hspt1s(:,1)=hspt1s(:,1)-1;hspt1s(:,end)=hspt1s(:,end)-1;
hspt2s(:,1)=hspt2s(:,1)-1;hspt2s(:,end)=hspt2s(:,end)-1;
hspt1=sum(hspt1s,1);
hspt2=sum(hspt2s,1);

load([name '.FFm_UP.mat'])
FFm=FFm_UP;
fr1=firingRateGen(DNIntss,hspt1,FFm(n1),tend);
fr2=firingRateGen(DNIntss,hspt2,FFm(n2),tend);

[spikeMat1, tVec] = poissonSpikeGen(fr1, tend, nTrials);
[spikeMat2, tVec] = poissonSpikeGen(fr2, tend, nTrials);

xaxf=[-duration/2:binSize:duration/2];
sp1=(find(spikeMat1)'-1)/1000;
sp2=(find(spikeMat2)'-1)/1000;
spikes_all1=[sp1 ones(length(sp1),1)];
spikes_all2=[sp2 ones(length(sp2),1)*2];
spikes_allss=sortrows([spikes_all1;spikes_all2],1);
CCSs=CCG64(spikes_allss(:,1),spikes_allss(:,2),'binSize',binSize,'duration',duration,'maxgroup',2);

figure('Position',[0 0 960 240])
subplot(1,3,1)
p1=plot([1:block],hspt1,'LineWidth',2)
hold on
p2=plot([1:block],hspt2,'LineWidth',2)
xlabel('Percentile in UP (%)')
ylabel('Count')
legend({'Neuron 1','Neuron 2'})
title('Spike count histogram')
p1.Color(4) = 0.5;
p2.Color(4) = 0.5;
yticks([0:500:1000])
subplot(1,3,2)
plot(xaxf,smoothdata(CCSs(:,1,2),'gaussian',10));
hold on
v=0;
plot([v v], ylim,'g--');
xlim([-1.5 1.5])
xlabel('Lag (s)')
ylabel('Count')
title('Simulation')
subplot(1,3,3)
xaxft=[-duration/2:binSize1:duration/2];
plot(xaxft,smoothdata(CCS_SWS(:,n1,n2),'gaussian',10));
xlim([-1.5 1.5])
hold on
v=0;
plot([v v], ylim,'g--');
xlabel('Lag (s)')
ylabel('Count')
title('nonREM')
sgtitle([num2str(n1) '-' num2str(n2)])
yticks([0:1500:3000])
print('Simulation poisson','-djpeg','-r300')
%%
function [fr]=firingRateGen(downti,hc,afr,tend)

fr=zeros(1,tend*1000);
for i=1:size(downti,1)-1
    binunit=round((downti(i+1,1)-1-downti(i,2)-1)/size(hc,2))+1;
    fr0=hc*afr/mean(hc);
    t0=reshape(ones(binunit,length(hc)).*fr0,1,[]);
    fr(downti(i,2)+1:downti(i+1,1)-1)=t0(1:downti(i+1,1)-1-downti(i,2));
end
fr=fr*mean(afr)/mean(fr);
end

function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials)
dt = 1/1000; % s
nBins = floor(tSim/dt);
spikeMat = rand(nTrials, nBins) < fr*dt;
tVec = 0:dt:tSim-dt;
end
