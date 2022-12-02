%used for local jitter shuffled CCG
function [dce,tce,dcetableform]=ttestdce(CCSo,CCSn,varargin)
UIDend=size(CCSn,2);
binTotNum=size(CCSn,1);
p = inputParser;
addParameter(p,'binSize',0.001);
addParameter(p,'UIDstart',1);
addParameter(p,'subset',[1:UIDend]);
addParameter(p,'sig',1e-20);
addParameter(p,'zerosth',1/3);
addParameter(p,'ccgBaselineTh',5);
addParameter(p,'slice',0.125); %sec
addParameter(p,'edgeDelay',floor(binTotNum/2)); % in bin number
parse(p,varargin{:})
binSize = p.Results.binSize;
UIDstart = p.Results.UIDstart;
fff = p.Results.subset;
sig = p.Results.sig;
zerosth=p.Results.zerosth;
ccgBaselineTh=p.Results.ccgBaselineTh;
slice=p.Results.slice;
edgeDelay=p.Results.edgeDelay;
dce=zeros(UIDend,UIDend);
tce=zeros(UIDend,UIDend);
for i=UIDstart:UIDend
    for j=i+1:UIDend
        firstNonzeros=find(CCSo(:,i,j)~=0,1,'first');
        lastNonZeros=find(CCSo(:,i,j)~=0,1,'last')-1;
        lastNonZeros=size(CCSo,1)-lastNonZeros;
        NonZeroi=min([firstNonzeros lastNonZeros]);
        NonZeroi=max([NonZeroi floor(binTotNum/2)-edgeDelay]);
        if NonZeroi==0
            NonZeroi=1;
        end
        ccgBaseline=(mean(CCSo(NonZeroi:NonZeroi+ceil(slice/binSize),i,j))+mean(CCSo(end-NonZeroi-ceil(slice/binSize)+2:end-NonZeroi+1,i,j)))/2;
        ccgedge=[CCSn(NonZeroi:ceil(slice/binSize),i,j);CCSn(end-NonZeroi-ceil(slice/binSize)+2:end-NonZeroi+1,i,j)];
        ccgcenter=[CCSn(ceil(size(CCSn,1)/2)-ceil(slice/binSize):ceil(size(CCSn,1)/2)+ceil(slice/binSize),i,j)];
        
        
        [hr,pr]=ttest2(ccgcenter,ccgedge,'Tail','right','Alpha',sig);
        [hl,pl]=ttest2(ccgcenter,ccgedge,'Tail','left','Alpha',sig);
        if (ccgBaseline>ccgBaselineTh & size(find(CCSo(:,i,j)==0),1)<zerosth*size(CCSo,1) & ismember(i,fff) & ismember(j,fff) & i>=UIDstart & j>=UIDstart)
            if (hr==1)
                findextreme=smoothdata(CCSn(:,i,j),'gaussian',0.2/binSize);
                [M,I]=max(findextreme);
                if(I>0.25/binSize & I<length(findextreme)-0.25/binSize)
                    ccgMid=mean(CCSn(I-ceil(0.02/binSize):I+ceil(0.02/binSize),i,j));
                    ccgMidMinusBase=ccgMid-mean(ccgedge);
                    dce(i,j)=ccgMidMinusBase;
                    tce(i,j)=1;
                end
            end
            if (hl==1)
                findextreme=smoothdata(CCSn(:,i,j),'gaussian',0.2/binSize);
                [M,I]=min(findextreme);
                if(I>0.25/binSize & I<length(findextreme)-0.25/binSize)
                    ccgMid=mean(CCSn(I-ceil(0.02/binSize):I+ceil(0.02/binSize),i,j));
                    ccgMidMinusBase=ccgMid-mean(ccgedge);
                    dce(i,j)=ccgMidMinusBase;
                    tce(i,j)=-1;
                end
            end
        end
        if (ccgBaseline<=ccgBaselineTh | size(find(CCSo(:,i,j)==0),1)>1/3*size(CCSo,1) | ~ismember(i,fff) | ~ismember(j,fff) | i<UIDstart | j<UIDstart)
            tce(i,j)=NaN;
        end
    end
end
dcetableform=[];
for j=1:UIDend
    for    k=1:UIDend
        if tce(j,k)==1 | tce(j,k)==-1
            dcetableform(:,end+1)=[j;k;tce(j,k);dce(j,k)];
        end
    end
end
tce=tce+tce';
dce=dce+dce';
tce(1:1+size(tce,1):end) = NaN;