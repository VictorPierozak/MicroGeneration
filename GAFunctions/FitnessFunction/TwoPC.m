%% Function to calculate the Two Point Correlation Function faster than Steve
function [hand]=TwoPC(hand)
MBytes=whos;MBytesA(1)=sum([MBytes.bytes])/1024^2;
startTime=tic;
[si(1),si(2),si(3)]=size(hand.Net_Or);
VoxDim_m=1e-9*[...
    str2double(get(hand.L1box,'String')),...
    str2double(get(hand.L2box,'String')),...
    str2double(get(hand.L3box,'String'))];

%% Pre-process data (pick a phase)

Check(1,1)=get(hand.Check_B1,'value');
Check(1,2)=get(hand.Check_B2,'value');
Check(1,3)=get(hand.Check_B3,'value');
Check(2,1)=get(hand.Check_G1,'value');
Check(2,2)=get(hand.Check_G2,'value');
Check(2,3)=get(hand.Check_G3,'value');
Check(3,1)=get(hand.Check_W1,'value');
Check(3,2)=get(hand.Check_W2,'value');
Check(3,3)=get(hand.Check_W3,'value');

NetVals=unique(hand.Net_Or);

Result=nan(3,3,max(size(hand.Net_Or)));
if length(size(hand.Net_Or))==3
    jdx_max=3;
else
    jdx_max=2;
end

for jdx=1:jdx_max % Directions
    for idx=1:3 % Colors
        if Check(idx,jdx)==1
            if jdx==1
                Net=hand.Net_Or;
            elseif jdx==2
                Net=permute(hand.Net_Or,[2,3,1]);
            else
                Net=permute(hand.Net_Or,[3,1,2]);
            end
            if length(hand.NetVals)==2 && idx==3
                Net=uint8(Net==NetVals(2));
            else
                Net=uint8(Net==NetVals(idx));
            end
            Result(idx,jdx,1:si(jdx))=SConv(Net);
            clear Net
        end
    end
end

%% Conv
    function Result=SConv(Net)
        [s(1),s(2),s(3)]=size(Net);
        for i=1:s(1)
            SumConv(i)=sum(Net(1:end-i+1,:,:).*Net(i:end,:,:),'all');
        end
        
        normVec=[[s(1):-1:1]]*(s(2)*s(3))*mean(Net(:));
        Result=SumConv./normVec;
    end
MBytes=whos;MBytesA(2)=sum([MBytes.bytes])/1024^2;
if length(hand.filename)>53
    hand.fil=[hand.filename(1:50)];
elseif length(hand.filename)>4
    hand.fil=hand.filename(1:find(hand.filename=='.')-1);
else
    hand.fil=['OutputResults'];
end
if ~isstrprop(hand.fil(1), 'alpha')
    hand.fil=['X',hand.fil];
end
if hand.InLineMode==0
    varname=[hand.fil,'.TwoPC'];
    assignin(hand.WoSpace,'temp',Result);
    evalin(hand.WoSpace,[varname,'=temp;']);
    evalin(hand.WoSpace,'clear temp');
else
    Pha='BGW';
    for i=1:3 %Phase
        for j=1:3 % Direction
            if Check(i,j)==1
                varname=['OutputResults.TwoPC','.',Pha(i),num2str(j)];
                eval(['hand.',varname,'= Result(i,j,:);']);
            end
        end
    end
end

if hand.InLineMode==0
    %% reporting
    if toc(startTime)<120
        timeStr=[num2str(ceil(toc(startTime))),' s'];
    elseif toc(startTime)<2*3600
        timeStr=[num2str(ceil(toc(startTime)/60)),' mins'];
    else
        timeStr=[num2str(ceil(toc(startTime)/3600)),' hours'];
    end
    if max(round(MBytesA))<1024
        disp(['Calculation for ',num2str(numel(hand.Net_Or)),' voxels (Memory < ',num2str(roundsf(max(MBytesA),2,'ceil')),' MB, Time < ',timeStr,').',char(10)])
    else
        disp(['Calculation for ',num2str(round(numel(hand.Net_Or)/1000^2)),' million voxels (Memory < ',num2str(roundsf(max(MBytesA)/1024,2,'ceil')),' GB, Time < ',timeStr,').',char(10)])
    end
    
    %% plotting
    
    TwoPCfig=figure(...
        'Name',['TF_2PC: ',hand.fil],...
        'Color',[1 1 1],...
        'renderer','painters',...
        'WindowStyle','normal',...
        'PaperPositionMode','auto',...
        'PaperOrientation','landscape',...
        'units','characters',...
        'position',[20 20 130 40]);
    [hAxis]=subplot(1,1,1);
    [hPlot]=plot(0:length(Result)-1,reshape(Result,9,max(size(hand.Net_Or))));
    set(hPlot,...
        'linewidth',1.5);
    
    set(hPlot(1:3:9),'Color',[0 0 0])
    set(hPlot(2:3:9),'Color',[0 1 0])
    set(hPlot(3:3:9),'Color',[0.5 0.5 .9])
    set(hPlot(1:3),'LineStyle','-')
    set(hPlot(4:6),'LineStyle','--')
    set(hPlot(7:9),'LineStyle',':')
    set(hAxis,...
        'TickLabelInterpreter','latex',...
        'linewidth',1.5,...
        'ylim',[0 1],...
        'fontsize',11);
    xlabel('Distance, $r$ / Voxels','interpreter','latex')
    ylabel('Normalised two point correlation, $S(r)/S(0)$','interpreter','latex')
    LegStr=...
        [{'Black phase, direction 1'},...
        {'Green phase, direction 1'},...
        {'White phase, direction 1'};...
        {'Black phase, direction 2'},...
        {'Green phase, direction 2'},...
        {'White phase, direction 2'};...
        {'Black phase, direction 3'},...
        {'Green phase, direction 3'},...
        {'White phase, direction 3'}]';
    leg=legend(hPlot(logical(Check(:))),LegStr(logical(Check)));
    set(leg,'interpreter','latex')
    if get(hand.Check_pdfSave,'value')==1
        print(TwoPCfig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_TF_2PC'],'-dpdf');
    end
end
end
