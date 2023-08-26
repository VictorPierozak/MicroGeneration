%%%%%%%%%%%%%%%%%%%%%%   TauFactor

%% Opis zmiennych %%

% -> NetVals - liczba faz %
% -> PhaDir - macierz w formacie
%   x y z
% B
% G
% W
% Jeśli 1 - sprawdza krętość w danym kierunku
% -> VoxDims - wymiary woxela [x,y,z]
% -> Solver mode / hand.Check_TauMode - od 1 do 8 ROZPISZ TRYBY!
% -> RVAmode/hand.Pop_RV - representative volume analysis, 0 is disabled
% -> Net_Or - microstructure


% TauFactor is a MatLab application for quantifying the Tortuosity Factor
% and several other microstructural properties based on image data. 
% It is suitable for 2D or 3D geometry data for materials containing up to 
% 3 distinct phases. 
% For support with this application please contact Sam Cooper at:
% camsooper@gmail.com

% Copyright (c) 2016, Samuel J Cooper
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function varargout = TauFactor(varargin)
%% TAUFACTOR MATLAB code for TauFactor.fig
%      TAUFACTOR, by itself, creates a new TAUFACTOR or raises the existing
%      singleton*.
%
%      H = TAUFACTOR returns the handle to a new TAUFACTOR or the handle to
%      the existing singleton*.
%
%      TAUFACTOR('CALLBACK',hObject,eventData,hand,...) calls the local
%      function named CALLBACK in TAUFACTOR.M with the given input arguments.
%
%      TAUFACTOR('Property','Value',...) creates a new TAUFACTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TauFactor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TauFactor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TauFactor

% Last Modified by GUIDE v2.5 11-Apr-2020 16:22:17

[varargout{1:nargout}] = InLine(varargin(1), varargin(2), varargin(3), varargin(4), varargin(5), varargin(6));

% End initialization code - DO NOT EDIT


% --- Outputs from this function are returned to the command line.
function varargout = TauFactor_OutputFcn(~, ~, hand)
varargout{1} = hand.output;

%% ZMIANY %%
% - Check boxy - zamiana z kontrolek na zmienne 0/1

function [OutputResults]=InLine(~,SolverMode, RVAmode, Net_Or, PhaDir, VoxDims)
%% Assign input variables to the global handle
% Inline mode allows the use to operate the functionsality of TauFactor
% from within other scripts or directly from the command line

% Input variable conversion %

Net_Or = cell2mat(Net_Or);
SolverMode = cell2mat(SolverMode);
RVAmode = cell2mat(RVAmode);
PhaDir = cell2mat(PhaDir);
VoxDims = cell2mat(VoxDims);

% End of conversion %

if nargin<6
    SolverMode=1;
    RVAmode=0;
    if ~exist('PhaDir','var')
        PhaDir=[1,0,0;0,0,0;0,0,0];
    end
    if ~exist('VoxDims','var')
        VoxDims=[1,1,1];
    end
    if ~exist('SolverMode','var')
        SolverMode=1;
    end
end
hand.InLineMode=1;
hand.mexCtrl='mex';
hand.Check_mex= 1; %uicontrol('Parent',[],'Value',1);
try AA=[1:2]';BB=uint32(AA);Mex3DTauIso_mex(2,3,[AA;AA],[AA;AA],1,AA,AA,...
        BB,BB,BB,BB,BB,BB,BB,BB,BB,BB,BB,BB);
catch
   % disp('MEX-file not compatible with the platform');
    hand.mexCtrl='mat';
    hand.Check_mex= 0; %uicontrol('Parent',[],'Value',0);
end
hand.WoSpace='caller';
hand.PhaDir=PhaDir;

%[~,sys] = memory;
%hand.RAM_availableMBs=sys.PhysicalMemory.Available/1024^2;
if length(VoxDims)==1
    VoxDims(2:3)=VoxDims(1);
elseif length(VoxDims)==2
    VoxDims(3)=1;
end
hand.L1box= num2str(VoxDims(1)); %uicontrol('Parent',[],'String',VoxDims(1));
hand.L2box= num2str(VoxDims(2)); %uicontrol('Parent',[],'String',VoxDims(2));
hand.L3box= num2str(VoxDims(3)); %uicontrol('Parent',[],'String',VoxDims(3));

hand.Edit_D_Black= 1; %uicontrol('Parent',[],'String',1);
hand.Edit_D_Green= 0; %uicontrol('Parent',[],'String',0);
hand.Edit_D_White= 0; %uicontrol('Parent',[],'String',0);

hand.Check_TauMode= 1; %uicontrol('Parent',[],'Value',1);

if length(SolverMode)==1
    hand.Check_TauMode = SolverMode; %set(hand.Check_TauMode,'Value',SolverMode);
else
    hand.Check_TauMode = 1;  %set(hand.Check_TauMode,'Value',1);
end
hand.Check_VaryD= 0; %uicontrol('Parent',[],'Value',0);

if RVAmode==0
    hand.Check_RVA= 0; %uicontrol('Parent',[],'Value',0);
    hand.Pop_RV= 1; %uicontrol('Parent',[],'Value',1);
else
    hand.Check_RVA= 1; %uicontrol('Parent',[],'Value',1);
    hand.Pop_RV= RVAmode; %uicontrol('Parent',[],'Value',RVAmode);
end

hObject= 0; %figure('Visible','off');
if ~exist('Name','var')
    Name='OutputResults';
end
hand.filename=[Name,'.zzz'];
hand.impCheck=0;
hand.Blocked=0;
hand.Check_FreqPlots=0;
hand.MBytesA=0;
hand.MemLoc={''};
hand.Aniso=0;
hand.NetVals=unique(Net_Or);


if length(SolverMode)==1
%     if SolverMode==2 %Dmap
%         hand.Net_Or=Net_Or<0;
%         hand.Dmap=Net_Or;
%         set(hand.Check_VaryD,'value',1)
%         set(hand.Edit_D_Black,'String',1000)
%         set(hand.Edit_D_Green,'String','Dmap')
%         set(hand.Edit_D_White,'String',0)
%     else % conventional
        if length(hand.NetVals)==1
            if hand.NetVals==0
                hand.Net_Or=false(size(Net_Or));
            elseif hand.NetVals==2
                hand.Net_Or=unit8(2*ones(size(Net_Or),'uint8'));
            else
                hand.Net_Or=ones(size(Net_Or),'logical');
            end
        elseif length(hand.NetVals)==2
            hand.Net_Or=Net_Or==hand.NetVals(2);
        elseif length(hand.NetVals)==3
            Net=uint8(Net_Or==hand.NetVals(2));
            Net(Net_Or==hand.NetVals(3))=2;
            hand.Net_Or=Net;
        else
            hand.Net_Or=Net_Or<0;
            hand.Dmap=Net_Or;
            hand.Check_VaryD = 1;  %set(hand.Check_VaryD,'value',1)
            hand.Edit_D_Black = '0';  %set(hand.Edit_D_Black,'String',0)
            hand.Edit_D_Green = 'Dmap';  %set(hand.Edit_D_Green,'String','Dmap')
            hand.Edit_D_White = '0'; %set(hand.Edit_D_White,'String',0)
        end
%     end
elseif length(SolverMode)>1 && length(SolverMode)==length(hand.NetVals) %% vary D by phase
    hand.Net_Or=Net_Or;
    if length(SolverMode)==2
        hand.Check_VaryD = 1;  %set(hand.Check_VaryD,'value',1)
        hand.Edit_D_Black = num2str(SolverMode(1));  %set(hand.Edit_D_Black,'String',SolverMode(1))
        hand.Edit_D_Green = '0';  %set(hand.Edit_D_Green,'String',0)
        hand.Edit_D_White = num2str(SolverMode(2)); %set(hand.Edit_D_White,'String',SolverMode(2))
    elseif length(SolverMode)==3
        hand.Edit_D_Black = num2str(SolverMode(1));  %set(hand.Edit_D_Black,'String',SolverMode(1))
        hand.Edit_D_Green = num2str(SolverMode(2)); %set(hand.Edit_D_Green,'String',SolverMode(2))
        hand.Edit_D_White = num2str(SolverMode(3));  %set(hand.Edit_D_White,'String',SolverMode(3))
    else
        return;
    end
else
    disp('Too many phases for this mode. Try loading a Dmap.')
end
clear Net Net_Or=1
%% !!! %%
%hand.Check_Reverse=uicontrol('Parent',[],'Visible','off');
%hand.Check_Blocked=uicontrol('Parent',[],'Visible','off');
%hand.Calculate_Button=uicontrol('Parent',[],'Visible','off');
%hand.TextB_Black=uicontrol('Parent',[],'Visible','off');
%hand.TextB_White=uicontrol('Parent',[],'Visible','off');
%hand.TextB_Green=uicontrol('Parent',[],'Visible','off');

%hand.Check_FluxMapSave=uicontrol('Parent',[],'Visible','off');
hand.Check_B1= PhaDir(1,1); %uicontrol('Parent',[],'Value',PhaDir(1,1));
hand.Check_B2= PhaDir(1,2); %uicontrol('Parent',[],'Value',PhaDir(1,2));
hand.Check_B3= PhaDir(1,3); %uicontrol('Parent',[],'Value',PhaDir(1,3));
hand.Check_W1= PhaDir(3,1); %uicontrol('Parent',[],'Value',PhaDir(3,1));
hand.Check_W2= PhaDir(3,2); %uicontrol('Parent',[],'Value',PhaDir(3,2));
hand.Check_W3= PhaDir(3,3); %uicontrol('Parent',[],'Value',PhaDir(3,3));
if length(hand.NetVals)>2
    hand.Check_G1= PhaDir(2,1); %uicontrol('Parent',[],'Value',PhaDir(2,1));
    hand.Check_G2= PhaDir(2,2); %uicontrol('Parent',[],'Value',PhaDir(2,2));
    hand.Check_G3= PhaDir(2,3); %uicontrol('Parent',[],'Value',PhaDir(2,3));
else
    hand.Check_G1= 0; %uicontrol('Parent',[],'Value',0);
    hand.Check_G2= 0; %uicontrol('Parent',[],'Value',0);
    hand.Check_G3= 0; %uicontrol('Parent',[],'Value',0);
end
hand.VFs(1)=mean(mean(mean(hand.Net_Or==0)));
hand.VFs(2)=mean(mean(mean(hand.Net_Or==1)));
hand.VFs(3)=mean(mean(mean(hand.Net_Or==2)));
[hand]=ExpectedTime(hObject, 1, hand);

switch SolverMode
    case 1
        [hand]=Tortuosity(hObject, 1, hand);
    case 2
        [hand]=Tortuosity(hObject, 1, hand);
    case 3
        [hand]=Tortuosity(hObject, 1, hand);
    case 4
        [hand]=Tortuosity(hObject, 1, hand);
    case 5
        [hand]=Tortuosity(hObject, 1, hand);
    case 6
        [hand]=Tortuosity(hObject, 1, hand);
    %% STILL UNDER DEVELOPMENT %%    
    %case 7
        %[hand]=Metrics(hObject, 1, hand);
    %case 8
     %   [hand]=TwoPC(hand);
end
% If all the simulations were on non percoalting volumes, no results with
% be generated, so a dummy variable is required
if ~isfield(hand,'OutputResults')
    hand.OutputResults.Reply=0;
    OutputResults=hand.OutputResults;
else
    OutputResults=hand.OutputResults;
end
try
 %   close gcf
catch
end

%% END OF THE MAIN FUNCTION %%

function [hand] = VoxelSizeFun(hObject, eventdata, hand)
%% Voxel dimension comparison function
[hand]=ExpectedTime(hObject, eventdata, hand);
if str2double(get(hand.L1box,'String'))~=str2double(get(hand.L2box,'String')) ||...
        str2double(get(hand.L2box,'String'))~=str2double(get(hand.L3box,'String'))
    hand.Aniso=1;
else
    hand.Aniso=0;
end

function [hand]=ExpectedTime(hObject, eventdata, hand)
%% Calculate_Button the expected end time of the simulation
%% Default set

%% GUI STATE DISABLED! %%
    % -> SetGUIstate(hObject, eventdata, hand)
    hand = SetGUIComponentsValue(hand);
%hand.Blocked=get(hand.Check_Blocked, 'Value');
%SimsNo=...
%    get(hand.Check_W1,'Value')+...
%    get(hand.Check_W2,'Value')+...
%    get(hand.Check_W3,'Value')+...
%    get(hand.Check_G1,'Value')+...
%    get(hand.Check_G2,'Value')+...
%    get(hand.Check_G3,'Value')+...
%    get(hand.Check_B1,'Value')+...
%    get(hand.Check_B2,'Value')+...
%    get(hand.Check_B3,'Value');

SimsNo= hand.Check_W1 +...
    hand.Check_W2 +...
    hand.Check_W3 +...
    hand.Check_G1 +...
    hand.Check_G2 +...
    hand.Check_G3 +...
    hand.Check_B1 +...
    hand.Check_B2 +...
    hand.Check_B3;

%if length(hand.NetVals)==2
%    set(hand.Check_G1,'Enable','off')
%    set(hand.Check_G2,'Enable','off')
%    set(hand.Check_G3,'Enable','off')
%end
%if length(size(hand.Net_Or))==2
%    set(hand.Check_B3,'Enable','off')
%    set(hand.Check_G3,'Enable','off')
%    set(hand.Check_W3,'Enable','off')
%end

if ~isfield(hand,'Net_Or')
    return
end
waitcoeff=1.6e-08;
if str2double(hand.L1box)~=str2double(hand.L2box) ||...
        str2double(hand.L2box)~=str2double(hand.L3box)
    waitcoeff=1.3*waitcoeff;
    hand.iter_max=max(size(hand.Net_Or))*40;
else
    hand.iter_max=max(size(hand.Net_Or))*30;
end
hand.time_approx=round(SimsNo*numel(hand.Net_Or)*hand.iter_max/4*waitcoeff/86400,2,'significant');

switch hand.Check_TauMode
    case 1
        hand.iter_max=hand.iter_max*1;
    case 2
        hand.iter_max=hand.iter_max*1;
    case 3
        hand.iter_max=hand.iter_max*1;
    case 4
        hand.iter_max=hand.iter_max*14;
        hand.time_approx=hand.time_approx*8;
    case 5
        hand.iter_max=hand.iter_max*60;
        hand.time_approx=hand.time_approx*8;
    case 6 % Electrode tau
        hand.iter_max=hand.iter_max*4;
        hand.time_approx=hand.time_approx*10;
    case 7 % Metrics
        hand.time_approx=round(numel(hand.Net_Or),2,'significant')*1e-10;
    case 8 % 2PC
        hand.time_approx=round(SimsNo*(numel(hand.Net_Or)^(4/3)),2,'significant')*6e-14;
end

if hand.Check_VaryD == 1
    hand.time_approx=hand.time_approx*1.1;
end

if hand.Check_RVA == 1
    hand.time_approx=hand.time_approx*5;
end

if hand.Check_mex == 1
    hand.time_approx=hand.time_approx/2;
end

%% !!! %%
%if isfield(hand,'batch_length')&& get(hand.Radio_BatchStack,'Value')==0
%    hand.time_approx=hand.time_approx*hand.batch_length;
%end

if hand.InLineMode==0
    disp('NO INLINE MODE');
    %[tocStr]=TimeString(hand.time_approx);
    %set(hand.TimeBox,'String',['Time ',tocStr]);
    %[hand]=ExpectedRAM(hObject, eventdata, hand);
end
%guidata(hObject, hand);

%% Used instead of the SetGUIState, it sets same values of fields as they were setted in the SetGUIState %%
function [hand] = SetGUIComponentsValue(hand)

%% TauMode cases
switch hand.Check_TauMode
    
    case 1 %    Tau Factor
        hand.impCheck=0;
        hand.Check_Blocked = 0;
        hand.Check_Reverse = 0;
    case 2 %    Periodic
        hand.impCheck=0;
        hand.Check_RVA = 0;
        hand.Check_VaryD = 0;
        hand.Check_Blocked = 0;
        hand.Check_Reverse = 0;
    case 3 %    Periodic
        hand.impCheck=0;
        hand.Check_RVA = 0;
        hand.Check_VaryD = 0;
        hand.Check_Blocked = 0;
        hand.Check_Reverse = 0;
    case 4 %    Diffusion Impedance
        hand.impCheck=1;
        hand.L1box = '1';
        hand.L2box = '1';
        hand.L3box = '1';
    case 5 %    Symmetrical cell
        %hand.Check_RVA = 0;
        %hand.Check_VaryD = 0;
        %% !!! %%
        %[hand]=SymmetericalPhaseGuesser(hand); 
    case 6 % Electrode tortuosity
        hand.L1box = '1';
        hand.L2box = '1';
        hand.L3box = '1';
    case 7 % Metrics
        hand.Check_FluxMapSave = 0;
    case 8 % 2PC
        hand.Check_RVA = 0;
        hand.Check_FluxMapSave = 0;
end


%% Calculate metrics %%
%% WARNING - NOT CLEARED %%
function [hand]=Metrics(hObject, ~, hand)
    disp('METRICS ARE NOT AVAIABLE');
%if length(hand.filename)>53
%    hand.fil=hand.filename(1:50);
%elseif length(hand.filename)>4
%    hand.fil=hand.filename(1:end-4);
%else
%    hand.fil='OutputResults';
%end
%hand.fil(~isstrprop(hand.fil,'alphanum'))='_';
%if ~isstrprop(hand.fil(1), 'alpha')
%    hand.fil=['X',hand.fil];
%end
%% Call the Metrics RVA function (multiple times if required)
%VoxDim(1)=str2double(get(hand.L1box,'String'));
%VoxDim(2)=str2double(get(hand.L2box,'String'));
%VoxDim(3)=str2double(get(hand.L3box,'String'));

% hand.Net_Or=flipud(rot90(rot90(rot90(hand.Net_Or))));
%hand.Dir=0;
%hand.RVArepeats=0;
%if get(hand.Check_RVA, 'Value')==1 && get(hand.Pop_RV, 'Value')>1
%    hand.RVArepeats=...
%        ((hand.Check_W1 + hand.Check_B1 + hand.Check_G1)>=1)+...
%        ((hand.Check_W2 + hand.Check_B2 + hand.Check_G2 )>=1)+...
%        ((hand.Check_W3 + hand.Check_B3 + hand.Check_G3 )>=1);
%    if ( hand.Check_W1 + hand.Check_B1 + hand.Check_G1 )>=1
%        hand.Dir=1;
%        hand.RVArepeats=hand.RVArepeats-1;
%        [hand]=MetricsRVA(VoxDim(1),VoxDim(2),VoxDim(3), hand.Check_RVA ,hand);
%        drawnow
%    end
%    if (hand.Check_W2 + hand.Check_B2 + hand.Check_G2)>=1
%        hand.Dir=2;
%        hand.RVArepeats=hand.RVArepeats-1;
%        [hand]=MetricsRVA(VoxDim(2),VoxDim(3),VoxDim(1), hand.Check_RVA, hand);
%        drawnow
%    end
%    if (hand.Check_W3 + hand.Check_B3 + hand.Check_G3 )>=1
%        hand.Dir=3;
%        hand.RVArepeats=hand.RVArepeats-1;
%        [hand]=MetricsRVA(VoxDim(3),VoxDim(1),VoxDim(2), hand.Check_RVA, hand);
%    end
%else
%    [hand]=MetricsRVA(VoxDim(1),VoxDim(2),VoxDim(3), hand.Check_RVA,hand);
%end
%guidata(hObject, hand);


%% Main function for tortuosity calculation %%
function [hand]=Tortuosity(hObject, eventdata, hand)
%xxx=tic;
hand.TauSet=zeros(3,3);
if length(hand.filename)>53
    hand.fil=hand.filename(1:50);
elseif length(hand.filename)>4
    hand.fil=hand.filename(1:end-4);
else
    hand.fil='OutputResults';
end
hand.fil(~isstrprop(hand.fil,'alphanum'))='_';
if ~isstrprop(hand.fil(1), 'alpha')
    hand.fil=['X',hand.fil];
end
%% Step through each phase/direction check box to run the tortuosity calculations
if str2double(hand.L1box)~=str2double(hand.L2box) ||...
        str2double(hand.L2box)~=str2double(hand.L3box)
    hand.Error_max=round(1e-4/(max(size(hand.Net_Or)))^1.2,3,'significant');
else
    hand.Error_max=round(1e-4/(max(size(hand.Net_Or)))^1.1,3,'significant');
end
if hand.iter_max<10
    hand.iter_max=10;
end
if hand.Check_TauMode < 4
    hand.check_f=round(hand.iter_max/200);
else
    hand.check_f=round(hand.iter_max/800);
end
hand.JacobiRat=200;
if hand.check_f<10
    hand.check_f=10;
end
if (hand.Check_B1 + hand.Check_B2 + hand.Check_B3 +...
        hand.Check_G1 + hand.Check_G2 + hand.Check_G3 +...
        hand.Check_W1 + hand.Check_W2 + hand.Check_W3)==0
    return
end

if hand.InLineMode==0
    disp('No inline mode!');
    return
end

if hand.Check_VaryD == 0
    if hand.Check_B1 == 1
        hand.Net=hand.Net_Or==0;
        hand.Dir=1;
        hand.Pha='Black';
        [hand]=RVA_Tau(hObject, hand);
        hand.TauSet(1,1)=hand.Results.Tau;
    end
    if hand.Check_B2 == 1
        [~,~,c]=size(hand.Net_Or);
        hand.Net=hand.Net_Or==0;
        if c>1
            hand.Net=logical(permute(hand.Net,[2 3 1]));
        else
            hand.Net=logical(rot90(hand.Net,3));
        end
        hand.Dir=2;
        hand.Pha='Black';
        [hand]=RVA_Tau(hObject, hand);
        hand.TauSet(1,2)=hand.Results.Tau;
    end
    if  hand.Check_B3 == 1
        hand.Net=hand.Net_Or==0;
        hand.Net=logical(flip(permute(hand.Net,[3 1 2]),3));
        hand.Net_Ag = hand.Net; %eSCM
        hand.Net_Cu = flip(permute(hand.Net,[3 1 2]),3); %eSCM
        hand.Dir=3;
        hand.Pha='Black';
        [hand]=RVA_Tau(hObject, hand);
        hand.TauSet(1,3)=hand.Results.Tau;
    end
    
    if  hand.Check_G1 == 1
        hand.Net=hand.Net_Or==1;
        hand.Dir=1;
        hand.Pha='Green';
        [hand]=RVA_Tau(hObject, hand);
        hand.TauSet(2,1)=hand.Results.Tau;
    end
    if  hand.Check_G2 == 1
        [~,~,c]=size(hand.Net_Or);
        hand.Net=hand.Net_Or==1;
        if c>1
            hand.Net=permute(hand.Net,[2 3 1]);
        else
            hand.Net=logical(rot90(hand.Net,3));
        end
        hand.Dir=2;
        hand.Pha='Green';
        [hand]=RVA_Tau(hObject, hand);
        hand.TauSet(2,2)=hand.Results.Tau;
    end
    if  hand.Check_G3 == 1
        hand.Net=hand.Net_Or==1;
        hand.Net=flip(permute(hand.Net,[3 1 2]),3);
        hand.Dir=3;
        hand.Pha='Green';
        [hand]=RVA_Tau(hObject, hand);
        hand.TauSet(2,3)=hand.Results.Tau;
    end
    
    if  hand.Check_W1 == 1
        if max(hand.Net_Or(:))>0
            hand.Net=hand.Net_Or==max(hand.Net_Or(:));
        else
            hand.Net=zeros(size(hand.Net_Or));
        end
        hand.Dir=1;
        hand.Pha='White';
        [hand]=RVA_Tau(hObject, hand);
        hand.TauSet(3,1)=hand.Results.Tau;
    end
    if  hand.Check_W2 == 1
        [~,~,c]=size(hand.Net_Or);
        if max(hand.Net_Or(:))>0
            hand.Net=hand.Net_Or==max(hand.Net_Or(:));
        else
            hand.Net=zeros(size(hand.Net_Or));
        end
        if c>1
            hand.Net=logical(permute(hand.Net,[2 3 1]));
        else
            hand.Net=logical(rot90(hand.Net,3));
        end
        hand.Dir=2;
        hand.Pha='White';
        [hand]=RVA_Tau(hObject, hand);
        
        hand.TauSet(3,2)=hand.Results.Tau;
    end
    if  hand.Check_W3 == 1
        if max(hand.Net_Or(:))>0
            hand.Net=hand.Net_Or==max(hand.Net_Or(:));
        else
            hand.Net=zeros(size(hand.Net_Or));
        end
        hand.Net=logical(flip(permute(hand.Net,[3 1 2]),3));
        hand.Dir=3;
        hand.Pha='White';
        [hand]=RVA_Tau(hObject, hand);
        hand.TauSet(3,3)=hand.Results.Tau;
    end

    if hand.InLineMode==0
        disp('NO INLINE MODE!');
        return;
    end

else % Variable D
    hand.Pha='Multi';
    hand.D_Strings{1}= hand.Edit_D_Black;
    hand.D_Strings{2}= hand.Edit_D_Green;
    hand.D_Strings{3}= hand.Edit_D_White;
    if length(hand.NetVals)<4
        if length(hand.NetVals)<3 % If two phases
            hand.Ds=[str2double(hand.D_Strings(1)), str2double(hand.D_Strings(3))];
        else
            hand.Ds=[str2double(hand.D_Strings(1)), str2double(hand.D_Strings(2)), str2double(hand.D_Strings(3))];
        end
        if ~isnan(sum(hand.Ds)) && max(hand.Ds)==0 % All Ds zero
            disp('At least one phase needs to have a non-zero D value.')
            return
        elseif length(unique(hand.Ds))==1
            disp('If all phases have the same D value, the system becomes analytical.')
            return
        else
            %disp([   'D_black = ',get(hand.Edit_D_Black,'String'),...
            %    '  || D_green = ', get(hand.Edit_D_Green,'String'),...
            %    '  || D_white = ', get(hand.Edit_D_White,'String')]);
            hand.Dmap=hand.Ds(1)*double(hand.Net_Or==0);
            hand.Dmap(hand.Net_Or==1)=hand.Ds(2);
            if length(hand.NetVals)==3
                hand.Dmap(hand.Net_Or==2)=hand.Ds(3);
            end
            if sum(isnan(hand.Ds))>0
                for i=1:length(hand.NetVals)
                    if isnan(hand.Ds(i))
                        Dist=1e-9*mean([str2double(hand.L1box),...
                            str2double(hand.L2box),...
                            str2double(hand.L3box)]);
                        if i==2 && length(hand.NetVals)==2
                            Str=hand.D_Strings{3};
                        else
                            Str=hand.D_Strings{i};
                        end
                        idx=find(Str=='v');
                        D=str2double(Str(1:idx-1));
                        v=str2double(Str(idx+1:end));
                        hand.Dmap(hand.Net_Or==(i-1))=0;
                        DM=bwdist(hand.Net_Or~=(i-1));
                        DM(DM>0)=1./(1/D+1./((DM(DM>0)-0.5)*Dist*v*sqrt(32/(3*pi))));
                        hand.Dmap=hand.Dmap+DM;
                    end
                end
            end
            hand.Dmap=padarray(hand.Dmap,[1,1,1],0);
            hand.D=mean(mean(mean(hand.Dmap(2:end-1,2:end-1,2:end-1))));
            hand.Dmap=1./hand.Dmap;
            if min(hand.Ds)>0
                [~,~,c]=size(hand.Net_Or);
                hand.Net=ones(size(hand.Net_Or),'uint8');
                if hand.Check_B1 == 1
                    hand.Dir=1;
                    [hand]=RVA_Tau(hObject, hand);
                end
                if hand.Check_B2 ==1
                    hand.Dir=2;
                    if c>1
                        hand.Net=logical(permute(ones(size(hand.Net_Or),'uint8'),[2 3 1]));
                        hand.Dmap=(permute(hand.Dmap,[2 3 1]));
                    else
                        hand.Net=logical(rot90(ones(size(hand.Net_Or),'uint8'),3));
                        hand.Dmap=(rot90(hand.Dmap,3));
                    end
                    [hand]=RVA_Tau(hObject, hand);
                    if hand.Check_B3 ==1
                        if c>1 %% undo Dmap adjustment
                            hand.Dmap=(permute(hand.Dmap,[3 1 2]));
                        else
                            hand.Dmap=rot90(hand.Dmap,1);
                        end
                    end
                end
                if hand.Check_B3 ==1
                    hand.Dir=3;
                    hand.Net=flip(permute(true(size(hand.Net_Or)),[3 1 2]),3);
                    hand.Dmap=(flip(permute(hand.Dmap,[3 1 2]),3));
                    [hand]=RVA_Tau(hObject, hand);
                end
            else
                NZD=find(hand.Ds);
                hand.Net=zeros(size(hand.Net_Or),'uint8');
                for i=1:length(NZD)
                    hand.Net=hand.Net+uint8(hand.Net_Or==(NZD(i)-1));
                end
                if hand.Check_B1 ==1
                    hand.Dir=1;
                    [hand]=RVA_Tau(hObject, hand);
                end
                if hand.Check_B2 ==1
                    hand.Net=zeros(size(hand.Net_Or),'uint8');
                    for i=1:length(NZD)
                        hand.Net=hand.Net+uint8(hand.Net_Or==(NZD(i)-1));
                    end
                    [~,~,c]=size(hand.Net_Or);
                    if c>1
                        hand.Net=logical(permute(hand.Net,[2 3 1]));
                        hand.Dmap=(permute(hand.Dmap,[2 3 1]));
                    else
                        hand.Net=logical(rot90(hand.Net,3));
                        hand.Dmap=(rot90(hand.Dmap,3));
                    end
                    hand.Dir=2;
                    [hand]=RVA_Tau(hObject, hand);
                    if hand.Check_B3 ==1
                        if c>1
                            hand.Dmap=(permute(hand.Dmap,[3 1 2]));
                        else
                            hand.Dmap=(rot90(hand.Dmap,1));
                        end
                    end
                end
                if hand.Check_B3 ==1
                    hand.Net=zeros(size(hand.Net_Or),'uint8');
                    for i=1:length(NZD)
                        hand.Net=hand.Net+uint8(hand.Net_Or==(NZD(i)-1));
                    end
                    hand.Net=logical(flip(permute(hand.Net,[3 1 2]),3));
                    hand.Dmap=(flip(permute(hand.Dmap,[3 1 2]),3));
                    hand.Dir=3;
                    [hand]=RVA_Tau(hObject, hand);
                end
            end
        end
    else
        if hand.InLineMode==0
            disp('NO INLINE MODE!')
            return;
        else
        end
        hand.Net=hand.Dmap>0;
        hand.D=mean(hand.Dmap(:));
        hand.Dmap=double(padarray(hand.Dmap,[1,1,1],0));
        hand.Dmap=1./hand.Dmap;
        if hand.Check_B1 == 1
            hand.Dir=1;
            [hand]=RVA_Tau(hObject, hand);
        end
        if  hand.Check_B2 == 1
            [~,~,c]=size(hand.Net_Or);
            if c>1
                hand.Net=logical(permute(hand.Dmap(2:end-1,2:end-1,2:end-1),[2 3 1]));
                hand.Dmap=(permute(hand.Dmap,[2 3 1]));
            else
                hand.Net=logical(rot90(hand.Dmap(2:end-1,2:end-1,2:end-1),3));
                hand.Dmap=(rot90(hand.Dmap,3));
            end
            hand.Dir=2;
            [hand]=RVA_Tau(hObject, hand);
            if hand.Check_B3 ==1
                if c>1
                    hand.Dmap=(permute(hand.Dmap,[3 1 2]));
                else
                    hand.Dmap=(rot90(hand.Dmap,1));
                end
            end
        end
        if hand.Check_B3 == 1
            hand.Net=logical(flip(permute(hand.Dmap(2:end-1,2:end-1,2:end-1),[3 1 2]),3));
            hand.Dmap=(flip(permute(hand.Dmap,[3 1 2]),3));
            hand.Dir=3;
            [hand]=RVA_Tau(hObject, hand);
        end
    end
end
if hand.InLineMode==0
    disp('NO INLINE MODE!')
    return;
end

function [hand]=PrintTauResults(hand)
%% Print the results from each direction
if hand.InLineMode==0 && hand.Check_TauMode ~=5
    if hand.Results.Tau~=inf
        disp([...
            'Tau Factor = ',num2str(round(hand.Results.Tau(end),4,'significant')),' in direction ',num2str(hand.Dir),' of the ',lower(hand.Pha),' phase.'... %,char(10),...
            ]);
    else
        disp([...
            'Tau Factor = inf in direction ',num2str(hand.Dir),' of the ',lower(hand.Pha),' phase.'... %,char(10),...
            ]);
    end
elseif hand.InLineMode==0 && get(hand.Check_TauMode, 'Value')==5  % eSCM
    
end


function [hand]=RVA_Tau(hObject, hand)
%% Prepare subvolume of data for simulation is necessary
if hand.Check_RVA == 1
    hand.Net_full=hand.Net;
    if hand.Check_VaryD == 1
        hand.Dmap_full=hand.Dmap(2:end-1,2:end-1,2:end-1);
    end
    steps=10;
    shrinkStep=1/steps;
    for i=steps:-1:1
        hand.shrink=(i-1)*shrinkStep;%volume percent of shrinkage
        switch hand.Pop_RV
            case 1
                [hand]=RVA_Net_Cube(hand);
            case 2
                [hand]=RVA_Net_Lconst(hand);
            case 3
                [hand]=RVA_Net_Aconst(hand);
            case 4
                [hand]=RVA_Net_Aconst(hand);
        end
        hand.Net=hand.Net_full(...
            hand.RVAdims(1,1):hand.RVAdims(1,2),...
            hand.RVAdims(2,1):hand.RVAdims(2,2),...
            hand.RVAdims(3,1):hand.RVAdims(3,2));
        if hand.Check_VaryD == 1
            hand.Dmap=padarray(hand.Dmap_full(...
                hand.RVAdims(1,1):hand.RVAdims(1,2),...
                hand.RVAdims(2,1):hand.RVAdims(2,2),...
                hand.RVAdims(3,1):hand.RVAdims(3,2)),[1,1,1],inf);
        end
        hand.NetVol(steps+1-i)=numel(hand.Net)/numel(hand.Net_full);
        [hand]=Initialise(hObject, hand);
        if sum(hand.Net_Perc(:))~=0
            hand.Results.Tau=mean([hand.TauFacTop(end),hand.TauFacBot(end)]);
        else
            hand.Results.Tau=inf;
        end
        hand.Tau_RV(steps+1-i)=hand.Results.Tau;
        hand.Pore_RV(steps+1-i)=hand.VolFrac(end);
        
        if hand.InLineMode==0
            disp('NO INLINE MODE')
            return;
        else
            varname=['OutputResults.TauRVA',num2str( hand.Pop_RV ),'_',num2str(hand.Pha(1)),num2str(num2str(hand.Dir))];
            eval(['hand.',varname,'.NetVol = hand.NetVol;']);
            eval(['hand.',varname,'.VolFrac = hand.Pore_RV;']);
            eval(['hand.',varname,'.Tau = real(hand.Tau_RV);']);
        end
        %%%plotting
        if hand.InLineMode==0
            disp('NO INLINE MODE');
           return;
        end
    end
    hand.Pore_RV=1;
    hand.Tau_RV=1;
    hand.NetVol=1;
    hand.RVfig=1;
else
    [a,b,c]=size(hand.Net);
    ErrorStr='Periodic does not function with there are an odd number of voxels in either of the in-plane directions';
    if hand.Check_TauMode == 2 || hand.Check_TauMode == 3
        if c~=1
            if rem(b,2)==1 || rem(c,2)==1
                disp(ErrorStr)
                return
            end
        elseif c==1
            if rem(b,2)==1
                disp(ErrorStr)
                return
            end
        end
    end
    if hand.Check_TauMode == 3
        if rem(a,2)==1
            disp(ErrorStr)
            return
        end
    end
    
    hand.RVAdims=[...
        1,a;...
        1,b;...
        1,c];
    [hand]=Initialise(hObject, hand);
    if isfield(hand,'VolFrac')
        if hand.InLineMode==0
            disp('NO INLINE MODE');
            return;
        else
            varname=['OutputResults.Tau','_',num2str(hand.Pha(1)),num2str(num2str(hand.Dir))];
            
            eval(['hand.',varname,'.VolFrac = hand.VolFrac;']);
            eval(['hand.',varname,'.Tau = real(hand.Results.Tau);']);
            eval(['hand.',varname,'.Perc = real(hand.VolFrac_Perc);']);
        end
        hand.Results.VolFrac=hand.VolFrac;
        if max(hand.Net_Perc(:))>0
            try
                hand=rmfield(hand,{'Tconv','XFlux','T1Top','T2Top','T1Bot','T2Bot','Map',});
            catch
            end
        end
    end
    hand.TauFacTop=0;
    hand.TauFacBot=0;
    hand.DeffTop=0;
    hand.DeffBot=0;
end
[hand]=PrintTauResults(hand);


function [hand]=Initialise(hObject, hand)
%% First step of preperation for tortuosity calculation in which the ordering of calls is set
%tic
%hand.Results.SimTime=0;
if sum(hand.Check_TauMode)==[4 5]
    hand.impCheck=1;
end
if hand.Check_Reverse ==1
    hand.Net=flipud(hand.Net);
end

hand.metric = 1; % if one want to quantify the dead-end pores, hand.metric = 1 otherwise hand.metric = 0 => only through pores are accounted

[hand]=Percolation(hand);
hand.Net_Perc=logical(hand.Net_Perc);
[a,~,~]=size(hand.Net_Perc);

hand.VolFrac_Perc=sum(hand.Net_Perc(:))/sum(hand.Net(:));
if sum(sum(hand.Net_Perc(end,:,:)))==0 && hand.Check_TauMode ~=5
    %hand.Results.SimTime='Does not percolate in this phase+direction';
    hand.VolFrac=mean(hand.Net(:));
    if hand.VolFrac>0
        hand.Results.Tau=inf;
    else
        hand.Results.Tau=nan;
    end
else
    hand.impCheck=0;
    [hand]=Preparation1(hand);
    if sum(hand.Check_TauMode==[1 2 3 4])
        [hand]=Preparation2(hand);
        [hand]=Preparation3(hand);
    elseif sum(hand.Check_TauMode==[5]) % for eSCM
        [hand]=Preparation2bis(hand);   % same task as Preparation2(hand) but for eSCM
        [hand]=Preparation3bis(hand);   % same task as Preparation3(hand) but for eSCM
    elseif sum(hand.Check_TauMode==[6]) % for tau_e
        [hand]=Preparation2(hand);
        hand.freq=complex(1e-6/a); %needs justification
        hand.TopStim=complex(0);
        hand.BotStim=complex(1);
        [hand]=Preparation3imp(hand);
    end
    switch hand.Check_TauMode
        case 1
            eval(['[hand]=Iterate_',hand.mexCtrl,'(hObject, hand);'])
        case 2 %D:D Periodic
            eval(['[hand]=Iterate_',hand.mexCtrl,'(hObject, hand);'])
        case 3 %P:P Periodic
            eval(['[hand]=Iterate_',hand.mexCtrl,'(hObject, hand);'])
        case 4
            hand.impCheck=0;
            if sum(sum(hand.Map(2,2:end-1,2:end-1)))~=0 && hand.Blocked==0
                eval(['[hand]=Iterate_',hand.mexCtrl,'(hObject, hand);'])
                hand.Tau=mean([hand.TauFacBot(end),hand.TauFacTop(end)]);
            else
                hand.Results.Tau=inf;
                hand.Tau=1;
                hand.Tconv=single(hand.Map);
            end
            hand.impCheck=1;
            hand.D=complex(1);
            hand.delta_x=complex(hand.L_X);
            %         freqChar=2*pi*hand.D/(a*hand.Tau(end)*hand.delta_x)^2
            if hand.PercFlag==1
                hand.freqChar=hand.D/(a*hand.delta_x)^2;%Hz?
            else
                hand.freqChar=hand.D/(hand.Max_Path*hand.delta_x)^2;
            end
            %         if sum(sum(sum(hand.Net_Perc)))==0
            %             freqChar=hand.D/(a*2*hand.delta_x)^2;
            %         else
            %             freqChar=hand.D/(a*hand.TauFacBot(end)*hand.delta_x)^2;
            %         end
            if sum(sum(hand.Map(2,2:end-1,2:end-1)))~=0 && hand.Blocked==0
                hand.y=-4:0.1:11;%-4:0.5:11;%-2:0.5:10
            else
                hand.y=-2:0.1:11;%-4:0.5:11;%-2:0.5:10
            end
            %tic
            %% Initialise figure
            hand.impFig=figure(...
                'Name',['TF_Impedance: ','p',hand.Pha(1),'d',num2str(hand.Dir),'_',hand.fil],...
                'units','characters',...
                'position',[277 35 99 40],...
                'renderer','painters',...
                'Color',[1 1 1]);
            
            %% Frequence steps
            hand.freqSet=complex(hand.freqChar*2.^hand.y);
            hand.ImpedanceBotConv=ones(length(hand.y),1)*nan;
            for freqNo=1:length(hand.y)
                %             hand.y(freqNo)
                hand.freqNo=freqNo;
                hand.freqSet=complex(hand.freqChar*2.^hand.y);
                hand.freq=complex(hand.freqSet(freqNo)); %times 2pi?
                if hand.Check_FreqPlots==1
                    try
                        [hand]=InitiatePlot1(hand);
                    catch
                        disp('Unable to plot.')
                    end
                end
                [hand]=Preparation3imp(hand);
                if hand.Check_FreqPlots==1
                    [hand]=InitiatePlot2(hand,(hand.Tconv(:,:,2)));
                end
                eval(['[hand]=Iterate_',hand.mexCtrl,'(hObject, hand);'])
                try
                    [hand]=ImpPlot(hand);
                catch
                    disp('Unable to plot.')
                end
                if (abs(imag(hand.ImpedanceBotConv(freqNo))/max(real(hand.ImpedanceBotConv(:))))<5*hand.conTol &&...
                        abs(real(hand.ImpedanceBotConv(freqNo))/max(real(hand.ImpedanceBotConv(:))))<5*hand.conTol)...
                        || hand.ConvErr==2
                    hand.ImpedanceBotConv=1; hand.ImpedanceBot=1;
                    hand.ImpedanceTopConv=1; hand.ImpedanceTop=1;
                    hand.TauFacBot=1; hand.TauFacTop=1;
                    hand.DeffBot=1; hand.DeffTop=1;
                    hand.freqSet=1;
                    hand.ConvErr=0;
                    % sound(1,1000);pause(0.02);sound(1,1000);
                    %hand.Results.SimTime=hand.Results.SimTime+toc;
                    [hand]=Saving(hand);
                    return
                end
            end
            [hand]=Saving(hand);
        case 5 %Symmetrical
            hand.delta_x=complex(hand.L_X);
            hand.C_dl=str2double(get(hand.edit_capacitance,'string'));
            hand.kappa=str2double(get(hand.edit_conductivity,'string'));
            
            hand.impCheck=1;
            
            hand.Results.Tau=nan;
            %             hand.Tau=1;
            hand.Tconv=single(hand.Map);
            
            %             if hand.PercFlag==1
            %                 hand.freqChar=hand.D/(a*hand.delta_x)^2;%Hz?
            %             else
            %                 hand.freqChar=hand.D/(hand.Max_Path*hand.delta_x)^2;
            %             end
            
            
            %% Initialise figure
            hand.impFig=figure(...
                'Name',['SymmetricCell_Impedance: ','p',hand.Pha(1)],...
                'units','characters',...
                'position',[277 35 99 40],...
                'renderer','painters',...
                'Color',[1 1 1]);
            
            %% Frequence steps
            hand.y=9:-0.1:5; % freq range
            hand.ImpedanceConv=ones(length(hand.y),1)*nan;
            rangef = length(hand.y) ;
            for freqNo=1:rangef
                %                 hand1 = cell2struct(hand13,name') ;
                %                 hand1.freqNo=freqNo;
                %                 hand1.freqSet=complex(10.^hand1.y);
                %                 hand1.freq=complex(hand1.freqSet(freqNo));
                %                 hand.Tconv=single(hand.Map);
                hand.freqNo=freqNo;
                hand.freqSet=complex(10.^hand.y);
                hand.freq=complex(hand.freqSet(freqNo));
                
                if hand.Check_FreqPlots==1
                    try
                        [hand]=InitiatePlot1(hand);
                    catch
                    end
                end
                [hand]=Preparation3imp_sym(hand);
                
                if hand.Check_FreqPlots==1
                    [hand]=InitiatePlot2(hand,(hand.Tconv(:,:,2)),hand.Map);
                end
                
                [hand]=Iterate_sym_mat(hObject,hand); % Iterations
                try
                    [hand]=ImpPlot(hand);
                catch
                    disp('Unable to plot.')
                end
            end
            [hand]=Saving(hand);
        case 6
            %             hand.mexCtrl='mat';
            eval(['[hand]=Iterate_',hand.mexCtrl,'(hObject, hand);'])
    end
    if isfield(hand,'NN_aV')
        hand=rmfield(hand,'NN_aV');
    end
    if isfield(hand,'R')
        hand=rmfield(hand,'R');
    end
end
%[hand]=MBytesRecord(hand,whos); %Memory

function [hand]=Saving(hand)
if get(hand.Check_pdfSave,'Value')==1
    if get(hand.Check_TauMode, 'Value')==5
        Zel = hand.ImpedanceConv ;
        Phase = 180*angle(Zel)/pi ;
        freqset = complex(10.^hand.y);
        Name=[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',hand.Pha(1)];
        save(Name,'freqset','Zel','Phase');
    elseif sum(get(hand.Check_TauMode, 'Value')==[4])
        Name=[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',hand.Pha(1),'d',num2str(hand.Dir)];
        GifWait=2/length(hand.y);
        if get(hand.Check_RVA, 'Value')==0
            print(hand.impFig,[Name,'_Imp'],'-dpdf');
            savefig(hand.impFig,[Name,'_Imp']);
        else
            print(hand.impFig,[Name,'_Imp_RVA',num2str(round(100*hand.NetVol(end)))],'-dpdf');
            savefig(hand.impFig,[Name,'_Imp','_Imp_RVA',num2str(round(100*hand.NetVol(end)))]);
        end
        print(hand.impFig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',hand.Pha(1),'d',num2str(hand.Dir),'_Imp'],'-dpdf');
        hand.im=0;
        for i=1:hand.freqNo
            if i==1
                try
                    imwrite(hand.impFigGif1(:,:,i),hand.cccmap2,[Name,'plot.gif'],'gif','DelayTime',GifWait,'LoopCount',inf)
                catch
                    disp('Unable to write gif.')
                end
                try
                    imwrite(hand.impFigGif2(:,:,i),hand.cccmap2,[Name,'real.gif'],'gif','DelayTime',GifWait,'LoopCount',inf)
                catch
                    disp('Unable to write gif.')
                end
            else
                try
                    imwrite(hand.impFigGif1(:,:,i),[Name,'plot.gif'],'gif','WriteMode','append','DelayTime',GifWait)
                catch
                    disp('Unable to write gif.')
                end
                try
                    imwrite(hand.impFigGif2(:,:,i),[Name,'real.gif'],'gif','WriteMode','append','DelayTime',GifWait)
                catch
                   disp('Unable to write gif.')
                end
            end
        end
    end
end

function [hand]=Percolation(hand)
[a,b,c]=size(hand.Net);
if sum(hand.Check_TauMode==[1 2 3 4 6])
    if sum(hand.Check_TauMode==[1 4 6])
        Mask=false(size(hand.Net));
        Mask(end,:,:)=hand.Net(end,:,:);
        Net_temp=logical(hand.Net);
    elseif hand.Check_TauMode == 2
        %Side percolation
        if c==1
            Net_temp=[hand.Net,hand.Net,hand.Net];
        else
            Net_temp=[hand.Net,hand.Net,hand.Net];
            Net_temp=cat(3,Net_temp,Net_temp,Net_temp);
        end
        Mask=false(size(Net_temp));
        Mask(end,:,:)=Net_temp(end,:,:);
    elseif hand.Check_TauMode == 3
        if c==1
            Net_temp=[hand.Net,hand.Net,hand.Net,hand.Net];
            Net_temp=cat(1,Net_temp,Net_temp);
        else
            Net_temp=[hand.Net,hand.Net,hand.Net];
            Net_temp=cat(3,Net_temp,Net_temp,Net_temp);
            Net_temp(end+1,:,:)=Net_temp(1,:,:);
        end
        Mask=false(size(Net_temp));
        Mask(end,:,:)=Net_temp(end,:,:);
    end
    DM=bwdistgeodesic(logical(Net_temp),logical(Mask),'cityblock');
    clear Net_temp
    
    hand.Min_Path=min(min(DM(1, :,:)))+1;
    hand.Max_Path=max(DM(isfinite(DM)))+1;
    hand.Net_Perc=isfinite(DM);
    
    if hand.impCheck==0 && hand.Check_TauMode ~=6
        Mask(end,:,:)=0;
        Mask(1,:,:)=hand.Net_Perc(1,:,:);
        if hand.Check_TauMode == 1 || hand.Check_TauMode == 4
            DM=bwdistgeodesic(logical(hand.Net_Perc),logical(Mask),'cityblock');
        elseif hand.Check_TauMode == 2
            DM=bwdistgeodesic(logical(hand.Net_Perc),logical(Mask),'cityblock');
            DM=isfinite(DM);
            if c==1
                DM=single(...
                    DM(:,1:b)+...
                    DM(:,b+1:2*b)+...
                    DM(:,2*b+1:3*b));
            else
                DM=single(...
                    DM(:,1:b      ,1:c      )+...
                    DM(:,1:b      ,c+1:2*c  )+...
                    DM(:,1:b      ,2*c+1:end)+...
                    DM(:,b+1:2*b  ,1:c      )+...
                    DM(:,b+1:2*b  ,c+1:2*c  )+...
                    DM(:,b+1:2*b  ,2*c+1:end)+...
                    DM(:,2*b+1:end,1:c      )+...
                    DM(:,2*b+1:end,c+1:2*c  )+...
                    DM(:,2*b+1:end,2*c+1:end));
            end
            DM(DM==0)=inf;
        elseif hand.Check_TauMode == 3
            DM=bwdistgeodesic(logical(hand.Net_Perc),logical(Mask),'cityblock');
            DM=isfinite(DM);
            if c==1
                DM=single(...
                    DM(1:a      ,1:b      )+...
                    DM(1:a      ,b+1:2*b  )+...
                    DM(1:a      ,2*b+1:3*b)+...
                    DM(1:a      ,3*b+1:end)+...
                    DM(a+1:2*a  ,1:b      )+...
                    DM(a+1:2*a  ,b+1:2*b  )+...
                    DM(a+1:2*a  ,2*b+1:3*b  )+...
                    DM(a+1:2*a  ,3*b+1:end));
            else
                DM=single(...
                    DM(:,1:b      ,1:c      )+...
                    DM(:,1:b      ,c+1:2*c  )+...
                    DM(:,1:b      ,2*c+1:end)+...
                    DM(:,b+1:2*b  ,1:c      )+...
                    DM(:,b+1:2*b  ,c+1:2*c  )+...
                    DM(:,b+1:2*b  ,2*c+1:end)+...
                    DM(:,2*b+1:end,1:c      )+...
                    DM(:,2*b+1:end,c+1:2*c  )+...
                    DM(:,2*b+1:end,2*c+1:end));
                
                Net_temp=[hand.Net;hand.Net;hand.Net];
                Net_temp(end+1,:,:)=Net_temp(1,:,:);
                Mask=false(size(Net_temp));
                Mask(end,:,:)=Net_temp(end,:,:);
                DM2=bwdistgeodesic(logical(Net_temp),logical(Mask),'cityblock');
                hand.Net_Perc=isfinite(DM2);
                Mask(end,:,:)=0;
                Mask(1,:,:)=hand.Net_Perc(1,:,:);
                DM2=bwdistgeodesic(logical(hand.Net_Perc),logical(Mask),'cityblock');
                DM2=DM2(a+1:2*a,:,:);
                DM=single(DM(1:end-1,:,:)+isfinite(DM2));
                clear DM2
            end
            DM(DM==0)=inf;
        end
        
        hand.Net_Perc=isfinite(DM);
    end
    % Mean Accessible Area per unit depth
    hand.MAA = uint16(sort(DM(isfinite(DM))));%
    hand.MAA = find([numel(hand.MAA);diff(hand.MAA);numel(hand.MAA)]);
    %     hand.MAAh=harmmean(diff(hand.MAA));  %Requires ML toolbox
    hand.MAAh=1/(mean(1./diff(hand.MAA)));
    hand.MAAa=mean(diff(hand.MAA));
    if sum(sum(hand.Net_Perc(1,:,:)))==0
        hand.PercFlag=0;
    else
        hand.PercFlag=1;
    end
elseif hand.Check_TauMode == 5
    if strcmp(hand.Edit_D_Black,'Elyt')
        OrPhaseOfInterest=0;
    elseif strcmp(hand.Edit_D_Green,'Elyt')
        OrPhaseOfInterest=1;
    elseif strcmp(hand.Edit_D_White,'Elyt')
        OrPhaseOfInterest=2;
    else
        disp('Choose an electrolyte phase')
    end
    
    if strcmp(hand.Edit_D_Black,'Sep')
        OrPhaseOfSep=0;
    elseif strcmp(hand.Edit_D_Green,'Sep')
        OrPhaseOfSep=1;
    elseif strcmp(hand.Edit_D_White,'Sep')
        OrPhaseOfSep=2;
    else
        disp('Choose an Separator phase')
    end
    
    hand.Net_Perc=zeros(size(hand.Net_Or));
    if length(hand.NetVals)==2
        hand.Net_Perc(hand.Net_Or==0)=1;
        hand.Net_Sep=[];
    else
        hand.Net_Perc(hand.Net_Or==OrPhaseOfInterest)=1;
        hand.Net_Sep=find(hand.Net_Or==OrPhaseOfSep);  % Find for separator
    end
    hand.PercFlag=1;
end
%% Plotting ditance maps with cut-off areas
% GrayPore=1-cat(3,Net_Perc,Net_Perc,Net_Perc);
% imagesc(DM);colormap(spring);
% hold on; h1=subimage(0*GrayPore);axis square;set(h1, 'AlphaData', 1-Net);hold off
% NonPerc=DM==inf;
% NonPerc=cat(3,NonPerc,NonPerc,NonPerc);
% hold on; h2=subimage(0.5*NonPerc);axis square;set(h2, 'AlphaData', NonPerc(:,:,1));hold off
% set(gca, 'XTick', [],'YTick',[]);

function [hand]=Preparation1(hand)
%% Second data preparation step where the checkerboarding and vectorisation of the is performed.
% This is also where the anisotropic weightings are accounted for
hand.VolFrac=sum(hand.Net(:))/numel(hand.Net);
[a,b,c]=size(hand.Net_Perc);
hand.Net=1;
switch hand.Dir
    case 1
        hand.L_X=1e-9*str2double(hand.L1box);
        hand.L_Y=1e-9*str2double(hand.L2box);
        hand.L_Z=1e-9*str2double(hand.L3box);
    case 2
        hand.L_X=1e-9*str2double(hand.L2box);
        hand.L_Y=1e-9*str2double(hand.L3box);
        hand.L_Z=1e-9*str2double(hand.L1box);
    case 3
        hand.L_X=1e-9*str2double(hand.L3box);
        hand.L_Y=1e-9*str2double(hand.L1box);
        hand.L_Z=1e-9*str2double(hand.L2box);
end
if str2double(hand.L1box)~=str2double(hand.L2box) ||...
        str2double(hand.L2box)~=str2double(hand.L3box)
    hand.Aniso=1;
else
    hand.Aniso=0;
end
if hand.Check_VaryD ==0
    front=0.5*(hand.L_X*hand.L_Z/hand.L_Y + hand.L_Z*hand.L_Y/hand.L_X + hand.L_X*hand.L_Y/hand.L_Z)^-1;
    hand.c_X=double(front*hand.L_Y*hand.L_Z/hand.L_X);
    hand.c_Y=double(front*hand.L_Z*hand.L_X/hand.L_Y);
    hand.c_Z=double(front*hand.L_X*hand.L_Y/hand.L_Z);
    
    hand.D=1;
end
hand.Qcv=hand.D*(hand.L_Y*b*hand.L_Z*c)/(hand.L_X*a);
%% !!! %%
%[hand]=MBytesRecord(hand,whos,'Preparation1 end'); %Memory

function [hand]=Preparation2(hand)
[a,b,c]=size(hand.Net_Perc);
% Generate maps of nearest neighbours
hand.Map=logical(padarray(hand.Net_Perc, [1,1,1],0));
if sum(hand.Check_TauMode == [2 3])
    if c==1
        hand.Map(:,1,2)=hand.Map(:,end-1,2);
        hand.Map(:,end,2)=hand.Map(:,2,2);
    else
        hand.Map(:,:,1)=hand.Map(:,:,end-1);
        hand.Map(:,:,end)=hand.Map(:,:,2);
        hand.Map(:,1,:)=hand.Map(:,end-1,:);
        hand.Map(:,end,:)=hand.Map(:,2,:);
    end
    if  hand.Check_TauMode == 3
        if c==1
            hand.Map(1,:,2)=hand.Map(end-1,:,2);
            hand.Map(end,:,2)=hand.Map(2,:,2);
        else
            hand.Map(1,:,:)=hand.Map(end-1,:,:);
            hand.Map(end,:,:)=hand.Map(2,:,:);
        end
    end
end
hand.Net_Perc=1;
% Calculate_Button adjusted map of adjusted nearest neighbours
if hand.Check_VaryD == 0
    hand.NN_a=zeros(size(hand.Map),'double');
    hand.NN_a(2:end-1,2:end-1,2:end-1)=...
        hand.c_X*double((hand.Map(1:end-2,2:end-1,2:end-1)+hand.Map(3:end  ,2:end-1,2:end-1)))+...
        hand.c_Y*double((hand.Map(2:end-1,1:end-2,2:end-1)+hand.Map(2:end-1,3:end  ,2:end-1)))+...
        hand.c_Z*double((hand.Map(2:end-1,2:end-1,1:end-2)+hand.Map(2:end-1,2:end-1,3:end  )));
    if sum(hand.Check_TauMode==[6])
        hand.NW_a=zeros(size(hand.Map),'double');
        if c>1
            hand.NW_a(2:end-1,2:end-1,2:end-1)=...
                hand.c_X*double(2-(hand.Map(1:end-2,2:end-1,2:end-1)+hand.Map(3:end  ,2:end-1,2:end-1)))+...
                hand.c_Y*double(2-(hand.Map(2:end-1,1:end-2,2:end-1)+hand.Map(2:end-1,3:end  ,2:end-1)))+...
                hand.c_Z*double(2-(hand.Map(2:end-1,2:end-1,1:end-2)+hand.Map(2:end-1,2:end-1,3:end  )));
        else
            hand.NW_a(2:end-1,2:end-1,2:end-1)=...
                hand.c_X*double(2-(hand.Map(1:end-2,2:end-1,2:end-1)+hand.Map(3:end  ,2:end-1,2:end-1)))+...
                hand.c_Y*double(2-(hand.Map(2:end-1,1:end-2,2:end-1)+hand.Map(2:end-1,3:end  ,2:end-1)));
        end
        if hand.Aniso==0
            hand.NW_a=hand.NW_a/hand.c_X;
        end
        
    end
    if hand.Blocked==1
        hand.NN_a(end-1,:,:)=double(hand.Map(end-1,:,:)).*(hand.NN_a(end-1,:,:)+(2*hand.c_X));
    else
        if  sum(hand.Check_TauMode ==[1 2 4 5])
            hand.NN_a([2 end-1],:,:)=double(hand.Map([2 end-1],:,:)).*(hand.NN_a([2 end-1],:,:)+(2*hand.c_X));
        elseif sum(hand.Check_TauMode ==[6])
            hand.NN_a([2 end-1],:,:)=double(hand.Map([2 end-1],:,:)).*(hand.NN_a([2 end-1],:,:)+(2*hand.c_X));
        end
    end
    if hand.Aniso==0
        hand.NN_a=hand.NN_a/hand.c_X;
    end
    hand.Dmap=1;
    %[hand]=MBytesRecord(hand,whos,'Preparation2 NN_A'); %Memory
else
    hand.Dmap(:,[1, end],:)=inf;hand.Dmap(:,:,[1, end])=inf;
    hand.Dmap([1, end],2:end-1,2:end-1)=0;% Easier than changing dx
    hand.NN_a=zeros(size(hand.Dmap));
    %[hand]=MBytesRecord(hand,whos,'Dmap'); %Memory
    hand.R.Xm=hand.NN_a;hand.R.Xp=hand.NN_a;
    hand.R.Ym=hand.NN_a;hand.R.Yp=hand.NN_a;
    if c>1
        hand.R.Zm=hand.NN_a;hand.R.Zp=hand.NN_a;
    end
    %[hand]=MBytesRecord(hand,whos,'R.XYZ'); %Memory
    hand.NN_a(2:end-1,2:end-1,2:end-1)=2*hand.Map(2:end-1,2:end-1,2:end-1).*(...
        hand.L_Y*hand.L_Z/hand.L_X*(...
        hand.Map(1:end-2,2:end-1,2:end-1)./(hand.Dmap(1:end-2,2:end-1,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1))+...
        hand.Map(3:end  ,2:end-1,2:end-1)./(hand.Dmap(3:end  ,2:end-1,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1)))...
        +hand.L_Z*hand.L_X/hand.L_Y*(...
        hand.Map(2:end-1,1:end-2,2:end-1)./(hand.Dmap(2:end-1,1:end-2,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1))+...
        hand.Map(2:end-1,3:end  ,2:end-1)./(hand.Dmap(2:end-1,3:end  ,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1)))...
        +hand.L_X*hand.L_Y/hand.L_Z*(...
        hand.Map(2:end-1,2:end-1,1:end-2)./(hand.Dmap(2:end-1,2:end-1,1:end-2)+hand.Dmap(2:end-1,2:end-1,2:end-1))+...
        hand.Map(2:end-1,2:end-1,3:end  )./(hand.Dmap(2:end-1,2:end-1,3:end  )+hand.Dmap(2:end-1,2:end-1,2:end-1))));
    hand.NN_a([2 end-1],2:end-1,2:end-1)=2*hand.Map([2 end-1],2:end-1,2:end-1).*(...
        hand.L_Y*hand.L_Z/hand.L_X*(...
        hand.Map([2 end-2],2:end-1,2:end-1)./(hand.Dmap([1 end-2],2:end-1,2:end-1)+hand.Dmap([2 end-1],2:end-1,2:end-1))+...
        hand.Map([3 end-1],2:end-1,2:end-1)./(hand.Dmap([3 end  ],2:end-1,2:end-1)+hand.Dmap([2 end-1],2:end-1,2:end-1)))...
        +hand.L_Z*hand.L_X/hand.L_Y*(...
        hand.Map([2 end-1],1:end-2,2:end-1)./(hand.Dmap([2 end-1],1:end-2,2:end-1)+hand.Dmap([2 end-1],2:end-1,2:end-1))+...
        hand.Map([2 end-1],3:end  ,2:end-1)./(hand.Dmap([2 end-1],3:end  ,2:end-1)+hand.Dmap([2 end-1],2:end-1,2:end-1)))...
        +hand.L_X*hand.L_Y/hand.L_Z*(...
        hand.Map([2 end-1],2:end-1,1:end-2)./(hand.Dmap([2 end-1],2:end-1,1:end-2)+hand.Dmap([2 end-1],2:end-1,2:end-1))+...
        hand.Map([2 end-1],2:end-1,3:end  )./(hand.Dmap([2 end-1],2:end-1,3:end  )+hand.Dmap([2 end-1],2:end-1,2:end-1))));
    
    
    hand.R.Xm(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_Y*hand.L_Z/hand.L_X./(hand.Dmap(1:end-2,2:end-1,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1));
    hand.R.Xp(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_Y*hand.L_Z/hand.L_X./(hand.Dmap(3:end  ,2:end-1,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1));
    
    hand.R.Ym(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_Z*hand.L_X/hand.L_Y./(hand.Dmap(2:end-1,1:end-2,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1));
    hand.R.Yp(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_Z*hand.L_X/hand.L_Y./(hand.Dmap(2:end-1,3:end  ,2:end-1)+hand.Dmap(2:end-1,2:end-1,2:end-1));
    if c>1
        hand.R.Zm(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_X*hand.L_Y/hand.L_Z./(hand.Dmap(2:end-1,2:end-1,1:end-2)+hand.Dmap(2:end-1,2:end-1,2:end-1));
        hand.R.Zp(2:end-1,2:end-1,2:end-1)= 2*hand.Map(2:end-1,2:end-1,2:end-1).*hand.L_X*hand.L_Y/hand.L_Z./(hand.Dmap(2:end-1,2:end-1,3:end  )+hand.Dmap(2:end-1,2:end-1,2:end-1));
    else
        
        hand.R.Zm=1;
        hand.R.Zp=1;
    end
end

hand.Cheq1.P=zeros([a,b,c],'uint8');
% Build checkerboard
hand.Cheq1.P(1:2:end)=1;
if rem(a,2)==0
    hand.Cheq1.P(:,2:2:end,:)=1-hand.Cheq1.P(:,2:2:end,:);
end
if rem(a*b,2)==0
    hand.Cheq1.P(:,:,2:2:end)=1-hand.Cheq1.P(:,:,2:2:end);
end
% Find checkboard neighbours
% N=North, S=South, E=East, W=West, U=Up, D=down
hand.Cheq2.P=1-hand.Cheq1.P;
hand.Cheq1.P=padarray(hand.Cheq1.P>0,[1,1,1]);
hand.Cheq1.P(hand.Map==0)=0;
if numel(hand.Map)<2^32
    hand.AddClass='uint32';
else
    hand.AddClass='uint64';
end
eval(['hand.Cheq1.P=',hand.AddClass,'(find(hand.Cheq1.P));'])
hand.Cheq1.P_Xm=hand.Cheq1.P-1;
hand.Cheq1.P_Xp=hand.Cheq1.P+1;
hand.Cheq1.P_Ym=hand.Cheq1.P-(a+2);
hand.Cheq1.P_Yp=hand.Cheq1.P+(a+2);

if c>1
    hand.Cheq1.P_Zm=hand.Cheq1.P-(a+2)*(b+2);
    hand.Cheq1.P_Zp=hand.Cheq1.P+(a+2)*(b+2);
else
    eval(['hand.Cheq1.P_Zm=',hand.AddClass,'(1);'])
    hand.Cheq1.P_Zp=hand.Cheq1.P_Zm;
end
%[hand]=MBytesRecord(hand,whos,'Prep2 Cheq'); %Memory

hand.Cheq2.P=padarray(hand.Cheq2.P>0,[1,1,1]);hand.Cheq2.P(hand.Map==0)=0;
eval(['hand.Cheq2.P=',hand.AddClass,'(find(hand.Cheq2.P));'])
hand.Cheq2.P_Xm=hand.Cheq2.P-1;
hand.Cheq2.P_Xp=hand.Cheq2.P+1;
hand.Cheq2.P_Ym=hand.Cheq2.P-(a+2);
hand.Cheq2.P_Yp=hand.Cheq2.P+(a+2);
if c>1
    hand.Cheq2.P_Zm=hand.Cheq2.P-(a+2)*(b+2);
    hand.Cheq2.P_Zp=hand.Cheq2.P+(a+2)*(b+2);
else
    eval(['hand.Cheq2.P_Zm=',hand.AddClass,'(1);'])
    hand.Cheq2.P_Zp=hand.Cheq2.P_Zm;
end

if hand.Check_TauMode ==2 || hand.Check_TauMode ==3
    
    Address=uint32(hand.Map).*reshape(uint32(1:(a+2)*(b+2)*(c+2)),a+2,b+2,c+2);
    hand.Cheq1.P_Ym(ismember(hand.Cheq1.P_Ym,Address(:,1,:)))=...
        hand.Cheq1.P_Ym(ismember(hand.Cheq1.P_Ym,Address(:,1,:)))+(a+2)*b;
    hand.Cheq1.P_Yp(ismember(hand.Cheq1.P_Yp,Address(:,end,:)))=...
        hand.Cheq1.P_Yp(ismember(hand.Cheq1.P_Yp,Address(:,end,:)))-(a+2)*b;
    
    hand.Cheq2.P_Ym(ismember(hand.Cheq2.P_Ym,Address(:,1,:)))=...
        hand.Cheq2.P_Ym(ismember(hand.Cheq2.P_Ym,Address(:,1,:)))+(a+2)*b;
    hand.Cheq2.P_Yp(ismember(hand.Cheq2.P_Yp,Address(:,end,:)))=...
        hand.Cheq2.P_Yp(ismember(hand.Cheq2.P_Yp,Address(:,end,:)))-(a+2)*b;
    
    if c>1
        hand.Cheq1.P_Zm(ismember(hand.Cheq1.P_Zm,Address(:,:,1)))=...
            hand.Cheq1.P_Zm(ismember(hand.Cheq1.P_Zm,Address(:,:,1)))+(a+2)*(b+2)*(c);
        hand.Cheq1.P_Zp(ismember(hand.Cheq1.P_Zp,Address(:,:,end)))=...
            hand.Cheq1.P_Zp(ismember(hand.Cheq1.P_Zp,Address(:,:,end)))-(a+2)*(b+2)*(c);
        
        hand.Cheq2.P_Zm(ismember(hand.Cheq2.P_Zm,Address(:,:,1)))=...
            hand.Cheq2.P_Zm(ismember(hand.Cheq2.P_Zm,Address(:,:,1)))+(a+2)*(b+2)*(c);
        hand.Cheq2.P_Zp(ismember(hand.Cheq2.P_Zp,Address(:,:,end)))=...
            hand.Cheq2.P_Zp(ismember(hand.Cheq2.P_Zp,Address(:,:,end)))-(a+2)*(b+2)*(c);
    end
    
    if hand.Check_TauMode ==3
        hand.PPtop1=uint32(find(ismember(hand.Cheq1.P_Xm,Address(1,:,:))));
        hand.PPtop2=uint32(find(ismember(hand.Cheq2.P_Xm,Address(1,:,:))));
        hand.PPbot1=uint32(find(ismember(hand.Cheq1.P_Xp,Address(end,:,:))));
        hand.PPbot2=uint32(find(ismember(hand.Cheq2.P_Xp,Address(end,:,:))));
        hand.PPmid1=uint32(find(~(...
            ismember(hand.Cheq1.P_Xm,Address(1,:,:))+...
            ismember(hand.Cheq1.P_Xp,Address(end,:,:)))));
        hand.PPmid2=uint32(find(~(...
            ismember(hand.Cheq2.P_Xm,Address(1,:,:))+...
            ismember(hand.Cheq2.P_Xp,Address(end,:,:)))));
        
        hand.Cheq1.P_Xm(ismember(hand.Cheq1.P_Xm,Address(1,:,:)))=...
            hand.Cheq1.P_Xm(ismember(hand.Cheq1.P_Xm,Address(1,:,:)))+a;
        hand.Cheq1.P_Xp(ismember(hand.Cheq1.P_Xp,Address(end,:,:)))=...
            hand.Cheq1.P_Xp(ismember(hand.Cheq1.P_Xp,Address(end,:,:)))-a;
        
        hand.Cheq2.P_Xm(ismember(hand.Cheq2.P_Xm,Address(1,:,:)))=...
            hand.Cheq2.P_Xm(ismember(hand.Cheq2.P_Xm,Address(1,:,:)))+a;
        hand.Cheq2.P_Xp(ismember(hand.Cheq2.P_Xp,Address(end,:,:)))=...
            hand.Cheq2.P_Xp(ismember(hand.Cheq2.P_Xp,Address(end,:,:)))-a;
    end
    clear Address
end

if ~isfield(hand,'MBytesA')
    hand.MBytesA=0;
    hand.MemLoc={''};
end
%[hand]=MBytesRecord(hand,whos,'Prep2 Cheq end'); %Memory
%%
% New method
% Find the index of all nodes on top and bottom of the padded volume
eval(['Top=',hand.AddClass,'([1:a+2:(a+2)*(b+2)*(c+2)]);']);
eval(['Base=',hand.AddClass,'([a+2:a+2:(a+2)*(b+2)*(c+2)]);']);

% Find relations that touch the top or bottom and what index they are
[LIA2t]=ismember(hand.Cheq2.P_Xm,Top);
[LIA2b]=ismember(hand.Cheq2.P_Xp,Base);
eval(['hand.T2Top=',hand.AddClass,'(find(LIA2t));']);
eval(['hand.T2Bot=',hand.AddClass,'(find(LIA2b));']);

% Find all the relations pointing at zeros and point them to the same place
[~,hand.Cheq2.P_Xm]=ismember(hand.Cheq2.P_Xm,hand.Cheq1.P);
eval(['hand.Cheq2.P_Xm=',hand.AddClass,'(hand.Cheq2.P_Xm);']);
eval(['hand.Cheq2.P_Xm(hand.Cheq2.P_Xm==0)=',hand.AddClass,'(length(hand.Cheq1.P)+1);']);
[~,hand.Cheq2.P_Xp]=(ismember(hand.Cheq2.P_Xp,hand.Cheq1.P));
eval(['hand.Cheq2.P_Xp=',hand.AddClass,'(hand.Cheq2.P_Xp);']);
eval(['hand.Cheq2.P_Xp(hand.Cheq2.P_Xp==0)=',hand.AddClass,'(length(hand.Cheq1.P)+1);']);
[~,hand.Cheq2.P_Ym]=(ismember(hand.Cheq2.P_Ym,hand.Cheq1.P));
eval(['hand.Cheq2.P_Ym=',hand.AddClass,'(hand.Cheq2.P_Ym);']);
eval(['hand.Cheq2.P_Ym(hand.Cheq2.P_Ym==0)=',hand.AddClass,'(length(hand.Cheq1.P)+1);']);
[~,hand.Cheq2.P_Yp]=(ismember(hand.Cheq2.P_Yp,hand.Cheq1.P));
eval(['hand.Cheq2.P_Yp=',hand.AddClass,'(hand.Cheq2.P_Yp);']);
eval(['hand.Cheq2.P_Yp(hand.Cheq2.P_Yp==0)=',hand.AddClass,'(length(hand.Cheq1.P)+1);']);
if c>1
    [~,hand.Cheq2.P_Zm]=ismember(hand.Cheq2.P_Zm,hand.Cheq1.P);
    eval(['hand.Cheq2.P_Zm=',hand.AddClass,'(hand.Cheq2.P_Zm);']);
    hand.Cheq2.P_Zm(hand.Cheq2.P_Zm==0)=length(hand.Cheq1.P)+1;
    [~,hand.Cheq2.P_Zp]=ismember(hand.Cheq2.P_Zp,hand.Cheq1.P);
    eval(['hand.Cheq2.P_Zp=',hand.AddClass,'(hand.Cheq2.P_Zp);']);
    hand.Cheq2.P_Zp(hand.Cheq2.P_Zp==0)=length(hand.Cheq1.P)+1;
end
eval(['hand.Cheq2.P_Xp(LIA2b)=',hand.AddClass,'(length(hand.Cheq1.P)+2);']);

[LIA1t]=ismember(hand.Cheq1.P_Xm,Top);
[LIA1b]=ismember(hand.Cheq1.P_Xp,Base);
eval(['hand.T1Top=',hand.AddClass,'(find(LIA1t));']);
eval(['hand.T1Bot=',hand.AddClass,'(find(LIA1b));']);
[~,hand.Cheq1.P_Xm]=(ismember(hand.Cheq1.P_Xm,hand.Cheq2.P));
eval(['hand.Cheq1.P_Xm=',hand.AddClass,'(hand.Cheq1.P_Xm);']);
eval(['hand.Cheq1.P_Xm(hand.Cheq1.P_Xm==0)=',hand.AddClass,'(length(hand.Cheq2.P)+1);']);
[~,hand.Cheq1.P_Xp]=(ismember(hand.Cheq1.P_Xp,hand.Cheq2.P));
eval(['hand.Cheq1.P_Xp=',hand.AddClass,'(hand.Cheq1.P_Xp);']);
eval(['hand.Cheq1.P_Xp(hand.Cheq1.P_Xp==0)=',hand.AddClass,'(length(hand.Cheq2.P)+1);']);
[~,hand.Cheq1.P_Ym]=(ismember(hand.Cheq1.P_Ym,hand.Cheq2.P));
eval(['hand.Cheq1.P_Ym=',hand.AddClass,'(hand.Cheq1.P_Ym);']);
eval(['hand.Cheq1.P_Ym(hand.Cheq1.P_Ym==0)=',hand.AddClass,'(length(hand.Cheq2.P)+1);']);
[~,hand.Cheq1.P_Yp]=(ismember(hand.Cheq1.P_Yp,hand.Cheq2.P));
eval(['hand.Cheq1.P_Yp=',hand.AddClass,'(hand.Cheq1.P_Yp);']);
eval(['hand.Cheq1.P_Yp(hand.Cheq1.P_Yp==0)=',hand.AddClass,'(length(hand.Cheq2.P)+1);']);
if c>1
    [~,hand.Cheq1.P_Zm]=ismember(hand.Cheq1.P_Zm,hand.Cheq2.P);
    eval(['hand.Cheq1.P_Zm=',hand.AddClass,'(hand.Cheq1.P_Zm);']);
    hand.Cheq1.P_Zm(hand.Cheq1.P_Zm==0)=length(hand.Cheq2.P)+1;
    [~,hand.Cheq1.P_Zp]=ismember(hand.Cheq1.P_Zp,hand.Cheq2.P);
    eval(['hand.Cheq1.P_Zp=',hand.AddClass,'(hand.Cheq1.P_Zp);']);
    hand.Cheq1.P_Zp(hand.Cheq1.P_Zp==0)=length(hand.Cheq2.P)+1;
end
eval(['hand.Cheq1.P_Xp(LIA1b)=',hand.AddClass,'(length(hand.Cheq2.P)+2);']);


if hand.Check_TauMode ==3
    % Make the top and bottom reference eachother...
end

function [hand]=Preparation2bis(hand)
[a,b,c]=size(hand.Net_Perc);
% Generate maps of nearest neighbours

hand.Map=logical(padarray(hand.Net_Perc, [1,1,1],0));

Map_sep = zeros([a,b,c]) ;
% Map_NN = hand.Map;

if length(hand.NetVals)==3
    hand.Net_Perc(hand.Net_Sep) = 1;
    Map_sep(hand.Net_Sep) = 1;
end

Map_NN_a=logical(padarray(hand.Net_Perc, [1,1,1],0));

hand.Net_t = ones(size(hand.Net_Perc));
Map_t=logical(padarray(hand.Net_t, [1,1,1],0));

Map_sep = logical(padarray(Map_sep, [1,1,1],0)); % Separator
Map_NN = Map_NN_a-Map_sep;

hand.Net_Perc=1;
hand.Net_t=1;

% Calculate nearest neighbours - separator voxels are considered as solid phase here

hand.NN=zeros(size(Map_NN),'double');
hand.NN(2:end-1,2:end-1,2:end-1)=...
    hand.c_X*double((Map_NN(1:end-2,2:end-1,2:end-1)+Map_NN(3:end  ,2:end-1,2:end-1)))+...
    hand.c_Y*double((Map_NN(2:end-1,1:end-2,2:end-1)+Map_NN(2:end-1,3:end  ,2:end-1)))+...
    hand.c_Z*double((Map_NN(2:end-1,2:end-1,1:end-2)+Map_NN(2:end-1,2:end-1,3:end  )));

% Calculate nearest neighbours adjusted - separator voxels are considered as liquid phase here

hand.NN_a=zeros(size(Map_NN_a),'double');
hand.NN_a(2:end-1,2:end-1,2:end-1)=...
    hand.c_X*double((Map_NN_a(1:end-2,2:end-1,2:end-1)+Map_NN_a(3:end  ,2:end-1,2:end-1)))+...
    hand.c_Y*double((Map_NN_a(2:end-1,1:end-2,2:end-1)+Map_NN_a(2:end-1,3:end  ,2:end-1)))+...
    hand.c_Z*double((Map_NN_a(2:end-1,2:end-1,1:end-2)+Map_NN_a(2:end-1,2:end-1,3:end  )));

% Calculate total nearest neighbours

hand.NN_tot=zeros(size(Map_t),'double');
hand.NN_tot(2:end-1,2:end-1,2:end-1)=...
    hand.c_X*double((Map_t(1:end-2,2:end-1,2:end-1)+Map_t(3:end  ,2:end-1,2:end-1)))+...
    hand.c_Y*double((Map_t(2:end-1,1:end-2,2:end-1)+Map_t(2:end-1,3:end  ,2:end-1)))+...
    hand.c_Z*double((Map_t(2:end-1,2:end-1,1:end-2)+Map_t(2:end-1,2:end-1,3:end)));

hand.NN_tot([2 end-1],:,:)=double(Map_t([2 end-1],:,:)).*(hand.NN_tot([2 end-1],:,:)+(hand.c_X));   % adding the double layer between the liquid phase and current collector at both ends


if hand.Aniso==0
    hand.NN=hand.NN/hand.c_X;
    hand.NN_a=hand.NN_a/hand.c_X;
    hand.NN_tot=hand.NN_tot/hand.c_X;
end

%[hand]=MBytesRecord(hand,whos,'Preparation2bis NN NN_a and NN_tot'); %Memory

hand.Cheq1.P=zeros([a,b,c],'uint8');

% Build checkerboard
hand.Cheq1.P(1:2:end)=1;
if rem(a,2)==0
    hand.Cheq1.P(:,2:2:end,:)=1-hand.Cheq1.P(:,2:2:end,:);
end
if rem(a*b,2)==0
    hand.Cheq1.P(:,:,2:2:end)=1-hand.Cheq1.P(:,:,2:2:end);
end
% Find checkboard neighbours
% N=North, S=South, E=East, W=West, U=Up, D=down
hand.Cheq2.P=1-hand.Cheq1.P;
hand.Cheq1.P=padarray(hand.Cheq1.P>0,[1,1,1]);
hand.Cheq1.P(hand.Map==0)=0;
hand.Cheq1.P=uint32(find(hand.Cheq1.P));
hand.Cheq1.P_Xm=hand.Cheq1.P-1;
hand.Cheq1.P_Xp=hand.Cheq1.P+1;
hand.Cheq1.P_Ym=hand.Cheq1.P-(a+2);
hand.Cheq1.P_Yp=hand.Cheq1.P+(a+2);
if c>1
    hand.Cheq1.P_Zm=hand.Cheq1.P-(a+2)*(b+2);
    hand.Cheq1.P_Zp=hand.Cheq1.P+(a+2)*(b+2);
else
    hand.Cheq1.P_Zm=uint32(1);
    hand.Cheq1.P_Zp=uint32(1);
end
%[hand]=MBytesRecord(hand,whos,'Prep2bis Cheq'); %Memory

hand.Cheq2.P=padarray(hand.Cheq2.P>0,[1,1,1]);hand.Cheq2.P(hand.Map==0)=0;
hand.Cheq2.P=uint32(find(hand.Cheq2.P));
hand.Cheq2.P_Xm=hand.Cheq2.P-1;
hand.Cheq2.P_Xp=hand.Cheq2.P+1;
hand.Cheq2.P_Ym=hand.Cheq2.P-(a+2);
hand.Cheq2.P_Yp=hand.Cheq2.P+(a+2);
if c>1
    hand.Cheq2.P_Zm=hand.Cheq2.P-(a+2)*(b+2);
    hand.Cheq2.P_Zp=hand.Cheq2.P+(a+2)*(b+2);
else
    hand.Cheq2.P_Zm=uint32(1);
    hand.Cheq2.P_Zp=uint32(1);
end

if ~isfield(hand,'MBytesA')
    hand.MBytesA=0;
    hand.MemLoc={''};
end
%[hand]=MBytesRecord(hand,whos,'Prep2bis Cheq end'); %Memory
%%
Map_lptop = logical(padarray(zeros([a,b,c]), [1,1,1],0));  % Top electrode
Map_lpbot = logical(padarray(zeros([a,b,c]), [1,1,1],0)); % Bot electrode
Map_lptop(1:floor(end/2),:,:) = hand.Map(1:floor(end/2),:,:);  % Top electrode
Map_lpbot(floor(end/2)+1:end,:,:) = hand.Map(floor(end/2)+1:end,:,:); % Bot electrode
Map_sep(floor(end/2),:,:) =  hand.Map(floor(end/2),:,:); % Liquid phase in Separator at x=floor(end/2)

Let = find(Map_lptop==1);  % Voxels index of Liquid phase of Top electrode
Leb = find(Map_lpbot==1);  % Voxels index of Liquid phase of Bot electrode
Sep = find(Map_sep==1);    % Voxels index of at where total inoic current within the cell is determined

[hand.Cheq1.Sep]=(ismember(hand.Cheq1.P,Sep));   % liquid phase voxels within separator (allow to calculate the total ionic current within the cell)
[hand.Cheq2.Sep]=(ismember(hand.Cheq2.P,Sep));   % liquid phase voxels within separator (allow to calculate the total ionic current within the cell)

[Cheq1_Let]=(ismember(hand.Cheq1.P,Let));  % member of top electrode
[Cheq1_Leb]=(ismember(hand.Cheq1.P,Leb));  % member of bot electrode

[Cheq2_Let]=(ismember(hand.Cheq2.P,Let)); % member of top electrode
[Cheq2_Leb]=(ismember(hand.Cheq2.P,Leb)); % member of bot electrode

hand.Cheq1.Pb(Cheq1_Let) = uint32(length(hand.Cheq2.P)+1);  % Boundary conditions of solid phase of top electrode
hand.Cheq1.Pb(Cheq1_Leb) = uint32(length(hand.Cheq2.P)+2);  % Boundary conditions of solid phase of bot electrode

hand.Cheq2.Pb(Cheq2_Let) = uint32(length(hand.Cheq1.P)+1);   % Boundary conditions of solid phase of top electrode
hand.Cheq2.Pb(Cheq2_Leb) = uint32(length(hand.Cheq1.P)+2);   % Boundary conditions of solid phase of bot electrode

%[hand]=MBytesRecord(hand,whos,'Prep2bis Cheq Top or Bot'); %Memory


[~,hand.Cheq2.P_Xm]=(ismember(hand.Cheq2.P_Xm,hand.Cheq1.P));
hand.Cheq2.P_Xm=uint32(hand.Cheq2.P_Xm);
hand.Cheq2.P_Xm(hand.Cheq2.P_Xm==0)=uint32(length(hand.Cheq1.P)+1);
[~,hand.Cheq2.P_Xp]=(ismember(hand.Cheq2.P_Xp,hand.Cheq1.P));
hand.Cheq2.P_Xp=uint32(hand.Cheq2.P_Xp);
hand.Cheq2.P_Xp(hand.Cheq2.P_Xp==0)=uint32(length(hand.Cheq1.P)+1);
[~,hand.Cheq2.P_Ym]=(ismember(hand.Cheq2.P_Ym,hand.Cheq1.P));
hand.Cheq2.P_Ym=uint32(hand.Cheq2.P_Ym);
hand.Cheq2.P_Ym(hand.Cheq2.P_Ym==0)=uint32(length(hand.Cheq1.P)+1);
[~,hand.Cheq2.P_Yp]=(ismember(hand.Cheq2.P_Yp,hand.Cheq1.P));
hand.Cheq2.P_Yp=uint32(hand.Cheq2.P_Yp);
hand.Cheq2.P_Yp(hand.Cheq2.P_Yp==0)=uint32(length(hand.Cheq1.P)+1);
if c>1
    [~,hand.Cheq2.P_Zm]=ismember(hand.Cheq2.P_Zm,hand.Cheq1.P);
    hand.Cheq2.P_Zm=uint32(hand.Cheq2.P_Zm);
    hand.Cheq2.P_Zm(hand.Cheq2.P_Zm==0)=length(hand.Cheq1.P)+1;
    [~,hand.Cheq2.P_Zp]=ismember(hand.Cheq2.P_Zp,hand.Cheq1.P);
    hand.Cheq2.P_Zp=uint32(hand.Cheq2.P_Zp);
    hand.Cheq2.P_Zp(hand.Cheq2.P_Zp==0)=length(hand.Cheq1.P)+1;
end

[~,hand.Cheq1.P_Xm]=(ismember(hand.Cheq1.P_Xm,hand.Cheq2.P));
hand.Cheq1.P_Xm=uint32(hand.Cheq1.P_Xm);
hand.Cheq1.P_Xm(hand.Cheq1.P_Xm==0)=uint32(length(hand.Cheq2.P)+1);
[~,hand.Cheq1.P_Xp]=(ismember(hand.Cheq1.P_Xp,hand.Cheq2.P));
hand.Cheq1.P_Xp=uint32(hand.Cheq1.P_Xp);
hand.Cheq1.P_Xp(hand.Cheq1.P_Xp==0)=uint32(length(hand.Cheq2.P)+1);
[~,hand.Cheq1.P_Ym]=(ismember(hand.Cheq1.P_Ym,hand.Cheq2.P));
hand.Cheq1.P_Ym=uint32(hand.Cheq1.P_Ym);
hand.Cheq1.P_Ym(hand.Cheq1.P_Ym==0)=uint32(length(hand.Cheq2.P)+1);
[~,hand.Cheq1.P_Yp]=(ismember(hand.Cheq1.P_Yp,hand.Cheq2.P));
hand.Cheq1.P_Yp=uint32(hand.Cheq1.P_Yp);
hand.Cheq1.P_Yp(hand.Cheq1.P_Yp==0)=uint32(length(hand.Cheq2.P)+1);
if c>1
    [~,hand.Cheq1.P_Zm]=ismember(hand.Cheq1.P_Zm,hand.Cheq2.P);
    hand.Cheq1.P_Zm=uint32(hand.Cheq1.P_Zm);
    hand.Cheq1.P_Zm(hand.Cheq1.P_Zm==0)=length(hand.Cheq2.P)+1;
    [~,hand.Cheq1.P_Zp]=ismember(hand.Cheq1.P_Zp,hand.Cheq2.P);
    hand.Cheq1.P_Zp=uint32(hand.Cheq1.P_Zp);
    hand.Cheq1.P_Zp(hand.Cheq1.P_Zp==0)=length(hand.Cheq2.P)+1;
end

function [hand]=Preparation3(hand)
%% Third preparation step for Tau calculation where the volume is initialised as linear
[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
hand.Area_top=sum(sum(hand.Map(2,2:end-1,2:end-1)));
hand.Area_bot=sum(sum(hand.Map(end-1,2:end-1,2:end-1)));
% Specify a relaxation factor
if hand.Check_TauMode ~=3
    hand.w=(2-(pi)/(a*1.5));
%else
%    hand.w=(2-(pi)/(a*1.5));
end
if hand.Check_VaryD ==1
    %     hand.w=1;
end
% hand.w=1.967;
% hand.w=1;
hand.omw=double(1-hand.w);
% Seperate nearest neighbours into checkerboard and precalculate relaxation
% factor and division
hand.NN_aV.w1=double(hand.w./double(hand.NN_a(hand.Cheq1.P)));
hand.NN_aV.w2=double(hand.w./double(hand.NN_a(hand.Cheq2.P)));
%[hand]=MBytesRecord(hand,whos,'Prep3'); %Memory
if hand.Check_TauMode ==1
    hand.NN_a=0;hand=rmfield(hand,'NN_a');
    hand.NN_tot=0;hand=rmfield(hand,'NN_tot');
end
%[hand]=MBytesRecord(hand,whos,'Prep3 remove NN_a'); %Memory
T=single(padarray(hand.Map(2:end-1,2:end-1,2:end-1),[1,1,1],0));
hand.TopStim=0;
hand.BotStim=1;
if hand.Check_VaryD ==1
    hand.R.Xm1=hand.R.Xm(hand.Cheq1.P).*hand.NN_aV.w1;
    hand.R.Xp1=hand.R.Xp(hand.Cheq1.P).*hand.NN_aV.w1;
    hand.R.Xm2=hand.R.Xm(hand.Cheq2.P).*hand.NN_aV.w2;
    hand.R.Xp2=hand.R.Xp(hand.Cheq2.P).*hand.NN_aV.w2;
    hand.R=rmfield(hand.R,{'Xm','Xp'});
    
    hand.R.Ym1=hand.R.Ym(hand.Cheq1.P).*hand.NN_aV.w1;
    hand.R.Yp1=hand.R.Yp(hand.Cheq1.P).*hand.NN_aV.w1;
    hand.R.Ym2=hand.R.Ym(hand.Cheq2.P).*hand.NN_aV.w2;
    hand.R.Yp2=hand.R.Yp(hand.Cheq2.P).*hand.NN_aV.w2;
    hand.R=rmfield(hand.R,{'Ym','Yp'});
    if c>1
        hand.R.Zm1=hand.R.Zm(hand.Cheq1.P).*hand.NN_aV.w1;
        hand.R.Zp1=hand.R.Zp(hand.Cheq1.P).*hand.NN_aV.w1;
        hand.R.Zm2=hand.R.Zm(hand.Cheq2.P).*hand.NN_aV.w2;
        hand.R.Zp2=hand.R.Zp(hand.Cheq2.P).*hand.NN_aV.w2;
    else
        hand.R.Zm1=1;
        hand.R.Zp1=1;
        hand.R.Zm2=1;
        hand.R.Zp2=1;
    end
    hand.R=rmfield(hand.R,{'Zm','Zp'});
    hand.NN_aV=1;
end
% T(1,2:end-1,2:end-1)=hand.Map(2,2:end-1,2:end-1)*2*hand.TopStim;
% T(end,2:end-1,2:end-1)=hand.Map(end-1,2:end-1,2:end-1)*2*hand.BotStim;
%% Initialise
if hand.TopStim==hand.BotStim || hand.Blocked==1 || sum(sum(hand.Map(2,2:end-1,2:end-1)))==0
    T=single(hand.Map*hand.BotStim);
elseif hand.Check_TauMode ==3
    %T(2:end-1,2:end-1,2:end-1)=0*hand.Map(2:end-1,2:end-1,2:end-1);
    step=(hand.BotStim-hand.TopStim)/a;
    V= (hand.TopStim+step/2:step:hand.BotStim);
    T(2:end-1,2:end-1,2:end-1)=single(ndgrid(V',1:b,1:c).*hand.Map(2:end-1,2:end-1,2:end-1));
else
    step=(hand.BotStim-hand.TopStim)/a;
    V= (hand.TopStim+step/2:step:hand.BotStim);
    T(2:end-1,2:end-1,2:end-1)=single(ndgrid(V',1:b,1:c).*hand.Map(2:end-1,2:end-1,2:end-1));
end

hand.T1=double(T(hand.Cheq1.P));
hand.T2=double(T(hand.Cheq2.P));
if hand.Check_VaryD ==0
    if hand.Check_TauMode ~=3
        hand.T1(end+2)=2*hand.BotStim;
        hand.T2(end+2)=2*hand.BotStim;
    else
        hand.T1(end+2)=0;
        hand.T2(end+2)=0;
    end
else
    hand.T1(end+2)=hand.BotStim;
    hand.T2(end+2)=hand.BotStim;
end
if hand.InLineMode==0
    try
        [hand]=InitiatePlot1(hand);
    catch
         disp('Unable to plot')
    end
    try
        [hand]=InitiatePlot2(hand,T(:,:,2));
    catch
         disp('Unable to plot')
    end
end
%[hand]=MBytesRecord(hand,whos,'Prep3'); %Memory

function [hand]=Preparation3imp(hand)
%% Alternative third preparation step if impedance is called
% hand.TopStim=0;
% hand.BotStim=1;
% T=complex(hand.Tconv);
% hand.Tconv=1;
% T(1,2:end-1,2:end-1)=2*hand.TopStim;
% T(end,2:end-1,2:end-1)=2*hand.BotStim;

[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));

if get(hand.Check_TauMode,'Value')==4
    hand.Tconv=complex((hand.Tconv));
elseif get(hand.Check_TauMode,'Value')==6
    hand.Tconv=complex(ones(size(hand.Map)));
    %     step=(hand.BotStim-hand.TopStim)/a;
    %     V= (hand.TopStim+step/2:step:hand.BotStim);
    %     hand.Tconv(2:end-1,2:end-1,2:end-1)=hand.Tconv(2:end-1,2:end-1,2:end-1)+1i*single(ndgrid(V',1:b,1:c).*hand.Map(2:end-1,2:end-1,2:end-1));
end
hand.T1=double(hand.Tconv(hand.Cheq1.P));
hand.T2=double(hand.Tconv(hand.Cheq2.P));
if get(hand.Check_VaryD,'value')==0
    hand.T1(end+2)=2*hand.BotStim;
    hand.T2(end+2)=2*hand.BotStim;
else
    hand.T1(end+2)=hand.BotStim;
    hand.T2(end+2)=hand.BotStim;
end
hand.Tconv=hand.Tconv(:,:,2);
hand.Area_top=complex(sum(sum(hand.Map(2,:,:))));
hand.Area_bot=complex(sum(sum(hand.Map(end-1,:,:))));
hand.w=complex(2-(pi)/(a*1.5));
if get(hand.Check_TauMode,'Value')==4
    if hand.InLineMode==0
        try
            [hand]=InitiatePlot2(hand,(hand.Tconv));
        catch
             disp('Unable to plot')
        end
    end
    if hand.PercFlag==1 && hand.Blocked==0
        if hand.y(hand.freqNo)>0.5
            hand.w=hand.w*0.90^(hand.y(hand.freqNo)-0.5);%-3
        end
    else
        if hand.y(hand.freqNo)>-1
            hand.w=hand.w*0.93^(hand.y(hand.freqNo)+1);%-3
        end
    end
    if hand.Aniso==1
        hand.w=0.95*hand.w;
    end
    hand.omw=complex(1-hand.w);
    hand.NN_aV.w1=complex(hand.w./ (complex((1i*hand.freq*hand.delta_x^2)/...
        hand.D+complex(double(hand.NN_a(hand.Cheq1.P)))) ) );
    hand.NN_aV.w2=complex(hand.w./ (complex((1i*hand.freq*hand.delta_x^2)/...
        hand.D+complex(double(hand.NN_a(hand.Cheq2.P)))) ) );
elseif get(hand.Check_TauMode,'Value')==6
    %     hand.w=1;
    DimNei=complex(2*length(size(hand.Net_Or)));
    hand.omw=complex(1-hand.w);
    hand.NN_aV.w1=complex(hand.w.*...
        complex(1./(double(hand.NN_a(hand.Cheq1.P)))).*...
        (1+(2i*(complex(double(hand.NN_a(hand.Cheq1.P)))-DimNei).*hand.freq)./...
        complex(double(hand.NN_a(hand.Cheq1.P)))));
    hand.NN_aV.w2=complex(hand.w.*...
        complex(1./(double(hand.NN_a(hand.Cheq2.P)))).*...
        (1+(2i*(complex(double(hand.NN_a(hand.Cheq2.P)))-DimNei).*hand.freq)./...
        complex(double(hand.NN_a(hand.Cheq2.P)))));
    %% Enforce mirror
    
    %% Correct the bottom
    hand.NN_aV.w1(hand.T1Bot)=complex(hand.w.*...
        complex(1./(double(hand.NN_a(hand.Cheq1.P(hand.T1Bot))))).*...
        (1+(2i*(complex(double(hand.NN_a(hand.Cheq1.P(hand.T1Bot))))-DimNei-1).*hand.freq)./...
        complex(double(hand.NN_a(hand.Cheq1.P(hand.T1Bot))-1))));
    hand.NN_aV.w2(hand.T2Bot)=complex(hand.w.*...
        complex(1./(double(hand.NN_a(hand.Cheq2.P(hand.T2Bot))))).*...
        (1+(2i*(complex(double(hand.NN_a(hand.Cheq2.P(hand.T2Bot))))-DimNei-1).*hand.freq)./...
        complex(double(hand.NN_a(hand.Cheq2.P(hand.T2Bot))-1))));
    %% Correct the top
    if hand.Blocked==0
        hand.NN_aV.w1(hand.T1Top)=complex(hand.w.*...
            complex(1./(double(hand.NN_a(hand.Cheq1.P(hand.T1Top))-2))).*...
            (1+(2i*(complex(double(hand.NN_a(hand.Cheq1.P(hand.T1Top))))-DimNei-1).*hand.freq)./...
            complex(double(hand.NN_a(hand.Cheq1.P(hand.T1Top))-2))));
        hand.NN_aV.w2(hand.T2Top)=complex(hand.w.*...
            complex(1./(double(hand.NN_a(hand.Cheq2.P(hand.T2Top))-2))).*...
            (1+(2i*(complex(double(hand.NN_a(hand.Cheq2.P(hand.T2Top))))-DimNei-1).*hand.freq)./...
            complex(double(hand.NN_a(hand.Cheq2.P(hand.T2Top))-2))));
        %     else %% If it is blocked, it's already corrrect!!
        %         hand.NN_aV.w1(hand.T1Top)=complex(hand.w.*...
        %             complex(1./(double(hand.NN_a(hand.Cheq1.P(hand.T1Top))))).*...
        %             (1+(2i*(complex(double(hand.NN_a(hand.Cheq1.P(hand.T1Top))))-DimNei).*hand.freq)./...
        %             complex(double(hand.NN_a(hand.Cheq1.P(hand.T1Top))))));
        %         hand.NN_aV.w2(hand.T2Top)=complex(hand.w.*...
        %             complex(1./(double(hand.NN_a(hand.Cheq2.P(hand.T2Top))))).*...
        %             (1+(2i*(complex(double(hand.NN_a(hand.Cheq2.P(hand.T2Top))))-DimNei).*hand.freq)./...
        %             complex(double(hand.NN_a(hand.Cheq2.P(hand.T2Top))))));
    end
    if hand.InLineMode==0
        try
            [hand]=InitiatePlot1(hand);
        end
        try
            [hand]=InitiatePlot2(hand,hand.Tconv);
        end
    end
end

%[hand]=MBytesRecord(hand,whos,'Prep3Imp'); %Memory




function [hand]=Preparation3bis(hand)
%% Third preparation step for Tau calculation where the volume is initialised as linear
[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));

T=single(hand.Map);

% Boundary conditions
hand.TopStim=0;  % Potential of top electrode (as Ref)
hand.BotStim=1;   % Potential of bottom electrode

if hand.TopStim==hand.BotStim || hand.Blocked==1 || sum(sum(hand.Map(2,2:end-1,2:end-1)))==0
    T=single(hand.Map*hand.BotStim);
else
    step=(hand.BotStim-hand.TopStim)/a;
    V= (hand.TopStim+step/2:step:hand.BotStim);
    T(2:end-1,2:end-1,2:end-1)=single(ndgrid(V',1:b,1:c).*hand.Map(2:end-1,2:end-1,2:end-1));
end

hand.T1=double(T(hand.Cheq1.P));
hand.T2=double(T(hand.Cheq2.P));

hand.T1(end+2)=hand.BotStim;
hand.T2(end+2)=hand.BotStim;


if hand.InLineMode==0
    T=T(:,:,2);
    try
        [hand]=InitiatePlot1(hand);
    end
    try
        [hand]=InitiatePlot2(hand,T,hand.Map);
    end
end
%[hand]=MBytesRecord(hand,whos,'Prep3bis'); %Memory

function [hand]=Preparation3imp_sym(hand)
%% Alternative third preparation step if impedance is called

% hand.Tconv=single(hand.Map);
hand.Tconv=complex((hand.Tconv));
hand.T1=double(hand.Tconv(hand.Cheq1.P));
hand.T2=double(hand.Tconv(hand.Cheq2.P));
hand.T1(end+2)=hand.BotStim;
hand.T2(end+2)=hand.BotStim;

if hand.InLineMode==0
    hand.Tconv=hand.Tconv(:,:,2);
    try
        [hand]=InitiatePlot2(hand,(hand.Tconv),hand.Map);
    end
end

hand.w = complex(0.85);   % fix omega

% hand.w=complex(2-(pi)/max(size(hand.Map))*0.5);

% hand.w = complex(1.3);   % fix omega

hand.omw=complex(1-hand.w);

hand.NN_aV.sp1 = complex((complex(double(hand.NN_tot(hand.Cheq1.P))-double(hand.NN_a(hand.Cheq1.P)))*1i*2*pi*hand.freq*hand.C_dl*hand.delta_x^2)/(hand.kappa*hand.delta_x)) ;
hand.NN_aV.sp2 = complex((complex(double(hand.NN_tot(hand.Cheq2.P))-double(hand.NN_a(hand.Cheq2.P)))*1i*2*pi*hand.freq*hand.C_dl*hand.delta_x^2)/(hand.kappa*hand.delta_x)) ;

hand.NN_aV.w1=complex(hand.w./ (complex(double(hand.NN_aV.sp1)+complex(double(hand.NN(hand.Cheq1.P))))));
hand.NN_aV.w2=complex(hand.w./ (complex(double(hand.NN_aV.sp2)+complex(double(hand.NN(hand.Cheq2.P))))));

%[hand]=MBytesRecord(hand,whos,'Prep3Imp_SymCell'); %Memory

function [hand]=Iterate_mat(hObject, hand)
%% Core iteration function for tortuosity factor calculation
hand.iter=0;
% Error=1;
% checkNo=1;
hand.whileFlag=1;
[~,~,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
%[hand]=MBytesRecord(hand,whos,'Iterate'); %Memory
if c>1 % if the volume is 3D
    if hand.Check_VaryD ==0
        if hand.Aniso==0 % if the volume has isotropic voxels
            if hand.Check_TauMode ~=3
                %% Parallel initialisation
                %                 p=4;
                %                 for n=1:p
                %                     b1(n,1)=floor(1+(n-1)*length(hand.T1(1:end-2))/p);
                %                     b1(n,2)=floor((n)*length(hand.T1(1:end-2))/p);
                %                     b2(n,1)=floor(1+(n-1)*length(hand.T2(1:end-2))/p);
                %                     b2(n,2)=floor((n)*length(hand.T2(1:end-2))/p);
                %                 end
                
                while hand.whileFlag>0 && hand.iter<hand.iter_max
                    hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                        hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp)+...
                        hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp)+...
                        hand.T2(hand.Cheq1.P_Zm)+hand.T2(hand.Cheq1.P_Zp));
                    hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                        hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp)+...
                        hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp)+...
                        hand.T1(hand.Cheq2.P_Zm)+hand.T1(hand.Cheq2.P_Zp));
                    %                     parfor n=1:p
                    %                         hand.T1(b1(n,1):b1(n,2))=hand.omw*hand.T1(b1(n,1):b1(n,2))+hand.NN_aV.w1(b1(n,1):b1(n,2)).*(...
                    %                             hand.T2(hand.Cheq1.P_Xm(b1(n,1):b1(n,2)))+hand.T2(hand.Cheq1.P_Xp(b1(n,1):b1(n,2)))+...
                    %                             hand.T2(hand.Cheq1.P_Ym(b1(n,1):b1(n,2)))+hand.T2(hand.Cheq1.P_Yp(b1(n,1):b1(n,2)))+...
                    %                             hand.T2(hand.Cheq1.P_Zm(b1(n,1):b1(n,2)))+hand.T2(hand.Cheq1.P_Zp(b1(n,1):b1(n,2))));
                    %                     end
                    %                     parfor n=1:p
                    %                         hand.T2(b2(n,1):b2(n,2))=hand.omw*hand.T2(b2(n,1):b2(n,2))+hand.NN_aV.w2(b2(n,1):b2(n,2)).*(...
                    %                             hand.T1(hand.Cheq2.P_Xm(b2(n,1):b2(n,2)))+hand.T1(hand.Cheq2.P_Xp(b2(n,1):b2(n,2)))+...
                    %                             hand.T1(hand.Cheq2.P_Ym(b2(n,1):b2(n,2)))+hand.T1(hand.Cheq2.P_Yp(b2(n,1):b2(n,2)))+...
                    %                             hand.T1(hand.Cheq2.P_Zm(b2(n,1):b2(n,2)))+hand.T1(hand.Cheq2.P_Zp(b2(n,1):b2(n,2))));
                    %                     end
                    
                    
                    
                    
                    hand.iter=hand.iter+1;
                    if rem(hand.iter,hand.check_f)==0
                        [hand]=Checks(hObject, hand);
                    end
                    
                end
            else %P:P 3D
                while hand.whileFlag>0 && hand.iter<hand.iter_max
                    hand.T1(hand.PPmid1)=hand.omw*hand.T1(hand.PPmid1)+hand.NN_aV.w1(hand.PPmid1).*(...
                        hand.T2(hand.Cheq1.P_Xm(hand.PPmid1))+hand.T2(hand.Cheq1.P_Xp(hand.PPmid1))+...
                        hand.T2(hand.Cheq1.P_Ym(hand.PPmid1))+hand.T2(hand.Cheq1.P_Yp(hand.PPmid1))+...
                        hand.T2(hand.Cheq1.P_Zm(hand.PPmid1))+hand.T2(hand.Cheq1.P_Zp(hand.PPmid1)));
                    hand.T1(hand.PPtop1(2:end))=hand.omw*hand.T1(hand.PPtop1(2:end))+hand.NN_aV.w1(hand.PPtop1(2:end)).*(...
                        hand.T2(hand.Cheq1.P_Xm(hand.PPtop1(2:end)))-1+hand.T2(hand.Cheq1.P_Xp(hand.PPtop1(2:end)))+...
                        hand.T2(hand.Cheq1.P_Ym(hand.PPtop1(2:end)))+hand.T2(hand.Cheq1.P_Yp(hand.PPtop1(2:end)))+...
                        hand.T2(hand.Cheq1.P_Zm(hand.PPtop1(2:end)))+hand.T2(hand.Cheq1.P_Zp(hand.PPtop1(2:end))));
                    hand.T1(hand.PPbot1)=hand.omw*hand.T1(hand.PPbot1)+hand.NN_aV.w1(hand.PPbot1).*(...
                        hand.T2(hand.Cheq1.P_Xm(hand.PPbot1))+1+hand.T2(hand.Cheq1.P_Xp(hand.PPbot1))+...
                        hand.T2(hand.Cheq1.P_Ym(hand.PPbot1))+hand.T2(hand.Cheq1.P_Yp(hand.PPbot1))+...
                        hand.T2(hand.Cheq1.P_Zm(hand.PPbot1))+hand.T2(hand.Cheq1.P_Zp(hand.PPbot1)));
                    
                    hand.T2(hand.PPmid2)=hand.omw*hand.T2(hand.PPmid2)+hand.NN_aV.w2(hand.PPmid2).*(...
                        hand.T1(hand.Cheq2.P_Xm(hand.PPmid2))+hand.T1(hand.Cheq2.P_Xp(hand.PPmid2))+...
                        hand.T1(hand.Cheq2.P_Ym(hand.PPmid2))+hand.T1(hand.Cheq2.P_Yp(hand.PPmid2))+...
                        hand.T1(hand.Cheq2.P_Zm(hand.PPmid2))+hand.T1(hand.Cheq2.P_Zp(hand.PPmid2)));
                    hand.T2(hand.PPtop2)=hand.omw*hand.T2(hand.PPtop2)+hand.NN_aV.w2(hand.PPtop2).*(...
                        hand.T1(hand.Cheq2.P_Xm(hand.PPtop2))-1+hand.T1(hand.Cheq2.P_Xp(hand.PPtop2))+...
                        hand.T1(hand.Cheq2.P_Ym(hand.PPtop2))+hand.T1(hand.Cheq2.P_Yp(hand.PPtop2))+...
                        hand.T1(hand.Cheq2.P_Zm(hand.PPtop2))+hand.T1(hand.Cheq2.P_Zp(hand.PPtop2)));
                    hand.T2(hand.PPbot2)=hand.omw*hand.T2(hand.PPbot2)+hand.NN_aV.w2(hand.PPbot2).*(...
                        hand.T1(hand.Cheq2.P_Xm(hand.PPbot2))+1+hand.T1(hand.Cheq2.P_Xp(hand.PPbot2))+...
                        hand.T1(hand.Cheq2.P_Ym(hand.PPbot2))+hand.T1(hand.Cheq2.P_Yp(hand.PPbot2))+...
                        hand.T1(hand.Cheq2.P_Zm(hand.PPbot2))+hand.T1(hand.Cheq2.P_Zp(hand.PPbot2)));
                    
                    hand.iter=hand.iter+1;
                    if rem(hand.iter,hand.check_f)==0
                        [hand]=Checks(hObject, hand);
                    end
                    
                end
            end
        elseif hand.Aniso==1 % if the volume has anisotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                    hand.c_X*(hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp))+...
                    hand.c_Y*(hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp))+...
                    hand.c_Z*(hand.T2(hand.Cheq1.P_Zm)+hand.T2(hand.Cheq1.P_Zp)));
                hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                    hand.c_X*(hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp))+...
                    hand.c_Y*(hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp))+...
                    hand.c_Z*(hand.T1(hand.Cheq2.P_Zm)+hand.T1(hand.Cheq2.P_Zp)));
                hand.iter=hand.iter+1;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hObject, hand);
                end
            end
        end
    else % Variable D
        while hand.whileFlag>0 && hand.iter<hand.iter_max
            hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+...
                hand.R.Xm1.*hand.T2(hand.Cheq1.P_Xm)+...
                hand.R.Xp1.*hand.T2(hand.Cheq1.P_Xp)+...
                hand.R.Ym1.*hand.T2(hand.Cheq1.P_Ym)+...
                hand.R.Yp1.*hand.T2(hand.Cheq1.P_Yp)+...
                hand.R.Zm1.*hand.T2(hand.Cheq1.P_Zm)+...
                hand.R.Zp1.*hand.T2(hand.Cheq1.P_Zp);
            hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+...
                hand.R.Xm2.*hand.T1(hand.Cheq2.P_Xm)+...
                hand.R.Xp2.*hand.T1(hand.Cheq2.P_Xp)+...
                hand.R.Ym2.*hand.T1(hand.Cheq2.P_Ym)+...
                hand.R.Yp2.*hand.T1(hand.Cheq2.P_Yp)+...
                hand.R.Zm2.*hand.T1(hand.Cheq2.P_Zm)+...
                hand.R.Zp2.*hand.T1(hand.Cheq2.P_Zp);
            hand.iter=hand.iter+1;
            if rem(hand.iter,hand.check_f)==0
                [hand]=Checks(hObject, hand);
            end
        end
    end
else %% 2D versions
    if hand.Check_VaryD ==0
        if hand.Aniso==0 % if the volume has isotropic voxels
            if hand.Check_TauMode ~=3
                while hand.whileFlag>0 && hand.iter<hand.iter_max
                    hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                        hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp)+...
                        hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp));
                    hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                        hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp)+...
                        hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp));
                    hand.iter=hand.iter+1;
                    if rem(hand.iter,hand.check_f)==0
                        [hand]=Checks(hObject, hand);
                    end
                end
            else %P:P 2D
                while hand.whileFlag>0 && hand.iter<hand.iter_max
                    %                     for i=1:1
                    hand.T1(hand.PPmid1)=hand.omw*hand.T1(hand.PPmid1)+hand.NN_aV.w1(hand.PPmid1).*(...
                        hand.T2(hand.Cheq1.P_Xm(hand.PPmid1))+hand.T2(hand.Cheq1.P_Xp(hand.PPmid1))+...
                        hand.T2(hand.Cheq1.P_Ym(hand.PPmid1))+hand.T2(hand.Cheq1.P_Yp(hand.PPmid1)));
                    hand.T1(hand.PPtop1(2:end))=hand.omw*hand.T1(hand.PPtop1(2:end))+hand.NN_aV.w1(hand.PPtop1(2:end)).*(...
                        hand.T2(hand.Cheq1.P_Xm(hand.PPtop1(2:end)))-1+hand.T2(hand.Cheq1.P_Xp(hand.PPtop1(2:end)))+...
                        hand.T2(hand.Cheq1.P_Ym(hand.PPtop1(2:end)))+hand.T2(hand.Cheq1.P_Yp(hand.PPtop1(2:end))));
                    hand.T1(hand.PPbot1)=hand.omw*hand.T1(hand.PPbot1)+hand.NN_aV.w1(hand.PPbot1).*(...
                        hand.T2(hand.Cheq1.P_Xm(hand.PPbot1))+1+hand.T2(hand.Cheq1.P_Xp(hand.PPbot1))+...
                        hand.T2(hand.Cheq1.P_Ym(hand.PPbot1))+hand.T2(hand.Cheq1.P_Yp(hand.PPbot1)));
                    
                    hand.T2(hand.PPmid2)=hand.omw*hand.T2(hand.PPmid2)+hand.NN_aV.w2(hand.PPmid2).*(...
                        hand.T1(hand.Cheq2.P_Xm(hand.PPmid2))+hand.T1(hand.Cheq2.P_Xp(hand.PPmid2))+...
                        hand.T1(hand.Cheq2.P_Ym(hand.PPmid2))+hand.T1(hand.Cheq2.P_Yp(hand.PPmid2)));
                    hand.T2(hand.PPtop2)=hand.omw*hand.T2(hand.PPtop2)+hand.NN_aV.w2(hand.PPtop2).*(...
                        hand.T1(hand.Cheq2.P_Xm(hand.PPtop2))-1+hand.T1(hand.Cheq2.P_Xp(hand.PPtop2))+...
                        hand.T1(hand.Cheq2.P_Ym(hand.PPtop2))+hand.T1(hand.Cheq2.P_Yp(hand.PPtop2)));
                    hand.T2(hand.PPbot2)=hand.omw*hand.T2(hand.PPbot2)+hand.NN_aV.w2(hand.PPbot2).*(...
                        hand.T1(hand.Cheq2.P_Xm(hand.PPbot2))+1+hand.T1(hand.Cheq2.P_Xp(hand.PPbot2))+...
                        hand.T1(hand.Cheq2.P_Ym(hand.PPbot2))+hand.T1(hand.Cheq2.P_Yp(hand.PPbot2)));
                    %                     end
                    %                     Tmap=nan(size(hand.Map));Tmap(hand.Cheq1.P)=hand.T1(1:end-2);Tmap(hand.Cheq2.P)=hand.T2(1:end-2);Tmap(2:end-1,2:end-1,2);
                    %                     imagesc(Tmap(:,:,2));colormap(spring);colorbar;Tmap(:,:,2)
                    
                    hand.iter=hand.iter+1;
                    if rem(hand.iter,hand.check_f)==0
                        [hand]=Checks(hObject, hand);
                    end
                end
            end
        else % if the volume has anisotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                    hand.c_X*(hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp))+...
                    hand.c_Y*(hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp)));
                hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                    hand.c_X*(hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp))+...
                    hand.c_Y*(hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp)));
                hand.iter=hand.iter+1;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hObject, hand);
                end
            end %iterations
        end
    else % Variable D
        while hand.whileFlag>0 && hand.iter<hand.iter_max
            hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+...
                hand.R.Xm1.*hand.T2(hand.Cheq1.P_Xm)+...
                hand.R.Xp1.*hand.T2(hand.Cheq1.P_Xp)+...
                hand.R.Ym1.*hand.T2(hand.Cheq1.P_Ym)+...
                hand.R.Yp1.*hand.T2(hand.Cheq1.P_Yp);
            hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+...
                hand.R.Xm2.*hand.T1(hand.Cheq2.P_Xm)+...
                hand.R.Xp2.*hand.T1(hand.Cheq2.P_Xp)+...
                hand.R.Ym2.*hand.T1(hand.Cheq2.P_Ym)+...
                hand.R.Yp2.*hand.T1(hand.Cheq2.P_Yp);
            hand.iter=hand.iter+1;
            if rem(hand.iter,hand.check_f)==0
                [hand]=Checks(hObject, hand);
            end
        end
    end
end
%[hand]=MBytesRecord(hand,whos,'Iterate complete'); %Memory
[hand]=AfterIterate(hand);

function [hand]=Iterate_mex(hObject, hand)
%% Core iteration function for tortuosity factor calculation
hand.iter=0;
% Error=1;
% checkNo=1;
hand.whileFlag=1;
[~,~,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
%[hand]=MBytesRecord(hand,whos,'Iterate'); %Memory
if hand.impCheck==0
    % ensure variable type is appropriate
else
end
%[hand]=MBytesRecord(hand,whos,'Iterate'); %Memory
if numel(hand.Map)<2^32
    if hand.impCheck==0 %Not impedance
        if hand.Check_VaryD ==0
            if hand.Aniso==0 % if the volume has isotropic voxels
                if sum(hand.Check_TauMode ==[1 2 4 5])
                    while hand.whileFlag>0 && hand.iter<hand.iter_max
                        [hand.T1,hand.T2]=Mex3DTauIso_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
                            hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                            hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp);
                        hand.iter=hand.iter+hand.check_f;
                        if rem(hand.iter,hand.check_f)==0
                            [hand]=Checks(hObject, hand);
                        end
                    end
                elseif hand.Check_TauMode ==3
                    while hand.whileFlag>0 && hand.iter<hand.iter_max
                        [hand.T1,hand.T2]=Mex3DTauIsoPP_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
                            hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                            hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp,...
                            hand.PPmid1,hand.PPmid2,hand.PPtop1,hand.PPtop2,hand.PPbot1,hand.PPbot2);
                        hand.iter=hand.iter+hand.check_f;
                        if rem(hand.iter,hand.check_f)==0
                            [hand]=Checks(hObject, hand);
                        end
                    end
                elseif hand.Check_TauMode ==6
                    while hand.whileFlag>0 && hand.iter<hand.iter_max
                        [hand.T1,hand.T2]=Mex3DTauEIso_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
                            hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                            hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp);
                        hand.iter=hand.iter+hand.check_f;
                        if rem(hand.iter,hand.check_f)==0
                            [hand]=Checks(hObject, hand);
                        end
                    end
                end
            elseif hand.Aniso==1 % if the volume has anisotropic voxels
                while hand.whileFlag>0 && hand.iter<hand.iter_max
                    [hand.T1,hand.T2]=Mex3DTauAni_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
                        hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                        hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp,...
                        hand.c_X,hand.c_Y,hand.c_Z);
                    hand.iter=hand.iter+hand.check_f;
                    if rem(hand.iter,hand.check_f)==0
                        [hand]=Checks(hObject, hand);
                    end
                end
            end
        else % Variable D
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                [hand.T1,hand.T2]=Mex3DTauVar_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,...
                    hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                    hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp,...
                    hand.R.Xm1,hand.R.Xp1,hand.R.Ym1,hand.R.Yp1,hand.R.Zm1,hand.R.Zp1,...
                    hand.R.Xm2,hand.R.Xp2,hand.R.Ym2,hand.R.Yp2,hand.R.Zm2,hand.R.Zp2);
                hand.iter=hand.iter+hand.check_f;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hObject, hand);
                end
            end
        end
    else %Impedance
        while hand.whileFlag>0 && hand.iter<hand.iter_max
            [hand.T1,hand.T2]=Mex3DTauIso_Imp_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,double(hand.NN_aV.w1),double(hand.NN_aV.w2),...
                hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp);
            hand.iter=hand.iter+hand.check_f;
            if rem(hand.iter,hand.check_f)==0
                [hand]=Checks(hObject, hand);
            end
        end
    end
else
    if hand.impCheck==0 %Not impedance
        if hand.Check_VaryD ==0
            if hand.Aniso==0 % if the volume has isotropic voxels
                while hand.whileFlag>0 && hand.iter<hand.iter_max
                    [hand.T1,hand.T2]=Mex3DTauIso_64_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
                        hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                        hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp);
                    hand.iter=hand.iter+hand.check_f;
                    if rem(hand.iter,hand.check_f)==0
                        [hand]=Checks(hObject, hand);
                    end
                end
            elseif hand.Aniso==1 % if the volume has anisotropic voxels
                while hand.whileFlag>0 && hand.iter<hand.iter_max
                    [hand.T1,hand.T2]=Mex3DTauAni_64_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
                        hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                        hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp,...
                        hand.c_X,hand.c_Y,hand.c_Z);
                    hand.iter=hand.iter+hand.check_f;
                    if rem(hand.iter,hand.check_f)==0
                        [hand]=Checks(hObject, hand);
                    end
                end
            end
        else % Variable D
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                [hand.T1,hand.T2]=Mex3DTauVar_64_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,...
                    hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                    hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp,...
                    hand.R.Xm1,hand.R.Xp1,hand.R.Ym1,hand.R.Yp1,hand.R.Zm1,hand.R.Zp1,...
                    hand.R.Xm2,hand.R.Xp2,hand.R.Ym2,hand.R.Yp2,hand.R.Zm2,hand.R.Zp2);
                hand.iter=hand.iter+hand.check_f;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hObject, hand);
                end
            end
        end
    else %Impedance
        while hand.whileFlag>0 && hand.iter<hand.iter_max
            [hand.T1,hand.T2]=Mex3DTauIso_Imp_64_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,double(hand.NN_aV.w1),double(hand.NN_aV.w2),...
                hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
                hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp);
            hand.iter=hand.iter+hand.check_f;
            if rem(hand.iter,hand.check_f)==0
                [hand]=Checks(hObject, hand);
            end
        end
    end
end
%[hand]=MBytesRecord(hand,whos,'Iterate complete'); %Memory
[hand]=AfterIterate(hand);

function [hand]=Iterate_uint64_mex(hObject, hand)
% %% Core iteration function for tortuosity factor calculation
% hand.iter=0;
% Error=1;
% checkNo=1;
% hand.whileFlag=1;
% [~,~,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
% [hand]=MBytesRecord(hand,whos,'Iterate'); %Memory
% if hand.impCheck==0
%     % ensure variable type is appropriate
% else
% end
% [hand]=MBytesRecord(hand,whos,'Iterate'); %Memory
% if hand.impCheck==0 %Not impedance
%     if get(hand.Check_VaryD,'value')==0
%         if hand.Aniso==0; % if the volume has isotropic voxels
%             while hand.whileFlag>0 && hand.iter<hand.iter_max
%                 [hand.T1,hand.T2]=Mex3DTauIso_uint64_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
%                     hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
%                     hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp);
%                 hand.iter=hand.iter+hand.check_f;
%                 if rem(hand.iter,hand.check_f)==0
%                     [hand]=Checks(hObject, hand);
%                 end
%             end
%         elseif hand.Aniso==1 % if the volume has anisotropic voxels
%             while hand.whileFlag>0 && hand.iter<hand.iter_max
%                 [hand.T1,hand.T2]=Mex3DTauAni_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,hand.NN_aV.w1,hand.NN_aV.w2,...
%                     hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
%                     hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp,...
%                     hand.c_X,hand.c_Y,hand.c_Z);
%                 hand.iter=hand.iter+hand.check_f;
%                 if rem(hand.iter,hand.check_f)==0
%                     [hand]=Checks(hObject, hand);
%                 end
%             end
%         end
%     else % Variable D
%         while hand.whileFlag>0 && hand.iter<hand.iter_max
%             [hand.T1,hand.T2]=Mex3DTauVar_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,...
%                 hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
%                 hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp,...
%                 hand.R.Xm1,hand.R.Xp1,hand.R.Ym1,hand.R.Yp1,hand.R.Zm1,hand.R.Zp1,...
%                 hand.R.Xm2,hand.R.Xp2,hand.R.Ym2,hand.R.Yp2,hand.R.Zm2,hand.R.Zp2);
%             hand.iter=hand.iter+hand.check_f;
%             if rem(hand.iter,hand.check_f)==0
%                 [hand]=Checks(hObject, hand);
%             end
%         end
%     end
% else %Impedance
%     while hand.whileFlag>0 && hand.iter<hand.iter_max
%         [hand.T1,hand.T2]=Mex3DTauIso_Imp_mex(c,hand.check_f,hand.T1,hand.T2,hand.omw,double(hand.NN_aV.w1),double(hand.NN_aV.w2),...
%             hand.Cheq1.P_Xm,hand.Cheq1.P_Xp,hand.Cheq1.P_Ym,hand.Cheq1.P_Yp,hand.Cheq1.P_Zm,hand.Cheq1.P_Zp,...
%             hand.Cheq2.P_Xm,hand.Cheq2.P_Xp,hand.Cheq2.P_Ym,hand.Cheq2.P_Yp,hand.Cheq2.P_Zm,hand.Cheq2.P_Zp);
%         hand.iter=hand.iter+hand.check_f;
%         if rem(hand.iter,hand.check_f)==0
%             [hand]=Checks(hObject, hand);
%         end
%     end
% end
%
% [hand]=MBytesRecord(hand,whos,'Iterate complete'); %Memory
% [hand]=AfterIterate(hand);

function [hand]=Iterate_sym_mat(hObject,hand)
%% Core iteration function for tortuosity factor calculation
hand.iter=0;
Error=1;
checkNo=1;
hand.whileFlag=1;
[~,~,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
%[hand]=MBytesRecord(hand,whos,'Iterate SymCell'); %Memory
if c>1 % if the volume is 3D
    if hand.Check_VaryD ==0
        if hand.Aniso==0 % if the volume has isotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                    hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp)+...
                    hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp)+...
                    hand.T2(hand.Cheq1.P_Zm)+hand.T2(hand.Cheq1.P_Zp)+...
                    hand.NN_aV.sp1.*(hand.T2(hand.Cheq1.Pb)));
                hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                    hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp)+...
                    hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp)+...
                    hand.T1(hand.Cheq2.P_Zm)+hand.T1(hand.Cheq2.P_Zp)+...
                    hand.NN_aV.sp2.*(hand.T1(hand.Cheq2.Pb)));
                hand.iter=hand.iter+1;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hObject, hand);
                end
            end
        end
    end
else %% 2D version
    if hand.Check_VaryD ==0
        if hand.Aniso==0 % if the volume has isotropic voxels
            while hand.whileFlag>0 && hand.iter<hand.iter_max
                hand.T1(1:end-2)=hand.omw*hand.T1(1:end-2)+hand.NN_aV.w1.*(...
                    hand.T2(hand.Cheq1.P_Xm)+hand.T2(hand.Cheq1.P_Xp)+...
                    hand.T2(hand.Cheq1.P_Ym)+hand.T2(hand.Cheq1.P_Yp)+...
                    hand.NN_aV.sp1.*(hand.T2(hand.Cheq1.Pb)));
                hand.T2(1:end-2)=hand.omw*hand.T2(1:end-2)+hand.NN_aV.w2.*(...
                    hand.T1(hand.Cheq2.P_Xm)+hand.T1(hand.Cheq2.P_Xp)+...
                    hand.T1(hand.Cheq2.P_Ym)+hand.T1(hand.Cheq2.P_Yp)+...
                    hand.NN_aV.sp2.*(hand.T1(hand.Cheq2.Pb)));
                hand.iter=hand.iter+1;
                if rem(hand.iter,hand.check_f)==0
                    [hand]=Checks(hObject, hand);
                end
            end
        end
    end
end
%[hand]=MBytesRecord(hand,whos,'Iterate 4 SymCell complete'); %Memory
[hand]=AfterIterate(hand);

function [hand]=Checks(hObject, hand)
hand.conDwell=5;
hand.conTol=0.005;
[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));

%% Function to check on the convergence status of the volume
checkNo=hand.iter/hand.check_f;
% Calculate_Button the apparent tortuosity factors at the top and bottom faces
hand=XFlux(hand);


if sum(hand.Check_TauMode ==[1 2 4 5])
    if hand.Check_TauMode ~=5
        if  hand.Check_VaryD ==0
            sumT_top=sum(hand.T1(hand.T1Top))+sum(hand.T2(hand.T2Top));
            sumT_bot=sum(hand.T1(hand.T1Bot))+sum(hand.T2(hand.T2Bot));
            
            hand.DeffTop(checkNo)=abs(hand.Area_top*hand.TopStim-sumT_top)*2*hand.D*(a/(b*c));
            hand.DeffBot(checkNo)=abs(hand.Area_bot*hand.BotStim-sumT_bot)*2*hand.D*(a/(b*c));
            
            deltaTB=abs((abs(hand.DeffTop(checkNo))-abs(hand.DeffBot(checkNo)))/mean([hand.DeffTop(checkNo),hand.DeffBot(checkNo)]));
            
            hand.TauFacTop(checkNo)=hand.D*hand.VolFrac./hand.DeffTop(checkNo);
            hand.TauFacBot(checkNo)=hand.D*hand.VolFrac./hand.DeffBot(checkNo);
        else
            F_top=sum([hand.T1(hand.T1Top);hand.T2(hand.T2Top)]./hand.Dmap([hand.Cheq1.P(hand.T1Top);hand.Cheq2.P(hand.T2Top)]))*(hand.L_Y*hand.L_Z/(0.5*hand.L_X));
            F_bot=sum((hand.BotStim-[hand.T1(hand.T1Bot);hand.T2(hand.T2Bot)])./hand.Dmap([hand.Cheq1.P(hand.T1Bot);hand.Cheq2.P(hand.T2Bot)]))*(hand.L_Y*hand.L_Z/(0.5*hand.L_X));
            deltaTB=abs((abs(F_top)-abs(F_bot))/mean([F_top,F_bot]));
            
            hand.DeffTop(checkNo)=abs(F_top/((hand.L_Y*b*hand.L_Z*c)/(hand.L_X*a)));
            hand.DeffBot(checkNo)=abs(F_bot/((hand.L_Y*b*hand.L_Z*c)/(hand.L_X*a)));
            
            hand.TauFacTop(checkNo)=hand.Qcv/abs(F_top);
            hand.TauFacBot(checkNo)=hand.Qcv/abs(F_bot);
        end
        % Check for convergence of tau at the two faces
    end
    if hand.impCheck~=1
        %     hand.SumT(checkNo)=sum(sum(sum(T(2:end-1,:,:))));
        if checkNo<3
            hand.whileFlag=hand.conDwell;
        else
            %         DsumT=abs((hand.SumT(checkNo)-hand.SumT(checkNo-1))/hand.SumT(checkNo-1))
            if      abs((hand.DeffTop(checkNo)-hand.DeffTop(checkNo-2))/hand.DeffTop(checkNo))<hand.conTol
                if      abs(hand.DeffTop(checkNo)-hand.DeffTop(checkNo-1))/hand.DeffTop(checkNo)<hand.conTol &&...
                        abs(hand.DeffTop(checkNo)-hand.DeffBot(checkNo))<hand.conTol &&...
                        deltaTB<0.04 &&...
                        (mean([hand.TauFacTop(checkNo),hand.TauFacBot(checkNo)])-mean([hand.TauFacTop(checkNo-1),hand.TauFacBot(checkNo-1)]))/hand.TauFacTop(checkNo)<hand.conTol
                    hand.whileFlag=hand.whileFlag-1;
                    %[hand]=MBytesRecord(hand,whos,'Checks start'); %Memory
                else
                    if hand.whileFlag~=hand.conDwell
                        hand.whileFlag=hand.conDwell;
                    end
                end
            end
            %%%%%%% Increase relaxation factor progressively?
            %         if      abs(hand.TauFacTop(checkNo)-hand.TauFacTop(checkNo-1))<10*hand.conTol &&...
            %                 abs(hand.TauFacTop(checkNo)-hand.TauFacBot(checkNo))<10*hand.conTol
            %             if hand.w==2-(pi)/max(size(T)*1.50);
            %                 %         [hand,~]=Preparation3(hand);
            %                 hand.w=2-(pi)/max(size(T)*1.55);
            %                 hand.omw=double(1-hand.w);
            %                 hand.NN_aV.w1=double(hand.w./double(hand.NN_a(hand.Cheq1.P)));
            %                 hand.NN_aV.w2=double(hand.w./double(hand.NN_a(hand.Cheq2.P)));
            %                 hand.w_increaseFlag=1;
            %             end
            %         end
        end
    else % Impedance mode
        if hand.Check_TauMode ==4
            TauRat=(hand.L_Y*hand.L_Z)/hand.L_X;
            hand.ImpedanceTop(hand.freqNo,checkNo)=-(0.5/(TauRat*(hand.Area_top*hand.TopStim-sumT_top)/(1/hand.D)));
            hand.ImpedanceBot(hand.freqNo,checkNo)= (0.5/(TauRat*(hand.Area_bot*hand.BotStim-sumT_bot)/(1/hand.D)));
            if checkNo<6
                hand.whileFlag=hand.conDwell-1;
            else
                hand.conTol=0.01;
                %         [mean(abs([hand.T1(:); hand.T2(:)])),...
                %             abs(real(hand.ImpedanceBot(hand.freqNo,checkNo))/real(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1),...
                %             abs(imag(hand.ImpedanceBot(hand.freqNo,checkNo))/imag(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1),...
                %             abs(abs(hand.ImpedanceBot(hand.freqNo,checkNo))/abs(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1),...
                %             abs(angle(hand.ImpedanceBot(hand.freqNo,checkNo))/angle(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1)];
                if      abs(real(hand.ImpedanceBot(hand.freqNo,checkNo))/real(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1)<hand.conTol &&...
                        abs(imag(hand.ImpedanceBot(hand.freqNo,checkNo))/imag(hand.ImpedanceBot(hand.freqNo,checkNo-5))-1)<hand.conTol% &&...
                    %                abs(real(hand.ImpedanceTop(hand.freqNo,checkNo))/real(hand.ImpedanceTop(hand.freqNo,checkNo-6))-1)<hand.conTol &&...
                    %                abs(imag(hand.ImpedanceTop(hand.freqNo,checkNo))/imag(hand.ImpedanceTop(hand.freqNo,checkNo-6))-1)<hand.conTol
                    hand.ImpedanceBotConv(hand.freqNo)=hand.ImpedanceBot(hand.freqNo,checkNo);
                    hand.ImpedanceTopConv(hand.freqNo)=hand.ImpedanceTop(hand.freqNo,checkNo);
                    hand.whileFlag=hand.whileFlag-1;
                    
                else
                    if hand.whileFlag~=hand.conDwell
                        hand.whileFlag=hand.conDwell;
                    end
                end
                if hand.iter>=hand.iter_max-hand.check_f
                    hand.ImpedanceBotConv(hand.freqNo)=hand.ImpedanceBot(hand.freqNo,checkNo);
                    hand.ImpedanceTopConv(hand.freqNo)=hand.ImpedanceTop(hand.freqNo,checkNo);
                    %             Error=1;
                    hand.ConvErr=hand.ConvErr+1;
                    hand.iter_max=hand.iter_max*1.5;
                    disp('Simulation not converged')
                    return
                end
            end
            %%%%%%     if unstable
            %         hand.w=hand.w*0.99;
            %         hand.omw=complex(1-hand.w);
            %
            %         hand.NN_aV.w1=complex(hand.w./ (complex((1i*hand.freq*hand.delta_x^2)/...
            %             hand.D+complex(double(hand.NN_a(hand.Cheq1.P)))) ) );
            %         hand.NN_aV.w2=complex(hand.w./ (complex((1i*hand.freq*hand.delta_x^2)/...
            %             hand.D+complex(double(hand.NN_a(hand.Cheq2.P)))) ) );
            %     end
        elseif hand.Check_TauMode ==5
            hand.cur = -(hand.kappa*hand.delta_x).*...
                (sum((hand.T1(hand.Cheq1.Sep).*(hand.Cheq1.P_Xp(hand.Cheq1.Sep)<=...
                length(hand.Cheq2.P))-hand.T2(hand.Cheq1.P_Xp(hand.Cheq1.Sep))))...
                +sum((hand.T2(hand.Cheq2.Sep).*(hand.Cheq2.P_Xp(hand.Cheq2.Sep)<=...
                length(hand.Cheq1.P))-hand.T1(hand.Cheq2.P_Xp(hand.Cheq2.Sep)))))/...
                (b*c*hand.delta_x^2);  % A/m²
            hand.Impedance(hand.freqNo,checkNo)= (hand.BotStim-hand.TopStim)./(hand.cur); % Ohm.m²
            if checkNo<6
                hand.whileFlag=hand.conDwell-1;
            else
                hand.conTol=0.01;
                if      abs(real(hand.Impedance(hand.freqNo,checkNo))/real(hand.Impedance(hand.freqNo,checkNo-5))-1)<hand.conTol &&...
                        abs(imag(hand.Impedance(hand.freqNo,checkNo))/imag(hand.Impedance(hand.freqNo,checkNo-5))-1)<hand.conTol% &&...
                    hand.ImpedanceConv(hand.freqNo)=hand.Impedance(hand.freqNo,checkNo);
                    hand.whileFlag=hand.whileFlag-1;
                else
                    if hand.whileFlag~=hand.conDwell
                        hand.whileFlag=hand.conDwell;
                    end
                end
                
                if hand.iter>=hand.iter_max-hand.check_f
                    hand.ImpedanceConv(hand.freqNo)=hand.Impedance(hand.freqNo,checkNo);
                    Error=1;
                    hand.ConvErr=hand.ConvErr+1;
                    hand.iter_max=hand.iter_max*1.5;
                    disp('Simulation not converged') %eSCM
                    return
                end
            end
        end
    end
    % Error=1;
    % hand.TauFacTop(checkNo)=TauRat*hand.VolFrac*hand.Qcv/(2*sum(sum(T(2,:,:)))); %Should be Perhand.VolFrac(counter)
    % hand.TauFacBot(checkNo)=TauRat*hand.VolFrac*hand.Qcv/(-2*(sum(sum( T(end-1,:,:)))-hand.Area_bot));
    
    %%%% Damping with Jacobi
    % if rem(hand.iter+1,hand.check_f*hand.JacobiRat)==0
    %     Tnew=T;
    %     Tnew(hand.Cheq1.P)=(hand.NN_aV.w1/hand.w).*(...
    %         hand.c_X*(T(hand.Cheq1.P_Xm)+T(hand.Cheq1.P_Xp))+...
    %         hand.c_Y*(T(hand.Cheq1.P_Ym)+T(hand.Cheq1.P_Yp))+...
    %         hand.c_Z*(T(hand.Cheq1.P_Zm)+T(hand.Cheq1.P_Zp)));
    %
    %     Tnew(hand.Cheq2.P)=(hand.NN_aV.w2/hand.w).*(...
    %         hand.c_X*(T(hand.Cheq2.P_Xm)+T(hand.Cheq2.P_Xp))+...
    %         hand.c_Y*(T(hand.Cheq2.P_Ym)+T(hand.Cheq2.P_Yp))+...
    %         hand.c_Z*(T(hand.Cheq2.P_Zm)+T(hand.Cheq2.P_Zp)));
    %     Error(checkNo)=max(max(max(T-Tnew)));
    %     T=Tnew;
    % else
    %     Error=nan;
    % end
    
elseif hand.Check_TauMode ==3%P:P check
    A_pore=length([hand.PPtop1',hand.PPtop2'])*hand.L_Y*hand.L_Z;
    L_pore=hand.L_X;
    
    deltaT=-mean(diff(...
        [[hand.T1(hand.PPtop1)',hand.T2(hand.PPtop2)']+1;...
        [hand.T2(hand.PPbot2)',hand.T1(hand.PPbot1)']]));
    Q_pore=A_pore*deltaT/L_pore;
    hand.DeffTop(checkNo)=Q_pore/hand.Qcv;
    hand.DeffBot(checkNo)=Q_pore/hand.Qcv;
    hand.TauFacTop(checkNo)=hand.D*hand.VolFrac./hand.DeffTop(checkNo);
    hand.TauFacBot(checkNo)=hand.D*hand.VolFrac./hand.DeffTop(checkNo);
    if ~isnan(deltaT)
        if checkNo>10
            if (abs(hand.DeffTop(checkNo)-hand.DeffTop(checkNo-7))/hand.DeffTop(checkNo))<hand.conTol/80
                if max(hand.XFlux_norm)-min(hand.XFlux_norm)<0.03
                    hand.whileFlag=0;
                end
            end
            if hand.iter>=hand.iter_max-hand.check_f
                hand.iter_max=hand.iter_max*1.5;
            end
            if sum(hand.TauFacBot(end-5:end)<0)==6 % if something is failing so it's negative
                hand.whileFlag=0;
                hand.DeffTop(checkNo)=0;
                hand.DeffBot(checkNo)=0;
                hand.TauFacTop(checkNo)=inf;
                hand.TauFacTop(checkNo)=inf;
            end
        end
    else % If the percolation checks failed
        hand.whileFlag=0;
        hand.DeffTop(checkNo)=0;
        hand.DeffBot(checkNo)=0;
        hand.TauFacTop(checkNo)=inf;
        hand.TauFacTop(checkNo)=inf;
    end
elseif hand.Check_TauMode ==6 %tau_e
    VF=mean(hand.Map(2:end-1,2:end-1,2:end-1),[1 2 3]);
    VolCorrect=VF*b*c/(2*a);
    hand.DeffTop(checkNo)=(hand.Area_bot-(sum(hand.T1(hand.T1Bot))+sum(hand.T2(hand.T2Bot))));
    hand.TauFacTop(checkNo)=3*VolCorrect/hand.DeffTop(checkNo);
    hand.DeffTop(checkNo)=VF/real(hand.TauFacTop(checkNo));
    %     if checkNo==1
    %         hand.tauefig=figure;
    %     elseif  checkNo>1
    %         hand.Tconv=zeros(size(hand.Map));
    %         hand.T1(end-1)=nan;
    %         hand.T2(end-1)=nan;
    %         hand.Tconv(hand.Cheq1.P)=single(hand.T1(1:end-2));
    %         hand.Tconv(hand.Cheq2.P)=single(hand.T2(1:end-2));
    %         hand.Tconv(~isfinite(hand.Tconv))=0;
    %         figure(hand.tauefig);
    %         subplot(1,2,1);imagesc(real(hand.Tconv(2:end-1,2:end-1,2))); axis image
    %         subplot(1,2,2);imagesc(imag(hand.Tconv(2:end-1,2:end-1,2))); axis image
    %     end
    if checkNo>10
        if abs(real((hand.DeffTop(checkNo)-hand.DeffTop(checkNo-7))/hand.DeffTop(checkNo)))<hand.conTol/1000 %1000
            hand.whileFlag=0;
        end
        if hand.iter>=hand.iter_max-hand.check_f
            hand.iter_max=hand.iter_max*1.5;
        end
    end
end

%% Update figure
if hand.InLineMode==0
    disp('NOT INLINE MODE');
    %if checkNo>2
        %if hand.impCheck==0
           % if get(hand.Check_RVA, 'Value')==0
               % try
               %     set(hand.ha.Flux(1),'XData',real(hand.XFlux_norm))
               % catch
               %                         disp('Unable to plot 1')
               % end
                %try
                    %set(hand.ha.Tau(1),...
                      %  'XData',hand.check_f*(1:length(hand.DeffTop)),...
                     %   'YData',real(hand.DeffTop));
                    %if sum(get(hand.Check_TauMode,'Value')==[1 2 4 5])
                   %     set(hand.ha.Tau(2),...
                  %          'XData',hand.check_f*(1:length(hand.DeffBot)),...
                 %           'YData',real(hand.DeffBot));
                %    end
               %     if rem(hand.iter/hand.check_f,5)==0
              %          set(hand.axesHandles(7),'XLim',[0 round(hand.iter+hand.check_f*5,2,'significant')]);
             %           set(hand.axesHandles(7),'XTick',...
            %                round((hand.iter+hand.check_f*5)*[0:0.3:0.9],2,'significant'));
           %         end
          %      catch
         %           disp('Unable to plot 2')
        %        end
       %         drawnow
      %      end
     %   else
    %    end
   % end
end

%% Stop?
if hand.InLineMode==0
    disp('NOT INLINE MODE');
   % if get(hand.Calculate_Button,'Value')==0
   %     hand.whileFlag=0;
   %     set(hand.Calculate_Button,...
   %         'String','Stop?',...
   %         'Value',1)
   %     drawnow
   % end
end
%guidata(hObject, hand);
%%%% Switch to double precision if not converging
% if checkNo>10 && mean(class(T)=='single')==1 &&...
%         abs(sum(diff(hand.TauFacTop(end-5:end)-hand.TauFacBot(end-5:end))))<0.001
% %     disp('Switching to double')
%     T=double(T);
% end

function hand=XFlux(hand)
C=zeros(size(hand.Map),'single');
em1=hand.T1(end-1);
hand.T1(end-1)=nan;
hand.T2(end-1)=nan;

C(hand.Cheq1.P)=single(hand.T2(hand.Cheq1.P_Xp)-hand.T1(1:end-2));
C(hand.Cheq2.P)=single(hand.T1(hand.Cheq2.P_Xp)-hand.T2(1:end-2));
C(~isfinite(C))=0;
if hand.Check_VaryD ==1
    hand.XFlux=sum(sum(C(1:end-1,:,:)./(hand.Dmap(1:end-1,:,:)+hand.Dmap(2:end,:,:)),3),2);
else
    hand.XFlux=sum(sum(C(1:end-1,:,:),3),2);
end
hand.XFlux_norm=flipud(hand.XFlux(2:end-1)/mean(hand.XFlux(2:end-1)));
hand.XFlux=flipud(hand.XFlux(2:end-1))/(hand.Qcv*(hand.L_X/(hand.L_Y*hand.L_Z)));
hand.T1(end-1)=em1;
hand.T2(end-1)=em1;


function [hand]=AfterIterate(hand)
%%% Future potential to remove voxels with one neighbour at steady state,
%%% then replace for visualisation
% T(hand.NN==1)=T(2:end-1,2:end-1,2:end-1)+T(2:end-1,2:end-1,2:end-1)+...
% T(2:end-1,2:end-1,2:end-1)+T(2:end-1,2:end-1,2:end-1)+...
% T(2:end-1,2:end-1,2:end-1)+T(2:end-1,2:end-1,2:end-1);
%[hand]=MBytesRecord(hand,whos,'Iterate complete'); %Memory
hand.Tconv=zeros(size(hand.Map),'single');
% Calculate_Button vertical flux
hand.T1(end-1)=nan;
hand.T2(end-1)=nan;
hand.Tconv(hand.Cheq1.P)=single(hand.T2(hand.Cheq1.P_Xp)-hand.T1(1:end-2));
hand.Tconv(hand.Cheq2.P)=single(hand.T1(hand.Cheq2.P_Xp)-hand.T2(1:end-2));
hand.Tconv(~isfinite(hand.Tconv))=0;

% if get(hand.Check_VaryD,'value')==1
%     hand.XFlux=sum(sum(hand.Tconv(1:end-1,:,:)./(hand.Dmap(1:end-1,:,:)+hand.Dmap(2:end,:,:)),3),2);
% else
%     hand.XFlux=sum(sum(hand.Tconv(1:end-1,:,:),3),2);
% end
% hand.XFlux_norm=flipud(hand.XFlux(2:end-1)/mean(hand.XFlux([2 end-1])));
% hand.XFlux=flipud(hand.XFlux(2:end-1)/hand.Qcv);
hand.Tconv=zeros(size(hand.Map),'double');
hand.Tconv(hand.Cheq1.P)=(hand.T1(1:end-2));
hand=rmfield(hand,'T1');
hand.Tconv(hand.Cheq2.P)=(hand.T2(1:end-2));
hand=rmfield(hand,'T2');
if hand.InLineMode==0
    if get(hand.Check_FluxMapSave,'Value')==1 && hand.impCheck==0
        TifSave3D(hand.Tconv(2:end-1,2:end-1,2:end-1),hand,'Cmap');
        z=whos('hand');
        if z.bytes/1024^3<3*2
            try
                save([hand.pathname,hand.fil,'_Cmap'],'-struct','hand','Tconv');
            catch
                disp('Insufficient permissions or space to save')
            end
        else
            try
                save([hand.pathname,hand.fil,'_Cmap'],'-struct','hand','Tconv', '-v7.3');
            catch
                disp('Insufficient permissions or space to save')
            end
        end
    end
end
hand.Tconv=single(hand.Tconv);
if hand.Check_VaryD ==0
    hand.Tconv(end,:,:)=single(2*hand.BotStim);
else
    hand.Tconv(end,:,:)=single(hand.BotStim);
end
if hand.Check_TauMode ==1
    hand.Cheq1=0;hand.Cheq2=0;hand=rmfield(hand,{'Cheq1','Cheq2'});
end
%hand.Results.SimTime=TimeString(toc/86400);
hand.Results.Iterations=hand.iter;
if hand.impCheck==0
    if sum(hand.Check_TauMode ==[1 2 4 5])
        hand.Results.Tau=mean([hand.TauFacBot(round(hand.iter/hand.check_f)),hand.TauFacTop(round(hand.iter/hand.check_f))]);
    elseif hand.Check_TauMode ==3
        hand.Results.Tau=hand.TauFacTop(end);
    elseif hand.Check_TauMode ==6
        hand.Results.Tau=hand.TauFacTop(end);
    end
end
%[hand]=MBytesRecord(hand,whos,'Tconv'); %Memory
hand.Results.MBytes=ceil(max(hand.MBytesA));
if hand.impCheck==0
    if hand.InLineMode==0
        try
            [hand]=InitiatePlot3(hand);
        catch
             disp('Unable to plot')
        end
    end
else
    if hand.Check_FreqPlots==1
        try
            [hand]=InitiatePlot3imp(hand);
        catch
             disp('Unable to plot')
        end
    end
end
%[hand]=MBytesRecord(hand,whos,'Iterate end'); %Memory



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function SaveResultsPDF(hand)
pdfstr=['_TFm',num2str(get(hand.Check_TauMode,'Value'))];
if get(hand.Check_Blocked, 'Value')==1
    pdfstr=[pdfstr,'CC'];
else
end
if get(hand.Check_Reverse, 'Value')==1
    pdfstr=[pdfstr,'F'];
else
end

if get(hand.Check_pdfSave,'Value')==1
    if get(hand.Check_RVA, 'Value')==0
        print(hand.ResFig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',...
            hand.Pha(1),'d',num2str(hand.Dir),pdfstr],'-dpdf');
    else
        print(hand.ResFig,[hand.pathname,hand.filename(1:find(hand.filename=='.')-1),'_p',...
            hand.Pha(1),'d',num2str(hand.Dir),'_RVA',num2str(round(100*hand.NetVol(end))),pdfstr],'-dpdf');
    end
end


function [Q_thro,Q_plane]=FluxHunter3ani(hand)
% hand.Map=single(hand.Map);
VoxDim=[hand.L_X, hand.L_Y, hand.L_Z];
DimFac(1)=VoxDim(2)*VoxDim(3)/VoxDim(1);
DimFac(2)=VoxDim(3)*VoxDim(1)/VoxDim(2);
DimFac(3)=VoxDim(1)*VoxDim(2)/VoxDim(3);
[a,b,c]=size(hand.Map(2:end-1,2:end-1,2:end-1));
if get(hand.Check_TauMode,'Value')==6
    %     hand.Tconv=imag(hand.Tconv);
end
if get(hand.Check_VaryD,'value')==0
    if c~=1 % && hand.MBytesA(end)*2<hand.RAM_availableMBs
        Q_thro=0.5*hand.Map(2:a+1,2:b+1,2:c+1)*DimFac(1).*...
            (hand.Map(1:a,  2:b+1,2:c+1).*abs(hand.Tconv(1:a,  2:b+1,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))...
            +hand.Map(3:a+2,2:b+1,2:c+1).*abs(hand.Tconv(3:a+2,2:b+1,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1)));
        Q_thro(1,:,:)=  Q_thro(1,:,:)+  DimFac(1)*hand.Map(2,2:end-1,2:end-1).*abs(hand.Tconv(2,2:end-1,2:end-1)-hand.TopStim);
        Q_thro(end,:,:)=Q_thro(end,:,:)+DimFac(1)*hand.Map(a+1,2:end-1,2:end-1).*abs(hand.BotStim-hand.Tconv(a+1,2:end-1,2:end-1));
        Q_plane=   0.5*hand.Map(2:a+1,2:b+1,2:c+1).*...
            (DimFac(2)*hand.Map(2:a+1,1:b,  2:c+1).*abs(hand.Tconv(2:a+1,1:b,2:c+1  )-hand.Tconv(2:a+1,2:b+1,2:c+1))...
            +DimFac(2)*hand.Map(2:a+1,3:b+2,2:c+1).*abs(hand.Tconv(2:a+1,3:b+2,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))...
            +DimFac(3)*hand.Map(2:a+1,2:b+1,1:c  ).*abs(hand.Tconv(2:a+1,2:b+1,1:c  )-hand.Tconv(2:a+1,2:b+1,2:c+1))...
            +DimFac(3)*hand.Map(2:a+1,2:b+1,3:c+2).*abs(hand.Tconv(2:a+1,2:b+1,3:c+2)-hand.Tconv(2:a+1,2:b+1,2:c+1)));
        if sum(get(hand.Check_TauMode,'Value')==[6])
            Q_thro(1,:,:)=  0.5*DimFac(1)*hand.Map(2,2:end-1,2:end-1).*hand.Map(3,2:end-1,2:end-1).*...
                abs(hand.Tconv(3,2:end-1,2:end-1)-hand.Tconv(2,2:end-1,2:end-1));
        end
        if sum(get(hand.Check_TauMode,'Value')==[2 3])
            Q_plane(:,1,:)=0.5*hand.Map(2:a+1,2,2:c+1).*...
                (DimFac(2)*hand.Map(2:a+1,b+1,2:c+1).*abs(hand.Tconv(2:a+1,b+1,2:c+1)-hand.Tconv(2:a+1,2  ,2:c+1))...
                +DimFac(2)*hand.Map(2:a+1,3  ,2:c+1).*abs(hand.Tconv(2:a+1,3  ,2:c+1)-hand.Tconv(2:a+1,2  ,2:c+1))...
                +DimFac(3)*hand.Map(2:a+1,2  ,[c+1,2:c]).*abs(hand.Tconv(2:a+1,2  ,[c+1,2:c])-hand.Tconv(2:a+1,2  ,2:c+1))...
                +DimFac(3)*hand.Map(2:a+1,2  ,[3:c+1,2]).*abs(hand.Tconv(2:a+1,2  ,[3:c+1,2])-hand.Tconv(2:a+1,2  ,2:c+1)));
            Q_plane(:,end,:)=0.5*hand.Map(2:a+1,b+1,2:c+1).*...
                (DimFac(2)*hand.Map(2:a+1,b  ,2:c+1).*abs(hand.Tconv(2:a+1,b  ,2:c+1)-hand.Tconv(2:a+1,b+1,2:c+1))...
                +DimFac(2)*hand.Map(2:a+1,2  ,2:c+1).*abs(hand.Tconv(2:a+1,2  ,2:c+1)-hand.Tconv(2:a+1,b+1,2:c+1))...
                +DimFac(3)*hand.Map(2:a+1,b+1,1:c  ).*abs(hand.Tconv(2:a+1,b+1,1:c  )-hand.Tconv(2:a+1,b+1,2:c+1))...
                +DimFac(3)*hand.Map(2:a+1,b+1,3:c+2).*abs(hand.Tconv(2:a+1,b+1,3:c+2)-hand.Tconv(2:a+1,b+1,2:c+1)));
            Q_plane(:,:,1)=0.5*hand.Map(2:a+1,2:b+1,2).*...
                (DimFac(2)*hand.Map(2:a+1,[b+1,2:b],2  ).*abs(hand.Tconv(2:a+1,[b+1,2:b],2  )-hand.Tconv(2:a+1,2:b+1,2  ))...
                +DimFac(2)*hand.Map(2:a+1,[3:b+1,2],2  ).*abs(hand.Tconv(2:a+1,[3:b+1,2],2  )-hand.Tconv(2:a+1,2:b+1,2  ))...
                +DimFac(3)*hand.Map(2:a+1,2:b+1,c+1).*abs(hand.Tconv(2:a+1,2:b+1,c+1)-hand.Tconv(2:a+1,2:b+1,2  ))...
                +DimFac(3)*hand.Map(2:a+1,2:b+1,3  ).*abs(hand.Tconv(2:a+1,2:b+1,3  )-hand.Tconv(2:a+1,2:b+1,2  )));
            Q_plane(:,:,end)=0.5*hand.Map(2:a+1,2:b+1,c+1).*...
                (DimFac(2)*hand.Map(2:a+1,[b+1,2:b],c+1).*abs(hand.Tconv(2:a+1,[b+1,2:b],c+1)-hand.Tconv(2:a+1,2:b+1,c+1))...
                +DimFac(2)*hand.Map(2:a+1,[3:b+1,2],c+1).*abs(hand.Tconv(2:a+1,[3:b+1,2],c+1)-hand.Tconv(2:a+1,2:b+1,c+1))...
                +DimFac(3)*hand.Map(2:a+1,2:b+1,c  ).*abs(hand.Tconv(2:a+1,2:b+1,c  )-hand.Tconv(2:a+1,2:b+1,c+1))...
                +DimFac(3)*hand.Map(2:a+1,2:b+1,2  ).*abs(hand.Tconv(2:a+1,2:b+1,2  )-hand.Tconv(2:a+1,2:b+1,c+1)));
            if get(hand.Check_TauMode,'Value')==3
                Q_thro(1,:,:)=0.5*hand.Map(2,2:b+1,2:c+1)*DimFac(1).*...
                    (hand.Map(a+1,  2:b+1,2:c+1).*abs(hand.Tconv(a+1,  2:b+1,2:c+1)-1-hand.Tconv(2,2:b+1,2:c+1))...
                    +hand.Map(3,2:b+1,2:c+1).*abs(hand.Tconv(3,2:b+1,2:c+1)-hand.Tconv(2,2:b+1,2:c+1)));
                Q_thro(end,:,:)=0.5*hand.Map(a+1,2:b+1,2:c+1)*DimFac(1).*...
                    (hand.Map(a,  2:b+1,2:c+1).*abs(hand.Tconv(a,  2:b+1,2:c+1)-hand.Tconv(a+1,2:b+1,2:c+1))...
                    +hand.Map(2,2:b+1,2:c+1).*abs(hand.Tconv(2,2:b+1,2:c+1)+1-hand.Tconv(a+1,2:b+1,2:c+1)));
            end
        end
    else %2D
        Q_thro=0.5*hand.Map(2:a+1,2:b+1,2).*...
            (DimFac(1)*hand.Map(1:a,2:b+1,2).*abs(  hand.Tconv(1:a,2:b+1,2)-hand.Tconv(2:a+1,2:b+1,2))...
            +DimFac(1)*hand.Map(3:a+2,2:b+1,2).*abs(hand.Tconv(3:a+2,2:b+1,2)-hand.Tconv(2:a+1,2:b+1,2)));
        Q_thro(1,:)  =Q_thro(1,:)+  DimFac(1)*hand.Map(2,2:end-1,2).*abs(hand.Tconv(2,2:end-1,2)-hand.TopStim);
        Q_thro(end,:)=Q_thro(end,:)+DimFac(1)*hand.Map(a+1,2:end-1,2).*abs(hand.Tconv(a+1,2:end-1,2)-hand.BotStim);
        
        Q_plane=0.5*hand.Map(2:a+1,2:b+1,2).*...
            (DimFac(2)*hand.Map(2:a+1,1:b,2).*abs(  hand.Tconv(2:a+1,1:b,2)-hand.Tconv(2:a+1,2:b+1,2))...
            +DimFac(2)*hand.Map(2:a+1,3:b+2,2).*abs(hand.Tconv(2:a+1,3:b+2,2)-hand.Tconv(2:a+1,2:b+1,2)));
        if sum(get(hand.Check_TauMode,'Value')==[6])
            Q_thro(1,:,:)=  0.5*DimFac(1)*hand.Map(2,2:end-1,2).*hand.Map(3,2:end-1,2).*...
                abs(hand.Tconv(3,2:b+1,2)-hand.Tconv(2,2:end-1,2));
        end
        if sum(get(hand.Check_TauMode,'Value')==[2 3])
            Q_plane(:,1)=0.5*hand.Map(2:a+1,2,2).*...
                (DimFac(2)*hand.Map(2:a+1,end-1,2).*abs(hand.Tconv(2:a+1,end-1,2)-hand.Tconv(2:a+1,2,2))...
                +DimFac(2)*hand.Map(2:a+1,3,2).*abs(hand.Tconv(2:a+1,3,2)-hand.Tconv(2:a+1,2,2)));
            Q_plane(:,end)=0.5*hand.Map(2:a+1,end-1,2).*...
                (DimFac(2)*hand.Map(2:a+1,end-2,2).*abs(hand.Tconv(2:a+1,end-2,2)-hand.Tconv(2:a+1,end-1,2))...
                +DimFac(2)*hand.Map(2:a+1,2,2).*abs(hand.Tconv(2:a+1,2,2)-hand.Tconv(2:a+1,end-1,2)));
            if get(hand.Check_TauMode,'Value')==3
                Q_thro(1,:,:)=0.5*hand.Map(2,2:b+1,2).*...
                    (DimFac(1)*hand.Map(a+1,2:b+1,2).*abs(  hand.Tconv(a+1,2:b+1,2)-1-hand.Tconv(2,2:b+1,2))...
                    +DimFac(1)*hand.Map(3,2:b+1,2).*abs(hand.Tconv(3,2:b+1,2)-hand.Tconv(2,2:b+1,2)));
                Q_thro(end,:,:)=0.5*hand.Map(a+1,2:b+1,2).*...
                    (DimFac(1)*hand.Map(a,2:b+1,2).*abs(  hand.Tconv(a,2:b+1,2)-hand.Tconv(a+1,2:b+1,2))...
                    +DimFac(1)*hand.Map(2,2:b+1,2).*abs(hand.Tconv(2,2:b+1,2)+1-hand.Tconv(a+1,2:b+1,2)));
            end
        end
    end
else
    Q_thro=2*hand.Map(2:a+1,2:b+1,2:c+1)*DimFac(1).*(...
        (hand.Map(1:a,  2:b+1,2:c+1).*abs(hand.Tconv(1:a,  2:b+1,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(1:a,  2:b+1,2:c+1)+ hand.Dmap(2:a+1,2:b+1,2:c+1)))+...
        (hand.Map(3:a+2,2:b+1,2:c+1).*abs(hand.Tconv(3:a+2,2:b+1,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(3:a+2,2:b+1,2:c+1)+ hand.Dmap(2:a+1,2:b+1,2:c+1))));
    Q_thro([1 end],:,:)=2*hand.Map([2 end-1],2:b+1,2:c+1)*DimFac(1).*(...
        (hand.Map([2 end-2],2:b+1,2:c+1).*abs(hand.Tconv([1 end-2],  2:b+1,2:c+1)-hand.Tconv([2 end-1],2:b+1,2:c+1))./...
        (hand.Dmap([1 end-2],  2:b+1,2:c+1)+ hand.Dmap([2 end-1],2:b+1,2:c+1)))+...
        (hand.Map([3 end-1],2:b+1,2:c+1).*abs(hand.Tconv([3 end  ],2:b+1,2:c+1)-  hand.Tconv([2 end-1],2:b+1,2:c+1))./...
        (hand.Dmap([3 end  ],2:b+1,2:c+1)+   hand.Dmap([2 end-1],2:b+1,2:c+1))));
    
    Q_plane=2*hand.Map(2:a+1,2:b+1,2:c+1).*(...
        (DimFac(2)*hand.Map(2:a+1,1:b  ,2:c+1).*abs(hand.Tconv(2:a+1,1:b,2:c+1)-  hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(2:a+1,1:b,2:c+1)+   hand.Dmap(2:a+1,2:b+1,2:c+1)))+...
        (DimFac(2)*hand.Map(2:a+1,3:b+2,2:c+1).*abs(hand.Tconv(2:a+1,3:b+2,2:c+1)-hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(2:a+1,3:b+2,2:c+1)+ hand.Dmap(2:a+1,2:b+1,2:c+1)))+...
        (DimFac(3)*hand.Map(2:a+1,2:b+1,1:c  ).*abs(hand.Tconv(2:a+1,2:b+1,1:c)-  hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(2:a+1,2:b+1,1:c)+   hand.Dmap(2:a+1,2:b+1,2:c+1)))+...
        (DimFac(3)*hand.Map(2:a+1,2:b+1,3:c+2).*abs(hand.Tconv(2:a+1,2:b+1,3:c+2)-hand.Tconv(2:a+1,2:b+1,2:c+1))./...
        (hand.Dmap(2:a+1,2:b+1,3:c+2)+ hand.Dmap(2:a+1,2:b+1,2:c+1))));
    %     else
    %         Q_thro=0.5*hand.Map(2:a+1,2:b+1,2).*...
    %             (DimFac(1)*hand.Map(1:a,2:b+1,2).*abs(  hand.Tconv(1:a,2:b+1,2)-hand.Tconv(2:a+1,2:b+1,2))...
    %             +DimFac(1)*hand.Map(3:a+2,2:b+1,2).*abs(hand.Tconv(3:a+2,2:b+1,2)-hand.Tconv(2:a+1,2:b+1,2)));
    %         Q_thro(1,:)  =Q_thro(1,:)+  DimFac(1)*hand.Map(2,2:end-1,2).*abs(hand.Tconv(2,2:end-1,2)-hand.TopStim);
    %         Q_thro(end,:)=Q_thro(end,:)+DimFac(1)*hand.Map(a+1,2:end-1,2).*abs(hand.Tconv(a+1,2:end-1,2)-hand.BotStim);
    %
    %         Q_plane=0.5*hand.Map(2:a+1,2:b+1,2).*...
    %             (DimFac(2)*hand.Map(2:a+1,1:b,2).*abs(  hand.Tconv(2:a+1,1:b,2)-hand.Tconv(2:a+1,2:b+1,2))...
    %             +DimFac(2)*hand.Map(2:a+1,3:b+2,2).*abs(hand.Tconv(2:a+1,3:b+2,2)-hand.Tconv(2:a+1,2:b+1,2)));
    %     end
end
Q_plane=real(Q_plane);
Q_thro=real(Q_thro);



function [hand]=RVA_Net_Cube(hand)
%% Calculated the idx for cubic RVA
[a,b,c]=size(hand.Net_full);
s=(1-hand.shrink)^(1/3);
d=(1-s)/2;
hand.RVAdims=[...
    floor(d*a+1),ceil(a*(1-d));...
    floor(d*b+1),ceil(b*(1-d));...
    floor(d*c+1),ceil(c*(1-d))];

function [hand]=RVA_Net_Lconst(hand)
%% Calculated the idx for Lconst RVA
[a,b,c]=size(hand.Net_full);
s=(1-hand.shrink)^(1/2);
d=(1-s)/2;
hand.RVAdims=[...
    1,a;...
    floor(d*b+1),ceil(b*(1-d));...
    floor(d*c+1),ceil(c*(1-d))];

function [hand]=RVA_Net_Aconst(hand)
%% Calculated the idx for Aconst RVA
[a,b,c]=size(hand.Net_full);
s=1-hand.shrink;
d=(1-s);
if get(hand.Pop_RV, 'Value')==3
    hand.RVAdims=[...
        1,ceil(a*(1-d));...
        1,b;...
        1,c];
else
    hand.RVAdims=[...
        floor(a*d+1),a;...
        1,b;...
        1,c];
end


function [hand]=MBytesRecord(hand,CurrentWoSpa,PlaceName)
hand.MBytesA=[hand.MBytesA sum([CurrentWoSpa.bytes])/1024^2]; %Memory
if nargin<3
    PlaceName='';
end
hand.MemLoc{end+1}=PlaceName;

%% USED IN THE CALCULATIONS %%
function [tocStr]=TimeString(tocking)
if tocking<1/86400
    tocStr = ['= 1 sec'];
elseif tocking<100/86400
    tocStr = ['= ',num2str(round(tocking*86400)),' sec'];
elseif tocking<100*60/86400
    tocStr = ['= ',num2str(round(tocking*86400/60)),' min'];
elseif tocking<36*3600/86400
    tocStr = ['= ',num2str(round(tocking*86400/3600*2)/2),' hrs'];
else
    tocStr = ['= ',num2str(ceil(tocking*2)/2),' days'];
end
