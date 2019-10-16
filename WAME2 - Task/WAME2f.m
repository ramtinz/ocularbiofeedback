function varargout = WAME2f(varargin)
% WAME2f MATLAB code for WAME2f.fig
%      WAME2f, by itself, creates a new WAME2f or raises the existing
%      singleton*.
%
%      H = WAME2f returns the handle to a new WAME2f or the handle to
%      the existing singleton*.
%
%      WAME2f('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAME2f.M with the given input arguments.
%
%      WAME2f('Property','Value',...) creates a new WAME2f or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before WAME2f_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to WAME2f_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help WAME2f

% Last Modified by GUIDE v2.5 07-Sep-2018 11:43:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @WAME2f_OpeningFcn, ...
    'gui_OutputFcn',  @WAME2f_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before WAME2f is made visible.
function WAME2f_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to WAME2f (see VARARGIN)

% Choose default command line output for WAME2f
handles.output = hObject;
% Initialization
set(handles.guifig, 'Units', 'pixels', 'Position', [1281+188 1+145  905   717]);
set(handles.axes1, 'Units', 'pixels', 'Position', [17.0666, 17.0666, 682.6660   682.6660]); %[left bottom width height]
set(handles.axes2, 'Units', 'pixels', 'Position', [717, 529.0661, 170.6665, 170.6665]);
set(handles.text4, 'Units', 'pixels', 'Position', [717, 480, 129, 39]);
set(handles.text3, 'Units', 'pixels', 'Position', [717, 425, 136, 57]);
% set(handles.checkbox3, 'Units', 'pixels', 'Position', [717, 5*34.1333, 169, 39]);
set(handles.startbutton, 'Units', 'pixels', 'Position', [717, 3*34.1333, 59, 20]);
set(handles.edit1, 'Units', 'pixels', 'Position', [717, 34.1333, 66, 21]);
set(handles.text4, 'Visible','off')

handles.UserData.fdbck_num=0; % initial values for feedback numbers
handles.UserData.label = zeros(1,12); % initial value of class prediction, 1 fatigue 0 alert


handles.UserData.acm_cyc_idx = []; % Accumulative cycle sample index
handles.UserData.Metrics = nan(5,12); % all variables (metrics) stored here
handles.UserData.BreakReqTime = nan(1,9);
% handles.UserData.tot = 1; % initial value for time on task number
handles.UserData.Bchoice = zeros(1,9); % initial value for number of breaks
global C_cnt etime start_WP start_MP start_RP stop_RP stop_WP c_clk smpl_cnt cyc_str_smpl cyc_stp_smpl BRK procXcuted sgmnt
C_cnt = 0;
start_WP = zeros(1,250);
start_MP = zeros(1,250);
start_RP = zeros(1,250);
stop_RP = zeros(1,250);
stop_WP = zeros(1,250);
%%%%%%%%%%%%%%%
cyc_str_smpl = nan(1,250);
cyc_stp_smpl = nan(1,250);
c_clk = 1;
etime = tic;
smpl_cnt = 1;
BRK=0;
procXcuted=0;
sgmnt=1;

handles.UserData.no_ratings = 0; % Number of ratings
 % Number of successive cycles (For KSS's self-assessment)
set(handles.KSS, 'Units', 'centimeters', 'Position',[8, 8, 10, 7]);
movegui(handles.guifig,'center');

handles.UserData.ic_clks_cycl(1) = 0; % incorrect clicks per cycle
handles.UserData.ic_clks_time(1).ttl_clks(1) = nan;  % incorrect clicks times
handles.UserData.c_clks_cycl(1) = 0; % correct clicks per cycle
handles.UserData.c_clks_time(1).ttl_clks(1) = nan; % correct clicks times
handles.ttlclks_per_cycle = 1; % total clicks per cycle

set(handles.edit1,'enable','on') % subject ID box
handles.UserData.c_index = 1;

tskmode_prompt = 'Choose manual 1, automatic 2, or training 0: ';
handles.UserData.tskmode = input(tskmode_prompt);
while isempty(handles.UserData.tskmode) || (handles.UserData.tskmode~=0 && handles.UserData.tskmode~=1 && handles.UserData.tskmode~=2)
    handles.UserData.tskmode = input('Invalid! Choose manual 1, automatic 2, or training 0: ');
end
Choose_day = input('Choose the day (1, 2 or 3): ');
while isempty(Choose_day) || (Choose_day~=1 && Choose_day~=2 && Choose_day~=3)
    Choose_day = input('Invalid! Choose the day (1, 2 or 3): ');
end
handles.UserData.Day = Choose_day;



if Choose_day == 1
    addpath(genpath('Day1_pnts'))
    handles.UserData.MaxCycNum = 61; % ~10 min (first day training)
elseif Choose_day == 2
    addpath(genpath('Day2_pnts'))
    if handles.UserData.tskmode==0
        handles.UserData.MaxCycNum = 45; % ~5 min (warm up task)
    else
        handles.UserData.MaxCycNum = 182;
    end
else
    addpath(genpath('Day3_pnts'))
    if handles.UserData.tskmode==0
        handles.UserData.MaxCycNum = 45; % ~5 min (warm up task)
    else
        handles.UserData.MaxCycNum = 182;
    end
end
if handles.UserData.tskmode==0 % train mode
    handles.UserData.predpoints = load('pts_med_train');
else % test mode
    handles.UserData.predpoints = load('pts_med');
end
set(handles.axes1,'hittest','off')
handles.UserData.ic_clks_cycl(1:250) = 0;
handles.UserData.c_clks_cycl(1:250) = 0;
handles.UserData.ic_clks_time(250).ttl_clks = 0;
handles.UserData.c_clks_time(250).ttl_clks = 0;

for i = 1:5
    hold(handles.axes2,'on')
    plot(rand(1,2),'LineWidth',1,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','visible','off','Parent',handles.axes2);
end
for i = 1:6
    hold(handles.axes1,'on')
    plot(rand(1,2),'LineWidth',1,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','visible','off','Parent',handles.axes1);
end
handles.markers = {'+','o','*','s','d','p','^'};
handles.markers_desc = {'Plus','Circle','Asterisk','Square','Diamond','Pentagram','triangle'};


hold(handles.axes1,'on')
handles.Hp5 = plot(5,5,'LineWidth',2,'MarkerSize',10,'Marker','x','MarkerFaceColor','k','MarkerEdgeColor','k','Parent',handles.axes1);
hold off
set(handles.axes1,'Xlim',[0 10]','Ylim',[0 10]','XTickLabel',[],'YTickLabel',[]);
set(handles.axes2,'Xlim',[0 10]','Ylim',[0 10]','XTickLabel',[],'YTickLabel',[],'Color','None');


set(get(handles.axes1,'Children'),'ButtonDownFcn',@(hObject,eventdata)WAME2f('axes1_ButtonDownFcn',hObject,eventdata,guidata(hObject)));

set(handles.KSS,'visible','off');
% set(handles.AskforBreak,'visible','off');
set(handles.Hp5,'visible','off');
% which says the program is set to relaiability mode when the checkbox is checked
% that is when this flag gets the one value and if not for the default value it gets
% the zero value instead

% mouse cursor adjustment
handles.pshap = get(handles.guifig,'PointerShapeCData');
set(handles.guifig, 'Pointer', 'hand')
% import java.awt.Robot;
% handles.mouse = Robot;

sel_option = get(handles.KSS,'SelectedObject');
set(sel_option, 'Value', 0)

handles.MP_tmr = timer(...
    'ExecutionMode', 'singleShot', ...
    'StartDelay', 0 ,... % 9.74-2*2.34
    'Busymode','queue',...
    'TimerFcn', {@MP_fcn, hObject}); % 5.06
handles.WP_tmr = timer(...
    'ExecutionMode', 'singleShot', ...
    'StartDelay', 2.34,...
    'Busymode','queue',...
    'TimerFcn', {@WP_fcn, hObject}); % 2.34
handles.RP_tmr = timer(...
    'ExecutionMode', 'singleShot', ...
    'StartDelay', 2.34, ...
    'Busymode','queue',...
    'TimerFcn', {@RP_fcn, hObject}); % 2.34
handles.KSS_tmr = timer(...
    'ExecutionMode', 'singleShot', ...
    'StartDelay', 0,...
    'Busymode','queue',...
    'TimerFcn', {@KSS_disp, hObject}); % 0
handles.KSSread_tmr = timer(...
    'ExecutionMode', 'singleShot', ...
    'StartDelay', 5,...
    'Busymode','queue',...
    'TimerFcn', {@KSS_read, hObject}); % 5
handles.BRKMSG = timer(...
    'ExecutionMode', 'singleShot', ...
    'StartDelay', 0,...
    'Busymode','queue',...
    'TimerFcn', {@BRKMSG_disp, hObject}); % 0


% The following lines are for communication between this computer and
if handles.UserData.tskmode~=0
    addpath(genpath('ET_RTcom'))
    % ExclDtItm=fullfile(CodDr,'DataItems2.xlsx');
    [~,handles.UserData.DtItmsTx,~]=xlsread('DataItems2.xlsx');
    SFct=handles.UserData.DtItmsTx(2:end,5); handles.UserData.SclFct=str2double(SFct); % deriving the scale factor
    handles.t_rzm = tcpip('169.254.181.51', 51000);
    set(handles.t_rzm,'ByteOrder','littleEndian');
    handles.t_rzm.InputBufferSize = 100000;
    guidata(hObject,handles);
    fopen(handles.t_rzm);
    pause(1);
    SOCKET_TYPE_SDATA_TCP  = 3; % Send Data to Remote via TCP / IP
    
    ET7_SetConnectType(handles.t_rzm, SOCKET_TYPE_SDATA_TCP);
    
    handles.t_2 = tcpip('169.254.181.51', 51000); % make a new tcpip object
    set(handles.t_2,'ByteOrder','littleEndian'); % set tcpip specifications
    handles.t_2.InputBufferSize = 100000; % set tcpip specifications
    guidata(hObject,handles);
    fopen(handles.t_2); % open tcpip port
    pause(1)
    ET7_OpenDataFile(handles.t_rzm); % opens a new data file on ET7 PC
    
    
    
    handles.data_rec_tmr = timer(...
        'ExecutionMode', 'fixedRate', ...
        'Period', 0.5, ...
        'Busymode','queue',...
        'TimerFcn', {@data_rec, hObject}); % Data Reseive
    
    handles.UserData.LftDt = [];
    
    
    handles.UserData.MsgBlkSize=86;
    handles.UserData.ET_data = nan(1000000,9); % where the eye data will be stored
    handles.UserData.acm_cyc_idx = []; % Accumulative cycle sample index
    C_cnt = 0;
    load Mdl_EnsTree % load fatigue state classification model
    handles.UserData.fmdl = Mdl;
    load ppl_fltr % load smoothing filter for the pupil diameter signal
    handles.UserData.P_fltr = d1;
    
    
end
% This function has codes to start recording on eye tracker which is
% provided by ASL (Finish recording is added to closing function)
% -----------

% Update handles structure
guidata(hObject,handles);

% UIWAIT makes WAME2f wait for user response (see UIRESUME)
% uiwait(handles.guifig);



% --- Outputs from this function are returned to the command line.
function varargout = WAME2f_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




function axes1_CreateFcn(hObject, eventdata, handles)

function axes1_DeleteFcn(hObject, eventdata, handles)






function startbutton_Callback(hObject, eventdata, handles)
global etime
% set(handles.checkbox3,'visible','off')
handles.UserData.Onset_date_time = datestr(now);
handles.UserData.task_starting_time = toc(etime);
% handles.UserData.absolute_etime = tic;
% The following three lines are for communication between this computer and
if handles.UserData.tskmode~=0
    ET7_StartDataFileRecording(handles.t_rzm)
    handles.UserData.ET_rec_start_time = toc(etime);
    
end
% This function has codes to start recording on eye tracker which is
% provided by ASL (Finish recording is added to closing function)
% -----------

SbjID = get(handles.edit1,'String'); % getting subject ID from subject ID box
set(handles.edit1,'Visible','Off'); % hiding subject ID box
set(handles.startbutton,'Visible','Off');  % hiding start button
handles.UserData.SbjID = SbjID;

guidata(hObject, handles);
start(handles.KSS_tmr);

function data_rec(~,~,hObject,~)
handles = guidata(hObject);
global C_cnt smpl_cnt cyc_str_smpl cyc_stp_smpl BRK procXcuted sgmnt

if handles.UserData.tskrning==1 % record data while the task is running (not during KSS or breaks)
    if handles.t_2.BytesAvailable>0
        data1=[handles.UserData.LftDt;fread(handles.t_2, handles.t_2.BytesAvailable)];
        [handles.UserData.LftDt,OutVrb]=ParsDt(data1,handles.UserData.MsgBlkSize,handles.UserData.DtItmsTx,handles.UserData.SclFct);
        [SmpNm,~]=size(OutVrb);
        handles.UserData.ET_data(smpl_cnt+1:smpl_cnt+SmpNm,:)=OutVrb(:,[7 8 9 10 11 12 13 14 15]); % with head-tracker
        smpl_cnt=smpl_cnt+SmpNm;
    end
end
if rem(C_cnt,20)==1 && C_cnt~=1 && procXcuted==0
    procXcuted=1; % to execute the processing one for each cycle after the cycle of 20
    handles.UserData.acm_cyc_idx = cyc_str_smpl(C_cnt-20):cyc_stp_smpl(C_cnt-1);
    EH_gaze_hcoord = fillmissing(handles.UserData.ET_data(handles.UserData.acm_cyc_idx,5),'linear');
    EH_gaze_vcoord = fillmissing(handles.UserData.ET_data(handles.UserData.acm_cyc_idx,6),'linear');
    EH_gaze_length = fillmissing(handles.UserData.ET_data(handles.UserData.acm_cyc_idx,4),'linear');
    PplDiam = handles.UserData.ET_data(handles.UserData.acm_cyc_idx,1);
    EH_gd_x = fillmissing(handles.UserData.ET_data(handles.UserData.acm_cyc_idx,7),'linear');
    EH_gd_y = fillmissing(handles.UserData.ET_data(handles.UserData.acm_cyc_idx,8),'linear');
    EH_gd_z = fillmissing(handles.UserData.ET_data(handles.UserData.acm_cyc_idx,9),'linear');
    
    cr_diam = handles.UserData.ET_data(handles.UserData.acm_cyc_idx,2); % corneal diameter
    %     When it goes to zero or nan it means it is not detecte
    scen_num = handles.UserData.ET_data(handles.UserData.acm_cyc_idx,3); % scene number
    % It has to be zero for the default scene defined
    EH_gaze_hcoord = medfilt1(EH_gaze_hcoord,3);
    EH_gaze_vcoord = medfilt1(EH_gaze_vcoord,3);
    EH_gaze_length = medfilt1(EH_gaze_length,3);
    EH_gd_x = medfilt1(EH_gd_x,3);
    EH_gd_y = medfilt1(EH_gd_y,3);
    EH_gd_z = medfilt1(EH_gd_z,3);
    model_pd = 77; % typical value for model_pd CHANGE
    PplDiam = PplDiam.*3.96/model_pd; % Pupil size in mm
    u_idx = find(PplDiam~=0 & isnan(PplDiam)==0);
    b_samples = find(PplDiam==0 | isnan(PplDiam)); %mis_idx';
    closd_smpl = b_samples; % closed eyes samples
    D = diff([0,diff(closd_smpl')==1,0]);
    ce_on = closd_smpl(D>0);
    ce_off = closd_smpl(D<0);
    ce_d = (1 + find(D<0) - find(D>0))/360;
    ce_on(ce_d<0.1 | ce_d>1)=[];ce_off(ce_d<0.1 | ce_d>1)=[];
    closd_smpl = arrayfun(@colon, ce_on, ce_off, 'Uniform', false);
    closeys_idx = cell2mat(closd_smpl');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    win=10;mis_idx=[];
    i=win+1;
    while i<length(PplDiam)
        k=1;
        while PplDiam(i+k)==0 && (i+k)<length(PplDiam)
            k=k+1;
        end
        if k>1 && i+k+win<length(PplDiam)
            mis_idx = [mis_idx i-win:i+k-1 i+k:i+k+win];
            i=i+win;
        end
        i=i+k;
    end
    mis_idx=unique([mis_idx find(PplDiam<(median(PplDiam(u_idx))-3*std(PplDiam(u_idx))))']);
    PplDiam(mis_idx)=nan;
    PplDiam=fillmissing(PplDiam,'linear');
    PplDiam = filtfilt(handles.UserData.P_fltr,PplDiam);
    PplDiam((PplDiam>median(PplDiam)+1) | (PplDiam<median(PplDiam)-1))=median(PplDiam);
    L = length(EH_gaze_hcoord);
    D = diff([0,diff(b_samples')==1,0]);
    b_on = b_samples(D>0);
    b_off = b_samples(D<0);
    b_d = (1 + find(D<0) - find(D>0))/360;
    b_outrng = b_d<0.1 | b_d>.4;
    undef_idx = arrayfun(@colon, b_on, b_off, 'Uniform', false);
    undef_idx = cell2mat(undef_idx');
    b_on(b_outrng)=[];b_off(b_outrng)=[];b_d(b_outrng)=[];
    undef_idx = unique([undef_idx find(scen_num~=0)']);
    
    b_samples = arrayfun(@colon, b_on, b_off, 'Uniform', false);
    b_samples = cell2mat(b_samples');
    v = [EH_gd_x';EH_gd_y';EH_gd_z'];
    theta = nan(1,L-1);
    for i=1:L-1
        if v(:,i)==zeros(3,1)
            v(:,i)=[0.0001;0.0001;0.0001];
        end
        if v(:,i+1)==zeros(3,1)
            v(:,i+1)=[0.0001;0.0001;0.0001];
        end
        theta(i)=acosd((v(:,i)'*v(:,i+1))/(norm(v(:,i))*norm(v(:,i+1))));
    end
    theta = abs(theta);
    [~,g] = sgolay(2,19);
    
    dx = zeros(L,2);
    dy = zeros(L,2);
    dt = 1/360;
    for p = 1:2
        dx(:,p) = conv(EH_gaze_hcoord, factorial(p)/(-dt)^p * g(:,p+1), 'same');
        dy(:,p) = conv(EH_gaze_vcoord, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    end
    for i=1:L-1
        Phi(i) = 2*atand(1/(2*(median(EH_gaze_length(i:i+1)))));
    end
    for i=1:L-1
        dtheta(i) = Phi(i)*sqrt(dx(i,1).^2+dy(i,1).^2);
        ddtheta(i) = Phi(i)*sqrt(dx(i,2).^2+dy(i,2).^2);
    end
    %--------------------------------------------------------------
    % disp(L)
    theta = theta';
    theta(L)=theta(L-1);
    dtheta(L)=dtheta(L-1);
    ddtheta(L)=ddtheta(L-1);
    %         err_smpl = 1:L;
    err_smpl = find(dtheta>500 | ddtheta>50000);
    int_err_smpl=[];
    for j=1:length(err_smpl)-1
        if err_smpl(j+1)-err_smpl(j)<20
            int_err_smpl = [int_err_smpl err_smpl(j):err_smpl(j+1)];
        end
    end
    err_smpl = sort([err_smpl int_err_smpl]);
    err_smpl(diff(err_smpl)==0)=[];
    D = diff([0,diff(err_smpl)==1,0]);
    err_smpl_on = err_smpl(D>0);
    err_smpl_off = err_smpl(D<0);
    err_smpl_d = (1 + find(D<0) - find(D>0))/360;
    long_err_smpl=[];
    err_smpl_on(err_smpl_d<0.020)=[];err_smpl_off(err_smpl_d<0.020)=[];
    long_err_smpl = arrayfun(@colon, err_smpl_on, err_smpl_off, 'Uniform', false);
    long_err_smpl = cell2mat(long_err_smpl);
    PT(1) = 100; % initial value for saccade onset detection
    PT(2) = mean(dtheta(dtheta<PT(1)))+6*std(dtheta(dtheta<PT(1)));
    i=3;
    while abs(PT(i-1)-PT(i-2))>1
        PT(i)= mean(dtheta(dtheta<PT(i-1)))+6*std(dtheta(dtheta<PT(i-1)));
        i=i+1;
    end
    s_onset_th = min(mean(dtheta(dtheta<PT(end))) + 3*std(dtheta(dtheta<PT(end))),45);
    %                 s_smpl = find(dtheta>s_onset_th);
    s_smpl = find(dtheta>max(s_onset_th,25));
    if isempty(s_smpl)==0
        s_smpl(ismembc(s_smpl,b_samples))=[];
    end
    if isempty(err_smpl)==0
        s_smpl(ismembc(s_smpl,err_smpl))=[];
    end
    s_smpl(ismembc(s_smpl,find(isnan(cr_diam) | cr_diam==0)))=[];
    %     D = [];
    D = diff([0,diff(s_smpl)==1,0]);
    s_on = s_smpl(D>0);
    s_off = s_smpl(D<0);
    s_d = (1 + find(D<0) - find(D>0))/360;
    s_rmv = find(s_d<0.020 | s_d>0.200);
    s_on(s_rmv)=[];s_off(s_rmv)=[];s_d(s_rmv)=[];
    s_on=s_on';s_off=s_off';
    s_pv = arrayfun(@(s, e) max(dtheta(s:e)), s_on, s_off);
    s_a = arrayfun(@(s,e) (mean(dtheta(s:e))).*length(s:e)/360, s_on, s_off);
    rmv_idx = find(sqrt(s_pv-(24+26*s_a))>5 | s_a>20);
    s_pv(rmv_idx)=[];s_a(rmv_idx)=[];s_on(rmv_idx)=[];s_off(rmv_idx)=[];s_d(rmv_idx)=[];
    
    undef_idx = unique([undef_idx rmv_idx' s_rmv]);  % index of undefined samples
    % Fixations
    f_samples = 1:L;
    cr_rmv = find(isnan(cr_diam) | cr_diam==0);
    f_rmv = unique([s_smpl,b_samples,long_err_smpl,undef_idx,cr_rmv']);
    f_samples(f_rmv)=[];
    
    D = diff([0,diff(f_samples)==1,0]);
    f_on = f_samples(D>0);
    f_off = f_samples(D<0);
    if isempty(f_on)==0
        f_rmv = find((f_on(2:end)-f_off(1:end-1))<5); % Minimum time between adjacent fixations to combine (11 ms) as the same as minimum saccade duration
        %         % But this says (75 ms) which is not accurate enough Ref: Olsen, A. (2012). The Tobii I-VT fixation filter. Tobii Technology.
        
        f1= f_on(1);f2=f_off(end);
        f_on(1)=[];f_off(end)=[];
        f_on(f_rmv)=[];f_off(f_rmv)=[];
        f_on = [f1 f_on];f_off = [f_off f2];
        
        f_samples = arrayfun(@colon, f_on, f_off, 'Uniform', false);
        f_samples = cell2mat(f_samples);
        
        D = diff([0,diff(f_samples)==1,0]);
        f_d = (1 + find(D<0) - find(D>0))/360;
        f_rmv2 = find(f_d<0.03 | f_d>2.5); % minimum fixation duration (40 ms)
        f_on(f_rmv2)=[];f_off(f_rmv2)=[];f_d(f_rmv2)=[];
        %     f_samples=[];
        f_samples = arrayfun(@colon, f_on, f_off, 'Uniform', false);
        f_samples = cell2mat(f_samples);
        
        
        rmv = [];
        for j=1:numel(f_on)
            C_GP_dist = sqrt((EH_gaze_hcoord(f_on(j):f_off(j))-median(EH_gaze_hcoord(f_on(j):f_off(j)))).^2+(EH_gaze_vcoord(f_on(j):f_off(j))-median(EH_gaze_vcoord(f_on(j):f_off(j)))).^2);
            if sum(C_GP_dist>1)>11
                rmv = [rmv j];
            end
        end
        %                     disp(numel(rmv))
        f_on(rmv)=[];f_off(rmv)=[];f_d(rmv)=[];
        f_samples=[];
        f_samples = arrayfun(@colon, f_on, f_off, 'Uniform', false);
        f_samples = cell2mat(f_samples);
        %     disp([6 10])
        undef_idx = unique([undef_idx cr_rmv' f_rmv f_rmv2]);
    else
        f_d=nan;%f_on=nan;f_off=nan;
    end
    
    if numel(s_on)>50
        [Line_p,~] = robustfit(s_a,s_pv,'welsch');
        Slope = Line_p(2);
    else
        Slope = 24;
    end
    BF=360*numel(b_on)/(L-numel(undef_idx)); % (cyc_stp_smpl-cyc_str_smpl)
    
    S_PVA=Slope;
    SF=360*(numel(s_on))/(L-numel(undef_idx));
    PDR=iqr(PplDiam(f_samples)); %mean(pdm,'omitnan');
    PERCLOS = numel(closeys_idx)/L;%(L-numel(undef_idx)); % Percentage of closed eyes to opened eyes
    
    
    if isempty(S_PVA) || isinf(S_PVA) || isnan(S_PVA)
        S_PVA=0;
    end
    if isempty(PDR) || isinf(PDR) || isnan(PDR)
        PDR=0;
    end
    if isempty(SF) || isinf(SF) || isnan(SF)
        SF=0;
    end
    if isempty(PERCLOS) || isinf(PERCLOS) || isnan(PERCLOS)
        PERCLOS=0;
    end
    if isempty(BF) || isinf(BF) || isnan(BF)
        BF=0;
    end
    
    X = [SF PERCLOS BF S_PVA PDR];
    
    % BF: blink frequency
    % SF: saccade frequency
    % S_PVA: Saccade Peak Velocity Amplitude Relationship
    % SCD: Saccade Duration
    % PDR: Pupil Dilation Interquartile Range
    % PERCLOS: Pecentage of eyes closed
    
    disp(X)
    lab = predict(handles.UserData.fmdl,X);
    BRK=0;
    disp(lab)
    if lab=='F'
        handles.UserData.label(sgmnt) = 1; % Fatigue
        BRK=1;
    end
    handles.UserData.Metrics(:,sgmnt) = X';
    disp(sgmnt)
    disp(X)
    sgmnt=sgmnt+1;
    
end
% end
guidata(hObject, handles);


function MP_fcn(obj,event,hObject,eventdata)
handles = guidata(hObject);
global C_cnt etime start_MP stop_RP c_clk smpl_cnt cyc_str_smpl cyc_stp_smpl
% global c_clk etime start_MP stop_RP sampl_cnt
set(handles.axes1,'hittest','off') % clicking prevention during WP and MP
set(handles.axes1.Children,'hittest','on')
% setappdata(handles,'C_cnt',C_cnt+1)
C_cnt = C_cnt+1;
% guidata(hObject, handles);
if C_cnt==handles.UserData.MaxCycNum % If number of tasks reaches MaxCycNum stop the task and close
    disp('close requested')
    %         guidata(hObject, handles);
    close(handles.guifig)
    return
end
disp(C_cnt)
cyc_str_smpl(C_cnt)=smpl_cnt;
if C_cnt>1
    cyc_stp_smpl(C_cnt-1)=cyc_str_smpl(C_cnt)-1;
end
if C_cnt>1
    stop_RP(C_cnt-1) = toc(etime);
%     disp(stop_RP(C_cnt-1)-start_MP(C_cnt-1))
end
delete(handles.axes1.Children(1:length(handles.axes1.Children)))
c_clk = 1;
handles.UserData.clkd_pnt{handles.ttlclks_per_cycle,C_cnt} = nan;
if (C_cnt-1)/20==handles.UserData.no_ratings
    start(handles.KSS_tmr);
    
else
    
    set(handles.guifig, 'Pointer', 'custom', 'PointerShapeCData', NaN(16,16))
    
    shuff_indx = circshift([1:7]',[rem(C_cnt-1,7) 0]);
    
    set(handles.text3,'String',handles.markers_desc{shuff_indx(1)});
    set(handles.text3,'Visible','On'); % showing "The starting point is:"
    % L = length(handles.axes1.Children);
    
    pnts=handles.UserData.predpoints;
    handles.x = pnts.pnts(C_cnt).iteration(1,:);
    handles.y = pnts.pnts(C_cnt).iteration(2,:);
    for i = 1:5
        hold(handles.axes2,'on')
        plot(handles.x(i),handles.y(i),'LineWidth',1,'Marker',handles.markers{shuff_indx(i)},'MarkerFaceColor','k','MarkerEdgeColor','k','visible','on', 'Parent', handles.axes2);
    end
    
    hold(handles.axes2,'on')
    
    line(handles.axes2,handles.x(1:5),handles.y(1:5),'LineWidth',0.5,'Color','k');
    
    
    set(handles.text4, 'Visible','on') % showing the name of the starting point
    
    start_MP(C_cnt) = toc(etime);
    handles.UserData.tskrning = 1;
end
guidata(hObject, handles);
if strcmp(get(handles.KSS_tmr,'Running'),'on') || strcmp(get(handles.KSSread_tmr,'Running'),'on')
    stop(handles.MP_tmr)
    set(handles.MP_tmr,'StartDelay', 0)
else
    stop(handles.MP_tmr)
    start(handles.WP_tmr);
end



function WP_fcn(obj,event,hObject,eventdata)
handles = guidata(hObject);
global C_cnt etime start_WP
% disp(C_cnt)
% global etime start_WP

set(handles.guifig, 'Pointer', 'custom', 'PointerShapeCData', NaN(16,16))
set(handles.text4, 'Visible','off')
movegui(handles.guifig,'center');
set(handles.text3,'Visible','Off'); % hiding text

delete(handles.axes2.Children(1:length(handles.axes2.Children))) % deleting lines from the current cycle
hold(handles.axes1,'on')
plot(5,5,'LineWidth',2,'MarkerSize',10,'Marker','x','MarkerFaceColor','k','MarkerEdgeColor','k','visible','on','Parent',handles.axes1);
% --------------------------- WP (Wash out Period) = MP ----------------
start_WP(C_cnt) = toc(etime);

guidata(hObject, handles);
stop(handles.WP_tmr);
start(handles.RP_tmr);


function RP_fcn(obj,event,hObject,eventdata)
handles = guidata(hObject);
global C_cnt etime stop_WP start_RP
% global stop_WP etime start_RP
% disp(C_cnt)
pos = get(handles.guifig, 'Position'); % recenter mouse cursor
set(groot,'PointerLocation',[pos(1)+pos(3)/2, pos(2)+pos(4)/2])
set(handles.MP_tmr,'StartDelay', 5.06)
delete(handles.axes1.Children(1:length(handles.axes1.Children)))
movegui(handles.guifig,'center');
set(handles.guifig,'Pointer', 'hand','PointerShapeCData',handles.pshap)
stop_WP(C_cnt) = toc(etime);
start_RP(C_cnt) = toc(etime);
shuff_indx = circshift([1:7]',[rem(C_cnt-1,7) 0]);
pnts=handles.UserData.predpoints;
handles.x = pnts.pnts(C_cnt).iteration(1,:);
handles.y = pnts.pnts(C_cnt).iteration(2,:);
for i = 1:6
    hold(handles.axes1,'on')
    plot(handles.x(i),handles.y(i),'LineWidth',1,'Marker',handles.markers{shuff_indx(i)},'MarkerFaceColor','k','MarkerEdgeColor','k','visible','on', 'Parent', handles.axes1);
end
% Both below ones are necessary
set(handles.axes1,'hittest','on') % enables click detection when graphs are shown on axes1
set(handles.axes1.Children,'hittest','off') % enables click detection when graphs are shown on axes1
guidata(hObject, handles);
stop(handles.RP_tmr);
start(handles.MP_tmr);



function KSS_disp(obj,event,hObject,eventdata)
handles = guidata(hObject);
global procXcuted
procXcuted = 0;
% Present the task even if it is the training
handles.UserData.tskrning = 0;
set(handles.axes1,'visible','off');
set(handles.axes2,'visible','off');
set(handles.KSS,'visible','on');
guidata(hObject, handles);
stop(handles.KSS_tmr);
start(handles.KSSread_tmr);



function KSS_read(obj,event,hObject,eventdata)
handles = guidata(hObject);
global C_cnt
handles.UserData.no_ratings = handles.UserData.no_ratings + 1;
sel_option = get(handles.KSS,'SelectedObject');
KSSchoice = get(sel_option,'tag');
if isempty(KSSchoice) % if nothing selected insert nan in KSS rating
    handles.UserData.KSSrate(handles.UserData.no_ratings) = nan;
else
    handles.UserData.KSSrate(handles.UserData.no_ratings) = sscanf(KSSchoice(1 + length('radiobutton'):end), '%g');
    set(sel_option, 'Value', 0)
end
set(handles.KSS,'visible','off');
guidata(hObject, handles);
% guidata(hObject, handles);
if handles.UserData.tskmode==0
    %  "Ask for Break"
    if handles.UserData.no_ratings>1
        
        disp(handles.UserData.Bchoice(handles.UserData.no_ratings-1))
        if handles.UserData.Bchoice(handles.UserData.no_ratings-1)==1
            % Condition: if there was at least one cycle with the fatigue label
            % detected in each segment there will be a rest break at the end of the
            % segment except for the last segment (ToT=12)
            
            handles.UserData.fdbck_num = handles.UserData.fdbck_num+1; % feedback numbers
            Fatigue_caut = figure(2);
            set(Fatigue_caut, 'Units', 'pixels', 'Position', [1281+188 1+145  905   717],'MenuBar', 'none','Color',[0 0.6 0]);
            movegui(Fatigue_caut,'center');
            i=1;
            plot(1,1)
            xlim([0 10])
            ylim([0 10])
            set(gca,'color',[0 0.6 0],'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
            while i<25
                str = sprintf('A short break,\n The task continues in %d seconds',25-i);
                h = text(5,5,str,'FontSize',24,'Color','w','HorizontalAlignment','center');
                i = i+1;
                pause(1);
                delete(h)
            end
            close(Fatigue_caut)
        end
        
    end
    
    
end
set(handles.axes1,'visible','on');
set(handles.axes2,'visible','on');


if handles.UserData.no_ratings>1
    C_cnt=C_cnt-1;
end
guidata(hObject, handles);
stop(handles.MP_tmr)
stop(handles.KSSread_tmr)
if handles.UserData.tskmode==0
    
    start(handles.MP_tmr);
else
    start(handles.BRKMSG);
end
%%%%%%%%%%%%%%%%%%%%%

% clear
function BRKMSG_disp(obj,event,hObject,eventdata)
handles = guidata(hObject);
global BRK sgmnt

%  "Ask for Break"
if handles.UserData.no_ratings>1
    if handles.UserData.tskmode==2

        disp(sgmnt)
        if handles.UserData.label(sgmnt)==1 || BRK==1% handles.UserData.Pred_lbl(ToT)>0 && C_cnt/20==ToT && ToT~=12
            % Condition: if there was at least one cycle with the fatigue label
            % detected in each segment there will be a rest break at the end of the
            % segment
            
            handles.UserData.fdbck_num = handles.UserData.fdbck_num+1; % feedback numbers
            Fatigue_caut = figure(2);
            set(Fatigue_caut, 'Units', 'pixels', 'Position', [1281+188 1+145  905   717],'MenuBar', 'none','Color',[0 0.6 0]);
            movegui(Fatigue_caut,'center');
            i=1;
            plot(1,1)
            xlim([0 10])
            ylim([0 10])
            set(gca,'color',[0 0.6 0],'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
            while i<25 % 25 seconds rest (micro-break)
                str = sprintf('A short break,\n The task continues in %d seconds',25-i);
                h = text(5,5,str,'FontSize',24,'Color','w','HorizontalAlignment','center');
                i = i+1;
                pause(1);
                delete(h)
            end
            close(Fatigue_caut)
        end
        %     end
    else
        
        disp(handles.UserData.Bchoice(handles.UserData.no_ratings-1))
        if handles.UserData.Bchoice(handles.UserData.no_ratings-1)==1
            % Condition: if there was at least one cycle with the fatigue label
            % detected in each segment there will be a rest break at the end of the
            % segment
            
            handles.UserData.fdbck_num = handles.UserData.fdbck_num+1; % feedback numbers
            Fatigue_caut = figure(2);
            set(Fatigue_caut, 'Units', 'pixels', 'Position', [1281+188 1+145  905   717],'MenuBar', 'none','Color',[0 0.6 0]);
            movegui(Fatigue_caut,'center');
            i=1;
            plot(1,1)
            xlim([0 10])
            ylim([0 10])
            set(gca,'color',[0 0.6 0],'XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
            while i<25
                str = sprintf('A short break,\n The task continues in %d seconds',25-i);
                h = text(5,5,str,'FontSize',24,'Color','w','HorizontalAlignment','center');
                i = i+1;
                pause(1);
                delete(h)
            end
            close(Fatigue_caut)
        end
        
    end
end
set(handles.axes1,'visible','on');
set(handles.axes2,'visible','on');

guidata(hObject, handles);

start(handles.MP_tmr);
stop(handles.data_rec_tmr);
flushinput(handles.t_2)
start(handles.data_rec_tmr);

stop(handles.BRKMSG)



% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% function axes1_ButtonDownFcn(obj,event,hObject,eventdata)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = guidata(hObject);
% global c_clk etime
% disp('here')

global C_cnt etime c_clk
% if strcmp(get(handles.RP_tmr,'Running'),'on')
handles.pnt = get(handles.axes1, 'CurrentPoint'); % the coordinates of clicking on axes1 (main axes)
pnts=handles.UserData.predpoints;
handles.x = pnts.pnts(C_cnt).iteration(1,:);
handles.y = pnts.pnts(C_cnt).iteration(2,:);
X = handles.pnt(1,1);
Y = handles.pnt(1,2);
if c_clk==1 % checking the first click that is starting point
    if pdist([X Y;handles.x(c_clk) handles.y(c_clk)],'euclidean') < 0.2 % comparing "roughly" the clicked and ploted coordinates
 % changing the marker size to be striking
        shuff_indx = circshift([1:7]',[rem(C_cnt-1,7) 0]);
        plot(handles.x(1),handles.y(1),'LineWidth',1,'Marker',handles.markers{shuff_indx(1)},'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',10,'visible','on', 'Parent', handles.axes1);
        c_clk = c_clk + 1;
        handles.UserData.c_clks_cycl(C_cnt) = handles.UserData.c_clks_cycl(C_cnt) + 1;
        handles.UserData.c_clks_time(C_cnt).ttl_clks(handles.UserData.c_clks_cycl(C_cnt)) = toc(etime);
    else
        handles.UserData.ic_clks_cycl(C_cnt) = handles.UserData.ic_clks_cycl(C_cnt) + 1;
        handles.UserData.ic_clks_time(C_cnt).ttl_clks(handles.UserData.ic_clks_cycl(C_cnt)) = toc(etime);
    end
elseif c_clk <= 5
    if pdist([X Y;handles.x(c_clk) handles.y(c_clk)],'euclidean') < 0.2
        line([handles.x(c_clk-1) handles.x(c_clk)],[handles.y(c_clk-1) handles.y(c_clk)],'LineWidth',1,'Color','k')
        % If they were approximaly equal then plot a line connecting the
        % two point together
        c_clk = c_clk + 1;
        handles.UserData.c_clks_cycl(C_cnt)=handles.UserData.c_clks_cycl(C_cnt)+1;
        handles.UserData.c_clks_time(C_cnt).ttl_clks(handles.UserData.c_clks_cycl(C_cnt)) = toc(etime);
    else
        handles.UserData.ic_clks_cycl(C_cnt)=handles.UserData.ic_clks_cycl(C_cnt)+1;
        handles.UserData.ic_clks_time(C_cnt).ttl_clks(handles.UserData.ic_clks_cycl(C_cnt)) = toc(etime);
    end
end
handles.ttlclks_per_cycle = handles.UserData.ic_clks_cycl(C_cnt)+...
    handles.UserData.c_clks_cycl(C_cnt);
handles.UserData.clkd_pnt{handles.ttlclks_per_cycle,C_cnt} = [X Y];
guidata(hObject, handles);



function guifig_CloseRequestFcn(hObject, eventdata, handles)
% global etime
% timer objects have to be found and then deleted
global C_cnt etime start_MP stop_RP start_WP stop_WP start_RP
handles= guidata(hObject);
handles.UserData.task_ending_time = toc(etime);

SbjID=handles.UserData.SbjID;
OutDr=[pwd '\Sbj' SbjID '\' ];
if exist(OutDr,'dir')~=7
    mkdir(pwd,['\Sbj' SbjID '\']);
end
OutDr=fullfile(pwd,['Sbj' SbjID]);

if handles.UserData.Day==1
    flNm = sprintf('s%s_data_trn_d1.mat',SbjID);
elseif handles.UserData.Day==2
    if handles.UserData.tskmode==0
        flNm = sprintf('s%s_data_trn_d2.mat',SbjID);
    else
        flNm = sprintf('s%s_data_fat_d2.mat',SbjID);
    end
else
    if handles.UserData.tskmode==0
        flNm = sprintf('s%s_data_trn_d3.mat',SbjID);
    else
        flNm = sprintf('s%s_data_fat_d3.mat',SbjID);
    end
end
handles.UserData.Fin_date_time = datestr(now);
handles.UserData.cycle_start_WP(1:length(start_WP)) = start_WP;
handles.UserData.cycle_stop_WP(1:length(stop_WP)) = stop_WP;
handles.UserData.cycle_start_RP(1:length(start_RP)) = start_RP;
handles.UserData.cycle_start_MP(1:length(start_MP)) = start_MP;
handles.UserData.cycle_stop_RP(1:length(stop_RP)) = stop_RP;
handles.UserData.cycle_No = C_cnt;

if handles.UserData.tskmode~=0
    ET7_StopDataFileRecording(handles.t_rzm)
    ET7_CloseDataFile(handles.t_rzm);
    handles.UserData.eye_rec_ending_time = toc(etime);

    fclose(handles.t_2)
    delete(handles.t_2)
    fclose(handles.t_rzm)
    delete(handles.t_rzm)
end
UsrDt = handles.UserData;
save(fullfile(OutDr,flNm),'UsrDt')
TmObj=timerfind;
if ~isempty(TmObj)
    TmObj1=TmObj(1);
    if isvalid(TmObj1)
        stop(TmObj)
        delete(TmObj)
    end
end
delete(gcp('nocreate'))
closereq
clear


% --- Executes on mouse motion over figure - except title and menu.
function guifig_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to guifig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global etime
if strcmp(get(handles.RP_tmr,'Running'),'on')
    c_index = handles.UserData.c_index;
    C = get (gca, 'CurrentPoint');
    handles.UserData.Crsr_xyt(:,c_index) = [C(1,1);C(1,2);toc(etime)];
    
    % disp(handles.UserData.Crsr_xyt(:,c_index))
    % cursor position on screen (x,y coordinates)
    % x - first row
    % y - second row
    % t - third row which is elapsed time from the begining (sec) it has to be
    % subtracted by starting time
    c_index = c_index + 1;
    handles.UserData.c_index = c_index;
end
if strcmp(get(handles.MP_tmr,'Running'),'on')
    b = get(handles.guifig,'selectiontype');
    % ToT = handles.UserData.tot;
    if (strcmpi(b,'alt') && isnan(handles.UserData.BreakReqTime(handles.UserData.no_ratings))) && handles.UserData.tskmode~=2
        % max one time per segment only when sync is off
        handles.UserData.BreakReqTime(handles.UserData.no_ratings) = toc(etime);
        handles.UserData.Bchoice(handles.UserData.no_ratings) = 1;
        disp('Asked for break')
    end
end

guidata(hObject,handles);





function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function KSS_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to KSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in KSS.
function KSS_SelectionChangedFcn(hObject, eventdata, handles)

% hObject    handle to the selected object in KSS
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disp('Questionnaire execution')
guidata(hObject, handles);
