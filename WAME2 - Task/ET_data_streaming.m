clear
close all hidden
clc

addpath(genpath('ET_RTcom'))
% ExclDtItm=fullfile(CodDr,'DataItems2.xlsx');
[~,DtItmsTx,~]=xlsread('DataItems2.xlsx');
SFct=DtItmsTx(2:end,5); SclFct=str2double(SFct); % deriving the scale factor
t_rzm = tcpip('169.254.181.51', 51000);
set(t_rzm,'ByteOrder','littleEndian');
t_rzm.InputBufferSize = 100000;
MsgBlkSize=86;
ET_data = nan(1000000,9); % where the eye data will be stored

fopen(t_rzm);
pause(1);
SOCKET_TYPE_SDATA_TCP  = 3; % Send Data to Remote via TCP / IP

ET7_SetConnectType(t_rzm, SOCKET_TYPE_SDATA_TCP);

t_2 = tcpip('169.254.181.51', 51000); % make a new tcpip object
set(t_2,'ByteOrder','littleEndian'); % set tcpip specifications
t_2.InputBufferSize = 100000; % set tcpip specifications
% guidata(hObject,handles);
fopen(t_2); % open tcpip port
pause(1)
ET7_OpenDataFile(t_rzm); % opens a new data file on ET7 PC
    ET7_StartDataFileRecording(t_rzm)
    etime = tic;
    ET_rec_start_time = toc(etime);
manualpause = 0;LftDt = [];
while toc(etime)<20
if t_2.BytesAvailable>0
    data1=[LftDt;fread(t_2, t_2.BytesAvailable)];
    [LftDt,OutVrb]=ParsDt(data1,MsgBlkSize,DtItmsTx,SclFct);
    [SmpNm,~]=size(OutVrb);
    ET_data(smpl_cnt+1:smpl_cnt+SmpNm,:)=OutVrb(:,[7 8 9 10 11 12 13 14 15]); % with head-tracker
    smpl_cnt=smpl_cnt+SmpNm;
end

% acm_cyc_idx = cyc_str_smpl(C_cnt-20):cyc_stp_smpl(C_cnt-1);
EH_gaze_hcoord = ET_data(:,5);
EH_gaze_vcoord = ET_data(:,6);
EH_gaze_length = ET_data(:,4);
PplDiam = ET_data(:,1);
EH_gd_x = ET_data(:,7);
EH_gd_y = ET_data(:,8);
EH_gd_z = ET_data(:,9);
cr_diam = ET_data(:,2); % corneal diameter
scen_num = ET_data(:,3); % scene number
% LftDt = [];
end


ET7_StopDataFileRecording(t_rzm)
ET7_CloseDataFile(t_rzm);
eye_rec_ending_time = toc(etime);
fclose(t_2)
delete(t_2)
fclose(t_rzm)
delete(t_rzm)