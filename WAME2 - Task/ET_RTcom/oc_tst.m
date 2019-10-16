clear all
close all
closereq
clc

t_rzm = tcpip('169.254.181.51', 51000, 'NetworkRole', 'client');
% pause(1)
set(t_rzm,'ByteOrder','littleEndian');
t_rzm.InputBufferSize = 1000000;
fopen(t_rzm);
pause(1)
ET7_OpenDataFile(t_rzm);
% pause(2);
ET7_SendXdat(t_rzm, 0);
% pause(1)

% flushinput(t_rzm)
ET7_SetConnectType(t_rzm, 3);
% pause(1)
ET7_StartDataFileRecording(t_rzm)


pause(10)

ET7_StopDataFileRecording(t_rzm)
% pause(1)
ET7_CloseDataFile(t_rzm);
pause(1)
fclose(t_rzm)
% pause(1)
delete(t_rzm)
clear