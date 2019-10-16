% Sample code demonstrating use of Matlab to open and close data files, 
% start and stop data recording, and set XDAT values on ASL EyeTrac 7 (ET7)
% via local network commands.
%
% IMPORTANT: Matlab Instrument Control Toolbox is required
% in order to run this sample
%
% The files "ET7_SendXdat.m", "ET7_OpenDataFile.m", "ET7_CloseDataFile.m", 
% "ET7_StartDataFileRecording.m", and "ET7_StopDataFileRecording.m" must 
% all be in the same directory with this script, or in the MatLab search 
% path.
%
% This sample was tested with Matlab version R2009b.
%
% Before running this program (before attempting connection) run
% the ET7 Application, click the "Configuration" button at the bottom 
% of the Control Table "Configuration" tab, in the group labeled 
% "Real-Time Input/Output", and note the IP address listed
% on the resulting dialog (see ET7 manual for details). Click the "Listen
% button".  
%
% When the scrip is run successfully, a data file, with date and time as the file
% name, will be opened on ET7; XDAT will be set to zero, and then recording
% will begin.  After 5 secondes, XDAT will be set to 100, then 150 and then
% 200, with a 5 second period between each XDAT change.  5 seconds after
% the last XDAT change recoding will pause and the file will close.
%
% ATTENTION: you must modify the next line to enter the name or ip address
% of the ET7 PC.  Find this address as described above.
remote_pc = '192.168.1.35';

%create  and configure TCP/IP object
t = tcpip(remote_pc, 51000)
set(t,'ByteOrder','littleEndian');

%connect to remote server
fopen(t);

%pause for 1 sec
pause(1);

ET7_OpenDataFile(t);
pause(2);
ET7_SendXdat(t, 0);
ET7_StartDataFileRecording(t);
pause(5);
ET7_SendXdat(t, 100);
pause(5);
ET7_SendXdat(t,150);
pause(5);
ET7_SendXdat(t,200);
pause(5);
ET7_StopDataFileRecording(t);
ET7_CloseDataFile(t);

pause(1);

%close copnnection and clean up
fclose(t);
delete(t);


