% Sample Matlab script to demonstrate reading real-time data values from
% ASL EyeTrac 7 (ET7) via local network connection.
%
% IMPORTANT: Matlab Instrument Control Toolbox is required
% in order to run this sample
%
% The file "ET7_GetDataItem.m" must be in the same directory with this script, 
% or in the MatLab search path.
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
% When the script is run, the latest horizontal gaze coordinate (Hpos),
% vertical gaze coordinate (Vpos), and pupil diameter value computed by ET7
% should appear on the Matlab command window.

% Note that this sample script shows the simplest means for getting data 
% values.  The data packet returned by ET7 contains a signature, 
% the ID for the command being responded to, the ID of the data value 
% requested, and a check sum.  None of these are checked for errors by this
% simple script; rather, the script simply extracts the expected data value
% from the input buffer.  (C++ and C# examples provided by ASL demonstrate
% full use of all possible error checking features).

% ATTENTION: you must modify the next line to enter the name or ip address
% of the ET7 PC.  Find this address as described above.
remote_pc = '169.254.181.51'

%create  and configure TCP/IP object
t = tcpip(remote_pc, 51000)
set(t,'ByteOrder','littleEndian')

%connect to remote server
fopen(t)

%pause for 1 sec
pause(1)

DataID_Horz_gaze_coord = 13;
DataID_Vert_gaze_coord = 14;
DataID_Pupil_diam = 8;

% Get latest horizontal gaze coordinate value from ET7
    flushinput(t); %clear input buffer
    ET7_GetDataItem(t,DataID_Horz_gaze_coord); %send command to ET7
    InputBuf = fread(t,40); %read 1st 40bytes sent back by ET7
    Hpos = fread(t,1,'single') %read floating point value starting at byte 41

% Get latest vertical gaze coordinate value from ET7
    flushinput(t); %clear input buffer
    ET7_GetDataItem(t,DataID_Vert_gaze_coord); %send command to ET7
    InputBuf = fread(t,40); %read 1st 40bytes sent back by ET7
    Vpos = fread(t,1,'single') %read floating point value starting at byte 41

% Get latest vertical gaze coordinate value from ET7
    flushinput(t); %clear input buffer
    ET7_GetDataItem(t,DataID_Pupil_diam)%send command to ET7
    InputBuf = fread(t,40); %read 1st 40bytes sent back by ET7
    PD = fread(t,1,'single')%read floating point value starting at byte 41

flushinput(t); %clear input buffer


pause(1)

%close copnnection and clean up
fclose(t)
delete(t)


