% Sample code showing how to send XDAT values to ASL EyeTrac 7 (ET7) via 
% local network connection
%
% IMPORTANT: Matlab Instrument Control Toolbox is required
% in order to run this sample
%
% The file "SendXdat.m" must be in the same directory with this script, 
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

% set xdat to any integer value between 0 and 65535
xdat = 500  

ET7_SendXdat(t,xdat)
  
%pause for 1 sec
pause(1)

%close copnnection and clean up
fclose(t)
delete(t)




