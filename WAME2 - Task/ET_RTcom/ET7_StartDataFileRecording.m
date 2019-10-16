%Subroutine to command ET7 to start recording on a currently opened data 
% file.
% If a data file is not currently opened on ET7, or is already recording, 
% no action will be taken.
% t is a TCPIP object.
% A connection must have been established between t and the ET7 before 
% calling this function.
function ET7_StartDataFileRecording(t)
     StartDataFileRecording_cmd = 1;
    % create an output array with chksum temporarily set to zero
    output_bytes = [76,83,65,32,16,0,0,0,StartDataFileRecording_cmd,0,0,0,239,0,0,0];
    % The first 5 bytes are the ASL ET7 "signature"; the 9th byte is the 
    % StartDataFileRecording_cmd; the 13th byte is the checksum.
    
    % write output array values
    fwrite(t,output_bytes);
    
end