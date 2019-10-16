%Subroutine to command ET7 to stop recording on a currently opened data 
%file.
% If a data file is not currently opened on ET7, or is not recording, 
% no action will be taken.
% t is a TCPIP object.
% A connection must have been established between t and the ET7 before 
% calling this function.
function ET7_StopDataFileRecording(t)
     StopDataFileRecording_cmd = 2;
    % create an output array with chksum temporarily set to zero
    output_bytes = [76,83,65,32,16,0,0,0,StopDataFileRecording_cmd,0,0,0,238,0,0,0];
    % The first 5 bytes are the ASL ET7 "signature"; the 9th byte is the 
    % StopDataFileRecording_cmd; the 13th byte is the checksum.
    
    % write output array values
    fwrite(t,output_bytes);
    
end