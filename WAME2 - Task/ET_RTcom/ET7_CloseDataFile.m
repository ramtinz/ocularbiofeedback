%Subroutine to command ET7 to close the currently opened data file
% If a data file is not currently opened on ET7, no action will be taken.
% t is a TCPIP object
% A connection must have been established between t and the ET7 before calling
% this function
function ET7_CloseDataFile(t)
     CloseDataFile_cmd = 4;
    % create an output array with chksum temporarily set to zero
    output_bytes = [76,83,65,32,16,0,0,0,CloseDataFile_cmd,0,0,0,236,0,0,0];
    % The first 5 bytes are the ASL ET7 "signature"; the 9th byte is the 
    % CloseDataFile_cmd; the 13th byte is the checksum.
    
    % write output array values
    fwrite(t,output_bytes);
    
end