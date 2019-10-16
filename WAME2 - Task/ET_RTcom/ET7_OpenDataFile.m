%Subroutine to command ET7 to open a data file
% File name will be current date and time
% t is a TCPIP object
% A connection must be established between t and the ET7 before calling
% this function
function ET7_OpenDataFile(t)
     OpenDataFile_cmd = 3;
    % create an output array with chksum temporarily set to zero
    output_bytes = [76,83,65,32,16,0,0,0,OpenDataFile_cmd,0,0,0,237,0,0,0];
    % The first 5 bytes are the ASL ET7 "signature"; the 9th byte is the 
    % OpenDataFile_cmd; the 13th byte is the checksum.
    
    % write output array values
    fwrite(t,output_bytes)
end

