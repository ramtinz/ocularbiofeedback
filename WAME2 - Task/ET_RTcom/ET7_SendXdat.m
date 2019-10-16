%Subroutine to send an "XDAT" value to ET7.
% xdat is an integer between 0 and 65535
% t is a TCPIP object
% A connection must be established between t and the ET7 before calling
% this function

function ET7_SendXdat(t,xdat)

    xdat_command = 5;

    % separate value into most significan byte and least significant byte 
      xdat_hi = bitand(xdat,65280);
      xdat_msb = bitshift(xdat_hi,-8);
      xdat_lsb = bitand(xdat,255);

    % create an output array with chksum temporarily set to zero
    output_bytes = [76,83,65,32,20,0,0,0,xdat_command,0,0,0,0,0,0,0,xdat_lsb,xdat_msb,0,0];
    % The first 5 bytes are the ASL ET7 "signature"; the 9th byte is the 
    % XDAT command; the 17th and 18th bytes are the xdat value lsb and msb.

    % compute check sum
    %  First, sum the bytes in the output array
    sum = int32(0);
    for count = 1:20; 
        incr = int32(output_bytes(count));
        sum = sum + incr;
    end
    %  Get the least significant byte of the sum
    sum16 = uint16(sum);
    sum_lsb = bitand(sum16,255);
    %  Take the "twos compliment" of the lsb
    chksm = bitcmp(sum_lsb,'uint16') + 1;
    %  insert the checksum in the output array  
    output_bytes(13) = chksm; 

    % write output array values
    fwrite(t,output_bytes);
end