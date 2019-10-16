%Subroutine to command ET7 to send latest value of specified data item. 

function ET7_GetDataItem(t, DataID)
  GetDataItem_cmd = 25;
  
  output_bytes = [76,83,65,32,20,0,0,0,GetDataItem_cmd,0,0,0,0,0,0,0,DataID,0,0,0];
  
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
    fwrite(t,output_bytes);
end
