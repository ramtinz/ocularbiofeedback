function [LftDt,OutVrb,InvlvVrbl]=ParsDt(data,MsgBlkSize,DtItmsTx,SclFct)

% Data parsing function 

% each data block consists of some header bytes and a data segment which
% includes the variables exported from the ET7

% LftDt: left data- extra bytes read from the TCP buffer and cannot be fit into one full data block
% OutVrb: Output variables two dimensional (time, [timestamp, involved variable names]) 
%           timestamp: apparently in units of 0.1 micro sec
% InvlvVrbl: the names of the involved variables

% data: all the read  bytes from the TCP buffer
% MsgBlkSize: the number of bytes in the data segment of each block
% DtItmsTx: content of the excel sheet (DataItems2) containing all the
%           information about the data segment content
% SclFct: scale factor


        % Byte 0: ASLSignature (0x2041534C) // ASL Signature: “ ASL”
        % Byte 4: MsgSize // Total message size: 56 + DataSize
        % Byte 8: Cmd (0x00000081) // Message Command: Data Message
        % Byte 12: Checksum (0x00000000) // Message Checksum: Reserved to 0 for Data Message
        % Byte 16: DataSize // Selected Data Size
        % Byte 20: FrameSize(0x00000000) // Video Frame Size: 0 since no video
        % Byte 24: FrameNo // Current Frame Number
        % Byte 28: Reserved // reserved bytes
        % Byte 32: TimeStamp // Time Stamp of this frame
        % Byte 40: UpdateRate // Frame Update Rate
        % Byte 44: Reserved // reserved bytes
        % Byte 48: CheckState // Data Selected State
        % Byte 56: BufStart // Data Buffer

        HdrInfo=[17 20; % data size
 25 28; % frame No
 33 40; % time stamp
 41 44; % update rate
 49 56; % Check state
 ]; % byte number start and stop

DtInx=[57 MsgBlkSize]; % data buffer in each data segment
TtlVrbNm=DtItmsTx(2:end,2); % the variable names (all possible entries)
BlkNm=fix(length(data)/MsgBlkSize); % per time stamp, there exists one data block
LftDt=data(BlkNm*MsgBlkSize+1:end); % if some extra bytes have been read which cannot fit into a full data block
BytCount=DtItmsTx(2:end,4); BytCount=str2double(BytCount); % the number of bytes contained in each item in the data segments
DtTyp=DtItmsTx(2:end,3); % the format of data type for the items in the data segment

for BlkCnt=1:BlkNm
    Inx=(BlkCnt-1)*MsgBlkSize+1:BlkCnt*MsgBlkSize;
    SgBlk=data(Inx); HdrNm=size(HdrInfo,1);
    
    HdrSg=nan(HdrNm,1);
    for HdrCnt=1:HdrNm
        TtlByts=SgBlk(HdrInfo(HdrCnt,1):HdrInfo(HdrCnt,2));
        Y=typecast(uint8(TtlByts),['uint' num2str(8*length(TtlByts))]);
        HdrSg(HdrCnt)=Y;
    end
    TmStmpVl=HdrSg(3); % Time stamp
    ChckSt=HdrSg(end); % check state bite (defining what variables are included in the data segment)
    %%%% to convert checkState into bitwise logical array
    StrChckSt=dec2bin(ChckSt); Lbn=length(StrChckSt);
    AppndStrChckSt=[StrChckSt(:) repmat(' ',Lbn,1)];
    InvlvIndx=logical(str2num(AppndStrChckSt)); % involved variables in the data buffer
    InvlvIndx=InvlvIndx(end:-1:1);
    
    if isequal(BlkCnt,1)
        InvlvVrbl=TtlVrbNm(InvlvIndx); % invovled variables int the data segment
    end
    DtTypSlc=DtTyp(InvlvIndx); BytCountSlc=BytCount(InvlvIndx); SclFctSlc=SclFct(InvlvIndx); % invovled data type, byte count and scale factor
    DtSegSz=length(DtTypSlc); RdDtSeg=0;
    DtSgBlk=SgBlk(57:MsgBlkSize);
    
    if ~exist('OutVrb','var')
        OutVrb=nan(BlkNm,DtSegSz+1); % time, varaible
    end
    
    OutVrb(BlkCnt,1)=TmStmpVl;
    for DtSegCnt=1:DtSegSz % parsing the data segment
        IndDtBlk=RdDtSeg+(1:BytCountSlc(DtSegCnt)); % byte counter
        DtSgBlkByts=uint8(DtSgBlk(IndDtBlk));
        DtFrmt=lower(DtTypSlc{DtSegCnt}); % data format
        DtFrmt(DtFrmt==32)=[]; % remove blank space
        switch DtFrmt
            case 'byte'
                OutDtFrmt='uint8';
            otherwise
                OutDtFrmt=DtFrmt;
        end
        SgVl=typecast(DtSgBlkByts,OutDtFrmt);
        OutVrb(BlkCnt,DtSegCnt+1)=SgVl*SclFctSlc(DtSegCnt);
        RdDtSeg=IndDtBlk(end);
    end
end