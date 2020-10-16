%Bitstream class:
%
%-Constructor: BitStream(filename)
%   Initialize bitstream from a file if a filename is specified.
%   Otherwise, constructs an empty bitstream.
%
%-Methods add_XXX, addFromFile -> Write data to the bitstream.
%   -obj.add_double(val)      from numeric value of type double.
%   -obj.add_single(val)      from numeric value of type single.
%   -obj.add_signedN(val,N)   from signed N bit integer:
%                               val type must have an integer type, 2<=N<=54).
%                             Limitation to 54-bit to ensure the integer can be represented with double type without loss.
%   -obj.add_unsignedN(val,N) from unsigned N bit integer (val type must have an integer type, 1<=N<=54).
%                             Limitation to 54-bit to ensure the integer can be represented with double type without loss.
%   -obj.add_bin(val)         from a single bin (true, false, 0, or 1).
%   -obj.addFromFile(fid,numBytes): from a file specified by identifier fid.
%                                   The maximum number of bytes to read from the file can be specified as numBytes (default -> read all the file).
%
%-Methods read_XXX -> Read data from the bitstream and advance the reading position accordingly.
%   -obj.read_double()      for numeric value of type double (-> return type = double).
%   -obj.read_single()      for numeric value of type single (-> return type = single).
%   -obj.read_signedN(N)    for signed N bit integer (-> return type = double):   2<=N<=54. Limitation to 54-bit to ensure the integer can be represented with double type without loss.
%   -obj.read_unsignedN(N)  for unsigned N bit integer (-> return type = double): 1<=N<=54. Limitation to 54-bit to ensure the integer can be represented with double type without loss.
%   -obj.read_bin()         for a single bin (-> return type = double with value 0 or 1).
%
%-Method writeToFile -> Write the bitstream in a file specifed by its identifier (without changing bitstream data).
%
%-Method isempty -> Returns true if bitstream is empty.
%
%-Method numbins -> %Returns the number of bins in the bitstream after the reading position.
%
%-Methods addedCount, readCount -> Returns the number bins added to/read from the bitstream since the last call to resetAddedCount/resetReadCount.
%
%-Methods resetAddedCount, resetReadCount -> %Reset the read/add count.
%

classdef BitStream < handle
    properties(Access=private)
        bins=''   % binary data
        bin_pos=1 % reading position
        
        added_num_lastReset=0;
        read_num_lastReset=1;
    end
    
    
    methods
        function obj = BitStream(filename)
            if(exist('filename','var'))
                fid = fopen(filename);
                try
                    obj.addFromFile(fid);
                %close file even if an error is raised during execution.
                    fclose(fid);
                catch exception
                    fclose(fid);
                    rethrow(exception);
                end
            end
        end
        
        %Returns true if bitstream is empty.
        function b = isempty(obj)
            b = obj.bin_pos == length(obj.bins)+1;
        end
        
        %Returns the number of bins in the bitstream after the reading position.
        function n = numbins(obj)
            n = length(obj.bins)+1- obj.bin_pos;
        end
        
        %Returns the number of bins added to the bitstream since the last call to resetAddedCount.
        function n = addedCount(obj)
            n = length(obj.bins) - obj.added_num_lastReset;
        end
        %Reset the added count.
        function resetAddedCount(obj)
            obj.added_num_lastReset = length(obj.bins);
        end
        
        %Returns the number bins read from the bitstream since the last call to resetReadCount.
        function n = readCount(obj)
            n = obj.bin_pos - obj.read_num_lastReset;
        end
        %Reset the read count.
        function resetReadCount(obj)
            obj.read_num_lastReset = obj.bin_pos;
        end
        
        %Write bitstream to a file specified by its file identifier (fid):
        %-Does not modify reading position.
        function writeToFile(obj, fid)
            bins_ = obj.bins(obj.bin_pos:end);
            nBits = length(bins_);
            nBytes = ceil(nBits/8);
            bins_(end+1:end+nBytes*8-nBits) = '0'; %add trailing zeros to have multiple of 8 number of bits.
            bytes = bin2dec( reshape(bins_,8,nBytes)' );
            
            fwrite(fid, bytes ,'uint8');
        end
        
        %Add file data to bitstream specified by its file identifier (fid):
        %-Optional input parameter numBytes can be specified to read only the
        % first bits 'numBytes' bytes of the file (default: numBytes=inf).
        %-Does not modify reading position.
        function addFromFile(obj, fid, numBytes)
            if(~exist('numBins','var') || isempty(numBytes))
                numBytes=inf;
            end
            bytes = fread(fid, numBytes, '*uint8');
            obj.bins = [obj.bins, reshape(dec2bin(bytes,8)',1,[])];
        end

%Methods for adding/reading numeric values to/from bitstream.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Single precision floating point value
        function add_single(obj,val)
            if(~isequal(class(val),'single')), error('expected datatype ''single''');end
            obj.bins = [obj.bins, dec2bin(typecast(val,'uint32'),32)];
        end
        function val = read_single(obj)
            val = typecast(uint32(bin2dec(obj.bins(obj.bin_pos:obj.bin_pos+31))),'single');
            obj.bin_pos = obj.bin_pos+32;
        end
%%%%%%%%
        %Double precision floating point value
        function add_double(obj,val)
            if(~isequal(class(val),'double')), error('expected datatype ''double''');end
            obj.bins = [obj.bins, BitStream.dec2bin64(typecast(val,'uint64'))];
        end
        function val = read_double(obj)
            val = typecast(uint64(BitStream.bin2dec64(obj.bins(obj.bin_pos:obj.bin_pos+63))),'double');
            obj.bin_pos = obj.bin_pos+64;
        end
%%%%%%%%
        %Signed integers with bitdepth N (2<=N<=54) (including sign bit).
        function add_signedN(obj,val,N)
            if(N < 2 || N > 54), error('N must be between 2 and 54');end
            if( ~isinteger(val) )
                error('expected integer type');
            end
            if( val > bitshift(1,N-1)-1 || val < -bitshift(1,N-1))
                error(['expected value between ' num2str(bitshift(1,N-1),'-%u') ' and ' num2str(uint64(bitshift(1,N-1))-1,'%u') '.']);
            end
            
            val_sgn = val < 0;
            if(val_sgn)
                obj.bins = [obj.bins, dec2bin(val_sgn), dec2bin(-(val+1),N-1)];
            else
                obj.bins = [obj.bins, dec2bin(val_sgn), dec2bin(val,N-1)];
            end
        end
        
        function val = read_signedN(obj,N)
            if(N < 2 || N > 54), error('N must be between 2 and 54');end
            
            val_sgn = bin2dec(obj.bins(obj.bin_pos));
            val_abs = bin2dec(obj.bins(obj.bin_pos+1:obj.bin_pos+N-1));
            if(val_sgn)
                val = -val_abs-1;
            else
                val = val_abs;
            end            
            obj.bin_pos = obj.bin_pos+N;
        end
%%%%%%%%
        %Unsigned integers with bitdepth N (1<=N<=54) (including sign bit).
        function add_unsignedN(obj,val,N)
            if(N < 1 || N > 54), error('N must be between 1 and 54');end
            if( ~isinteger(val) )
                error('expected integer type');
            end
            if( val > bitshift(1,N)-1 || val < 0)
                error(['expected value between 0 and ' num2str(uint64(bitshift(1,N))-1,'%u') '.']);
            end
            obj.bins = [obj.bins, dec2bin(val,N)];
        end
        
        function val = read_unsignedN(obj,N)
            if(N < 1 || N > 54), error('N must be between 1 and 54');end
            val = bin2dec(obj.bins(obj.bin_pos:obj.bin_pos+N-1));
            obj.bin_pos = obj.bin_pos+N;
        end
%%%%%%%%
        %Single bin.
        function add_bin(obj,val)
            if(val~=1 && val~=0)
                error('expected binary value (0,1 or true,false)');
            end
            obj.bins = [obj.bins, dec2bin(val)];
        end
        
        function val = read_bin(obj)
            val = bin2dec(obj.bins(obj.bin_pos));
            obj.bin_pos = obj.bin_pos+1;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    methods (Static)
        %Helper functions to overcome limitation of matlab's dec2bin/bin2dec functions for 64 bit integers.
        function bins = dec2bin64(val)
            if(~isequal(val,uint64(val))), error('val must be a non-negative integer smaller than 2^64.');end
            ab_32 = typecast(uint64(val),'uint32');
            bins = [dec2bin(ab_32(2),32), dec2bin(ab_32(1),32)];
        end
        
        function val = bin2dec64(bins)
            if(length(bins)~=64),error('bins must be a char array of 64 elements');end
            a_32 = uint32(bin2dec(bins(1:32)));
            b_32 = uint32(bin2dec(bins(33:64)));
            val = typecast([b_32, a_32],'uint64');
        end
    end
    
end