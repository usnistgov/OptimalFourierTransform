classdef Flags 
%The flag system is for alerting the user to the possibilities that an MFT or OFT is either:
%  1. Corrupted by round-off error (ususally caused by a frequency being too close to an edge).
%  2. Mathmetically legitimate but maybe physically absurd (cause unknown).    
    properties
        flagString
    end
    
    methods
        function obj = Flags  
            obj = obj.ClearFlag;
         end
        
        function obj = ClearFlag (obj)
            obj.flagString = '-';
        end
        
        function FlagIsSet = CheckFlag(obj,char)
            k = strfind(obj.flagString,char);
            TF = isempty(k);
            FlagIsSet = ~TF;
        end
        
        function obj=AddFlag(obj,char)
            if length(char)~=1 || ~ischar(char)
                error('Flags.addFlag fail, must be a single string character');
            end
            if obj.flagString == '-'
                obj.flagString = '';
            end
            obj.flagString=strcat(obj.flagString,char);
        end
        
        function obj=FlagLowElementsInMatrix(obj, aa, smallChar, largeChar)
            n = length(aa);
            for  i = 1:n
                for j = 1:n
                    x = aa(i,j);
                    if x < 0 
                        x = -x;
                    end
                    if x < 1e-10
                        obj=obj.AddFlag(smallChar);
                    end
                    if x < 1e-12
                        obj = obj.AddFlag(largeChar);
                    end
                end
            end
        end
        
        function obj=FlagHighAmplitudeSinusoids(obj, ts, n, freqs, extent, cosPart, sinPart)
            maxDataValue = -9.9e21;
            for i = n-1:-1:0
                x = abs(ts(i));
                if maxDataValue < x
                    maxDataValue = x;
                end
            end
            loThreshold = 2 * maxDataValue;
            hiThreshold = 100 * maxDataValue;
            
            for i = 1:length(freqs)
               c = cosPart(i);
               s = sinPart(i);
               amp = sqrt(c^2 + s^2);
               if amp > loThreshold
                   nu = extent * freq(i);
                   if nu > 0.0001
                       obj = obj.AddFlag('a');
                        if amp > hiThreshold
                            obj = obj.AddFlag('A');
                        end
                   end
               end
            end     
        end
        
    end
end

