function A=merge(D,direction)
% this function merge multicell matrix to one matrix
%
% direction:direction ==1 ?? direction ==2 -> (default=1)
%
%
% HORIKAWA tomoyasu 09/09/03
% modified by tomoyasu horikawa 090917; direction was inversed
len=numel(D);
if len == 0;
    A = [];
else
    A=D{1};
    if exist('direction','var')==1
        if direction==1
            for i=2:len
                A=[A;D{i}];
            end
        else
            for i=2:len
                A=[A,D{i}];
            end
        end
    else
        for i=2:len
            A=[A;D{i}];
        end
        
    end
end