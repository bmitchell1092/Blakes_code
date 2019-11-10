function CSD = calcCSD(EVP)

totchan = size(EVP,2);
if ndims(EVP) == 3
    maxtr = size(EVP,3);
else
    maxtr = 1;
end

%%% SET CSD PARAMS %%%
% electrode contact spacing in mm:
el_pos = [0.1:0.1:totchan/10];

N = length(el_pos); d = mean(diff(el_pos));
for i=1:N-2
    for j=1:N
        if (i == j-1)
            out(i,j) = -2/d^2;
        elseif (abs(i-j+1) == 1)
            out(i,j) = 1/d^2;
        else
            out(i,j) = 0;
        end
    end
end
%%%%%% END of CSD Params %%%%%%


for tr=1:maxtr
    CSD(:,:,tr)=(-1*out*squeeze(EVP(:,:,tr))');
end

