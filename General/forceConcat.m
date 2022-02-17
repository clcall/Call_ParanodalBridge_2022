function [a,b] = forceConcat(a,b,dim)
    if size(a) == size(b)
        return
    end
    if nargin==2 && isvector(a) && isvector(b)
        if size(a,1) > size(a,2)
            flipflag = true;
            a = a';
            b = b';
            dim = 2;
        else
            flipflag = false;
        end
    elseif dim==1
        biggest = max(size(a,dim),size(b,dim));
        if biggest > size(a,dim) % b has more rows than a
            r = size(b,1);
            c = size(a,2);
            temp = NaN(r,c); % since b has more rows than a, make a matrix with the number of rows in b and keep the number of columns from a
            temp(1:size(a,1),1:size(a,2)) = a;
            a = temp;
        else % a has more rows than b
            r = size(a,1);
            c = size(b,2);
            temp = NaN(r,c);
            temp(1:size(b,1),1:size(b,2)) = b;
            b = temp;
        end
        return
    else
        flipflag = false;
    end
    biggest = max(size(a,dim),size(b,dim));
    % Force the sizes of the matrices/vectors to be the same length for concatenation
    if biggest > size(a,dim)
        temp = NaN(size(b));
        temp(1:size(a,1),1:size(a,2)) = a;
        a = temp;
    else
        temp = NaN(size(a));
        temp(1:size(b,1),1:size(b,2)) = b;
        b = temp;
    end
    if flipflag
        a = a';
        b = b';
    end
end
%% ORIGINAL CODE FROM calcProcLosses.m
%     biggest = max(length(ALL.stats.procLosses.ctrl(:,end)) , length(ALL.stats.procLosses.cupr(:,end)));
%     % Force the sizes of cupr/ctrl vectors to same length for concatenation
%     if biggest > length(ALL.stats.procLosses.ctrl(:,end))
%         temp = NaN(size(ALL.stats.procLosses.cupr));
%         temp(1:size(ALL.stats.procLosses.ctrl,1),1:end) = ALL.stats.procLosses.ctrl;
%         ALL.stats.procLosses.ctrl = temp;
%     else
%         temp = NaN(size(ALL.stats.procLosses.ctrl));
%         temp(1:size(ALL.stats.procLosses.cupr,1),1:end) = ALL.stats.procLosses.cupr;
%         ALL.stats.procLosses.cupr = temp;
%     end