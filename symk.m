function logKvec = symk(logKh,Ntypes)

if(rem(Ntypes,2)==0)
    % even N
    %if(length(logKh)==Ntypes/2)
    logKvec = [-fliplr(abs(logKh(1:Ntypes/2))) abs(logKh(1:Ntypes/2))];
else
   logKvec = [-fliplr(abs(logKh(1:((Ntypes-1)/2)))) 0 abs(logKh(1:((Ntypes-1)/2)))];
end
    