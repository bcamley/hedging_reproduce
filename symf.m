function fout = symf(f,Ntypes)

if(rem(Ntypes,2)==0)
    % even N
    %if(length(logKh)==Ntypes/2)
    fout = [abs(f(1:Ntypes/2)) fliplr(abs(f(1:Ntypes/2)))];
else
   fout = [abs(f(1:((Ntypes-1)/2))) abs(f((Ntypes+1)/2)) fliplr(abs(f(1:((Ntypes-1)/2))))];
end
    