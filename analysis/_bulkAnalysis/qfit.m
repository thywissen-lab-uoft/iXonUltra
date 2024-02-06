function sigma = qfit(N,T,G,w)
    
    global Rvalues;
    global energies;
    
    hbar = 6.626e-34/(2*pi);
    kb = 1.380649e-23;
    aL = 532e-9;
    beta = 1/(kb*T);

    ss = 0;
    for loop1 = 1:length(Rvalues)
        for loop2 = 1:length(Rvalues)
            ss = ss + (abs(Rvalues(loop1,loop2))^2)*(exp(-beta*energies(loop1))-exp(-beta*energies(loop2)))/(w - (energies(loop1)-energies(loop2))/hbar + 1j*G);
        end
    end

    sigma = -1j*(N/(aL^2))*w*ss;

end