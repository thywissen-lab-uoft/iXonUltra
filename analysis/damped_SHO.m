function damped_SHO

% Damped Harmonic Oscillator
% The equation of motion for a damped harmonic oscillator is given by
% ma = -kx - bv +F sin(omega*t)
% 
% 
% omega_0 = sqrt(k/m)
%
% The relaxation time is
% The quality factor relates to the coefficients via
% Q = omega_0*tau = omega_0 * m /b

% ==> tau = m/b
% 
% So that the equation of motion becomes
%
% a + x*k/m + v/tau = F_0/m * sin(omega*t)
% a + x*omega0^2 + v/tau = F_0/m * sin(omega*t)
%
% We often express the force in terms of the drive amplitude and the
% harmonic trapping force F_0=XD * omega_D^2 * m
%
% So that the EOM reduces to 
%
% a + x*omega0^2 + v/tau = XD omega_D^2*sin(omega*t)
%
% This leads to coupled differential equations of 
% dot(x) = v
% dot(v) = -x*omega_0^2-v/tau+XD*omegaD^2*sin(omega*t)

% Experimentally we ramp the amplitude XD=XD(t) to smoothly turn on the
% driving force. How does the finite ramp affect the measured signal x(t)?
% 
% Expect in the long time limit for the presence of the ramp on to become
% neglible. But how long is long? And how does it disrupt measurement of
% the lineshape function to extract tau?



amu         = 1.66054e-27; % atomic mass unit [kg]
m           = 40*amu;        % 40K mass


omega0      = 2*pi*57;     % Natural trap Frequency in Hz [1/s]
omegaD      = 2*pi*42;     % Trap frequency of XDT [1/s]
T           = .1;          % Ramp up time     [s]
tau         = 2e-3;        % Relaxation Time [s]
gamma       = 1/tau

x0          = 1*1e-6;       % Drive amplitude [m]


fvec = 1:1:120;

kB=1.380649e-23;  % boltzmman consntat [J/K]

aL=0.532e-6;

hbar = 1.05457182e-34;

sigma0=aL^2/hbar;
% hF=figure(20);
% hF.Color='w';
% clf(hF)

cc=cool(length(fvec));
Svec=zeros(length(fvec),1);
Cvec=zeros(length(fvec),1);

sigmaR=Svec;
sigamI=Cvec;
injected_energy = zeros(length(fvec),1);
for kk=1:length(fvec)
    [t,y]=timeEvolve(2*pi*fvec(kk),tau);

    force_curve = m*XD(t).*omegaD^2.*sin(2*pi*fvec(kk)*t);

    % Power = Force * velocity
    dt=t(2)-t(1);
    power_curve = force_curve.*y(:,2);
    injected_energy(kk) = 1e9*trapz(dt,power_curve)/kB;


    % plot(t,y(:,1),'-','color',cc(kk,:));
    % hold on


    Tdrive = 1/fvec(kk);
    Tinspect = t(end)-2*Tdrive;
    ind = find(t>Tinspect,1);
    tSub = t(ind:end);
    xSub = y(ind:end,1);
    myfit = fittype(@(C,S,t) C*cos(2*pi*fvec(kk)*t)+S*sin(2*pi*fvec(kk)*t),'independent','t',...
        'coefficients',{'C','S'});
    opt=fitoptions(myfit);
    opt.Start=[0 1];
    fout = fit(tSub,xSub,myfit,opt);

    Svec(kk) = fout.S;
    Cvec(kk) = fout.C;

    sigma_R(kk) = -fout.C*2*pi*fvec(kk)./(m*x0.*omegaD^2);
    sigma_I(kk) = -fout.S*2*pi*fvec(kk)./(m*x0.*omegaD^2);

end
% xlim([0 100]);
% 
% colormap(hF,cc);
% colorbar


hF=figure(21);
hF.Color='w';
clf(hF)

subplot(121);
plot(fvec,sigma_R/sigma0);
hold on
plot(fvec,sigma_I/sigma0);
xlabel('frequency (Hz)');
ylabel('conductivity (\sigma_0)')

subplot(122);
plot(fvec,injected_energy,'k-');
ylabel('injected energy (nK)');
xlabel('frequency (Hz)');
%% Functions

function x=XD(t)
    tt = t/T;
    x=x0*((6*tt.^5-15*tt.^4+10*tt.^3).*(tt<=1) + 1*(tt>1));
end

    function [t,y]=timeEvolve(myomega,mytau)
        tVec = 0:1e-4:(T+0.05);
        [t,y] = ode45(@foo,tVec,[0; 0]);

        % [t,y] = ode45(@foo,[0 T+0.05],[0; 0]);

        function dydt = foo(t,y)%
        dydt = [y(2);
            -y(1)*omega0^2-y(2)/mytau+XD(t)*omegaD^2*sin(myomega*t)];
        end
    end




end

