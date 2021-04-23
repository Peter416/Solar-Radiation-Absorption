%% ALTITUDE VS TRANSMITTED SOLAR FLUX
% This code was written by Peter Gemayel on Tuesday May 15, 2018
% This script shows: 1) A graph of CO2 Number Density vs Altitude
%2) A graph of CO2 Cross sections vs Frequency
%3) A graph of Solar Flux vs Altitude
%4) A graph of the Surface Solar Flux vs Frequency
clear all
R = 8.314;      %Gas Constant in J/(mol.K)
T = 288.7;      %Average Earth temperature in Kelvin
Ma = 29e-3;     %Average molar mass of air in kg
g = 9.81;       %Gravitational acceleration in m/s^2
H=R*T/(Ma*g);   %Characteristic Atmospheric Height in meters
%--------------------------------------------------------------------------
%Set boundaries for altitude
ztop = 100000;              %top of atmosphere in meters
zsurface = 0;               %Surface height of atmosphere in meters
Nalt = 1000;                %Number of bins for integral of # density
dz=(ztop - zsurface)/Nalt;  %Height step size
z = linspace(1e5,0,Nalt);   %Height array in units of 100 m
 
%--------------------------------------------------------------------------
%Calculate CO2 number density vs height
%CO2 concentration is 400 ppm in 2017
CO2surf=400;                        %CO2 surface concentration in ppm
co2surface=(CO2surf/400)*1e16;      %CO2 density Earth surface molecules/cm^3
nco2_400 = co2surface*exp(-z./H);   %Number density function over altitudes
%Plot the Density Function vs altitude
figure(1)
plot(nco2_400,z)
set(findall(gca, 'Type', 'Line'),'LineWidth',2)
xlabel('Number Density (molecules/cm^3)')
ylabel('Altitude (km)')
yticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
title('Altitude VS Number Density of CO2')
legend('400 ppm')
legend boxoff
%--------------------------------------------------------------------------
%Calculate frequency dependence of CO2 light absorption cross section
%First consider CO2 bending mode
v1=linspace(500,850,100);               %Frequency interval for CO2 bending
sigmab = 1.26e-19;                      %Cross-Section in cm^2
lambdab = 0.0812;                       %Wavelength in cm
vb = 669.2;                             %Peak requency for CO2 bending in cm^-1
sig1=sigmab*exp(-lambdab.*abs(v1 - vb));%Cross-section in cm^2
%Next consider CO2 asymmetric stretching mode
v2=linspace(2300,2340,1000);            %Interval for first line
v3=linspace(2340,2360,1000);            %Interval for second line
v4=linspace(2360,2400,1000);            %Interval for third line
sig2=v2*1.3333e-19 -3.079923e-16;       %Equation for first line in cm^2
sig3(1,1:1000)=4e-18;                   %Equation for second line in cm^2
sig4=-v4*2.2222e-19+5.2843916e-16;      %Equation for third line in cm^2
vtot=[v1 v2 v3 v4];                     %Create array for all frequencies
sigtot=[sig1 sig2 sig3 sig4];           %Create array for all cross sections
%Plot the cross section of CO2 vs frequency
figure(2)
plot(vtot,sigtot)
grid on
xlabel('Frequency v (cm^{-1})')
xlim([100 2700])
ylabel('Cross section cm^2')
ylim([0 4.5e-18])
set(findall(gca, 'Type', 'Line'),'LineWidth',2)
text(600,2.5e-19,'CO_2 bending','fontSize',8)
text(1900,4.2e-18,'CO_2 asymmetric stretching','fontsize',8)
title('CO_2 cross section dependence on frequency')
%--------------------------------------------------------------------------
%Find solar flux at earth's surface
h=6.62607004e-27;       %Planck's Constant in erg*sec
K=1.38064852e-16;       %Boltzmann constant in erg/K
c=2.99792458e10;        %Speed of light cm/s
TSOL=6000;              %Temperature of the Sun in Kelvin
fred=pi*(288.7/TSOL)^4; %Reduction factor
nus1=0;                 %Start frequency
nus2=100000;            %Stop frequency
Delnu=1;                %Frequency step
nsloop=(nus2-nus1)/Delnu;   %Frequency intervals
for i=1:Nalt
    TOT1(i)=0;  %Initialize total intensity for 0 ppm
    TOT2(i)=0;  %Initialize total intensity for 200 ppm
    TOT3(i)=0;  %Initialize total intensity for 400 ppm
    TOT4(i)=0;  %Initialize total intensity for 800 ppm
end
%Loop for integrating over all frequencies and altitudes
for i1 = 1:nsloop
    nus=nus1+Delnu*i1;  %Frequency value for every iteration
    SURFNU(i1)=nus;     %Create array for all frequencies
    Sigma=0;            %Initialize cross section
    if nus>=500 && nus<=850
        Sigma=sigmab*exp(-lambdab*abs(nus-vb));
    elseif nus>=2310 && nus<2340
        Sigma=nus*1.3333e-19 -3.079923e-16;
        elseif nus>=2340 && nus<2360
            Sigma=4e-18;
            elseif nus>=2360 && nus<=2378
                Sigma=-nus*2.2222e-19 +5.2843916e-16;
    end
    %calculate incident solar radiation at frequency nus
    ISOL= fred*(2*h*nus^3*c^2)/(exp((h*nus*c)/(K*TSOL))-1)/1000;%W/m^2/sr/cm^-1
    %Calculate solar intensity at frequency nus vs. altitude
    INU0(1)=ISOL;   %Initialize the intensity for 0 ppm
    INU200(1)=ISOL; %Initialize the intensity for 200 ppm
    INU400(1)=ISOL; %Initialize the intensity for 400 ppm
    INU800(1)=ISOL; %Initialize the intensity for 800 ppm
    for i2 = 2:Nalt
        i3 = i2 - 1;
        INU0(i2) = ISOL;
        INU200(i2) = INU200(i3)*exp(-0.5*Sigma*nco2_400(i2)*dz);
        INU400(i2) = INU400(i3)*exp(-Sigma*nco2_400(i2)*dz);
        INU800(i2) = INU800(i3)*exp(-2*Sigma*nco2_400(i2)*dz);
    end
    SURFINT400(i1)=INU400(1000); %Store surface solar flux for 400 ppm
    %Integrate solar flux over all intensities
    for i4=1:Nalt
        TOT1(i4)=TOT1(i4)+INU0(i4)*Delnu;
        TOT2(i4)=TOT2(i4)+INU200(i4)*Delnu;
    end
end
TOT3(i4)=TOT3(i4)+INU400(i4)*Delnu;
TOT4(i4)=TOT4(i4)+INU800(i4)*Delnu;
%Plot the solar flux vs altitude
figure(3)
plot(TOT1,z,TOT2,z,TOT3,z,TOT4,z)
set(findall(gca, 'Type', 'Line'),'LineWidth',2)
xlabel('Solar Flux (W/m^2)')
ylabel('Altitude (km)')
yticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
title('Altitude VS Solar Flux')
legend('0 ppm','200 ppm','400 ppm','800 ppm','Location','northwest')
legend boxoff
%Plot the surface flux vs frequency
figure(4)
plot(SURFNU,SURFINT400)
xlim([0 5e4])
xlabel('Frequency v (cm^{-1})')
ylabel('Surface Flux per wavenumber(W/m^2/cm^{-1})')
title('Surface Solar Flux')
set(findall(gca,'Type','Line'),'LineWidth',2)
%Find the solar flux at earth's surface for CO2 = 0, 200, 400, 800 ppm
MIN0=min(TOT1);
MIN200=min(TOT2);
MIN400=min(TOT3);
MIN800=min(TOT4);
Table1 = table(MIN0,MIN200,MIN400,MIN800)
