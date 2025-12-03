%% UT25 Mini-Module Simulation

%Cell body parameters
Cp_cell = 732;                      %Cell heat capacity in J/(kg*K)
m_cell = 235/1000;                  %Cell mass in kilograms
C_cell = Cp_cell * m_cell;          %Cell thermal mass in J/K
DCR_20 = 2.90;                      %Cell internal resistance at 20 deg C and 50% SOC in milliOhms
DCR_35 = 1.70;                      %Cell internal resistance at 35 deg C and 50% SOC in milliOhms

% Air Cooling parameters
T_init = 25;                        %Initial cell temperature
T_amb = 25;                         %Ambient temperature
A_sink = 700;                       %Surface area of the heatsink in mm^2
h = 120;                            %Guess for forced convective heat transfer coefficient
htc = h * A_sink * (1/1000)^2;      %Heat transfer coefficient in W/K

% Tab Thermal model assuming the following specific heat capacities
% - 470 J/(kg*K) for steel
% - 900 J/(kg*K) for 6063 aluminum and 1100 aluminum
% - 385 J/(kg*K) for copper
m_m3x10_bolt    = 0.1004;           %Mass of M3x10 bolt used in stackup
m_m3x8_bolt     = 0.0891;           %Mass of M3x8 bolt used in stackup
m_nut           = 0.0713;           %Mass of nut

m_heatsink      = 1.42;             %Mass of 6063 aluminum heatsink used in stackup
m_washer        = 0.483;            %Mass of aluminum washer

m_cathode_tab   = 0.263;            %Mass of cathode tab estimated with Al1100
m_anode_tab     = 0.869;            %Mass of anode tab estimated with C101, ignoring the nickel plating

C_tab = (m_m3x8_bolt + m_m3x10_bolt + 2 * m_nut) * 470 / 1000 + (m_heatsink + m_washer + m_cathode_tab) * 900 / 1000 + (m_anode_tab) * 385 / 1000;

%Finding the tab thermal resistance assuming the following thermal conductivities
% - 222 W/(m*K) for Al1100
% - 391.1 W/(m*K) for C101 copper
l = 15;                                                 %Length of exposed tab to bottom of heatsink
w = 18;                                                 %Tab width
t = 0.2;                                                %Tab thickness
R_Al_tab = l / (222 / 1000 * w * t);                    %Thermal resistance of aluminum tab
R_Cu_tab = l / (391.1 / 1000 * w * t);                  %Thermal resistance of copper tab
R_tab = R_Al_tab * R_Cu_tab / (R_Al_tab + R_Cu_tab);    %As one of each type of tab connects to a thermistor, finding the equivalent thermal resistance



% Initializing temperature parameters
T_cell = T_init;                    %Cell Temperature
T_sens = T_init;                    %Thermistor Temperature


thermal_data = importdata("25-12-02 25 Celsius 0% Fan 3C Analysis.csv");

simulated_cell = zeros(length(thermal_data),5);
time = 0;


for i=1:length(thermal_data.data(:,1))-1
    
    step_length = (thermal_data.data(i+1,1) - thermal_data.data(i,1)) / 1000;
    
    if T_cell > 27.5
        Cell_DCR = DCR_35 / 1000;
    else
        Cell_DCR = DCR_20 / 1000;
    end

    %Heat generated in the cell body
    Q_gen = (thermal_data.data(i,2)) ^ 2 * Cell_DCR * step_length;

    %Heat conducted out of the cell to the tab stackup
    Q_cond = 1 / R_tab * (T_cell - T_sens) * step_length;
    
    %Heat convected out of the tab
    Q_cool = htc * (T_sens - T_amb) * step_length;

    %Heat generated at the tab interface due to contact resistance, not used at the moment
    Q_res = 0;      %(thermal_data.data(i,2)) ^ 2 * 187 * 10^(-6) * step_length;

    T_cell = T_cell + (Q_gen - Q_cond) / C_cell;

    T_sens = T_sens + (Q_cond + Q_res - Q_cool) / C_tab;
    
    simulated_cell(i,1) = Q_gen;
    
    simulated_cell(i,2) = Q_cool;

    simulated_cell(i,3) = T_cell;

    simulated_cell(i,4) = time;

    simulated_cell(i,5) = T_sens;

    time = time + step_length;
end

hold on
plot(simulated_cell(:,4),simulated_cell(:,5))
plot(thermal_data.data(:,7),thermal_data.data(:,6))
hold off

writematrix(simulated_cell,'Output.csv')