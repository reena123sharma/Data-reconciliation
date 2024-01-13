%%-----------------SOLAR CELL SPECS------%%
Pmp_ref = 49;
Imp_ref = 2.88;
Vmp_ref = 17;
Isc_ref = 3.11;
Voc = 21.8;
a_ref = 1.2;
G_ref = 1000;
Iph_ref = Isc_ref;
Io_ref = Isc_ref*exp(-Voc/a_ref);
Rp_ref = 3244.6;%%calculated from the previous code
Rs_ref = 0.571;%%calculatedfrom the previous code

%%--------------DATA-----------%%
load('ts_30.mat');

%%----GENERATING MEASURED DATA-------------%%
MV_solar = V_solar + 0.5*randn(m,1) ;%% normal noise of standard deviation 0.2
MI_solar = I_solar + 0.2*randn(m,1); %% normal noise of standard deviation 0.1
MP_solar = MV_solar.*MI_solar;

%%----------PARAMETERS FOR STATE ESTIMATION-------%%

                    %--SOLAR

Var_solar = [1.0 , 0; 0, 0.01];
A_solar = zeros(2,1);
Iph = G.*Iph_ref/G_ref; %%Assume the cell temp. = ref cell temp. = 25C
Io = Io_ref;
a = a_ref;
Rs = Rs_ref;%%assumption
Rp = Rp_ref*G_ref./G;

tic();
%% 
for i=1:m

               %%----SOLAR------%%
y = [MI_solar(i);MV_solar(i)]  ;             
y_cap = y;
 for j=1:100
  I = y_cap(1);   
  V = y_cap(2);
  A_solar = [ -(Io*Rs/a)*exp((V+I*Rs)/a)-Rs/Rp(i)-1 ,  -(Io/a)*exp((V + I*Rs)/a)-1/Rp(i) ];
  f = Iph(i) - Io*(exp((V+I*Rs)/a) -1) -(V+I*Rs)/Rp(i) - I;
  b = A_solar*y_cap - f;
  y_cap = y - Var_solar*(A_solar')*(inv(A_solar*Var_solar*A_solar'))*(A_solar*y - b);
 end

EP_solar(i) = y_cap(1)*y_cap(2);
EI_solar(i) = y_cap(1);
EV_solar(i) = y_cap(2);

end
RMSE_PEsolar = (immse(P_solar,EP_solar'))^0.5
RMSE_PMsolar = (immse(P_solar,MP_solar))^0.5
RMSE_IMsolar = (immse(I_solar,MI_solar))^0.5
RMSE_IEsolar = (immse(I_solar,EI_solar'))^0.5
RMSE_VMsolar = (immse(V_solar,MV_solar))^0.5
RMSE_VEsolar = (immse(V_solar,EV_solar'))^0.5