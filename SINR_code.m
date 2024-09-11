clc;
clear all;
close all;
% hold on 
 
%% =================== RF Parameters ======================================
alpha =2.5;                       %Path Loss Exponent better to have indoor small
PR = 1;
fcRF=2.1*10^9;
GT=1;
GR=1;
gammaI=(3*10^8)^2*GT*GR/(16*pi^2*fcRF^2); %used in sinr 
 
%% ===================== THz Parameters ===========================
PT = 1;                             % Tranmitted Power
kf = 0.05;                           % Absorbtion loss
fcTH=1.0*10^12;
GTT=316.2;
GRR=316.2;
% GTT=31.62;
% GRR=31.62;

thetabs=pi/6;%%%in degrees
thetamt=pi/6;%%%in degrees
FBS=thetabs/(2*pi);
FMT=thetamt/(2*pi);
prob= FBS*FMT;
gammaII=(3*10^8)^2*GTT*GRR/(16*pi^2*fcTH^2);
 

R_max = 100;
Nit = 10000;
%lambdas=10;

%% ===================== Rate and SINR theshold calculation ========================================
Rate = [5 ]*10^9;
Wt=5*10^8;
Wr=40*10^6;
SINRthRF =2.^(Rate./(Wr))-1; %thresholds
SINRthTH =2.^(Rate./(Wt))-1;

Bias=[1000 100 1 0.001 0.0001 0.0001 0.00001];%%%0.05
%Bias=[ 10^6  10^5 10^4  10^4 10^3 10^3 10^3];%%%0.2


      
      
      for i=1:Nit
        
        UoI = 1;
        lambda=2;
        % The index of User of Interest
        NRF = 2;
        NTHzT = 2;
        Nu = 16;

        
        fadeRand = exprnd(1,NRF,UoI);
        
      
            rue = R_max*sqrt(rand(Nu,1));
            thetau = 2*pi*rand(Nu,1);
            Xu = rue.*cos(thetau);
            Yu = rue.*sin(thetau);
            
            
            rTHzbs = R_max*sqrt(rand(NTHzT,1));
            thetaTb = 2*pi*rand(NTHzT,1);
            Xb = rTHzbs.*cos(thetaTb);
            Yb = rTHzbs.*sin(thetaTb);
            
            rRFzbs = R_max*sqrt(rand(NRF,1));
            thetaRb = 2*pi*rand(NRF,1);
            Xrb = rRFzbs.*cos(thetaRb);
            Yrb = rRFzbs.*sin(thetaRb);
            
            
            % Distacnces from THz BSs
            
            [Xmp_Tmat, Xp_Tmat] = meshgrid(Xu,Xb);
            [Ymp_Tmat, Yp_Tmat] = meshgrid(Yu,Yb);
            D_ue_Tbs = sqrt((Xmp_Tmat-Xp_Tmat).^2 + (Ymp_Tmat-Yp_Tmat).^2);
            
            % Distacnces from RF BSs
            
            [Xmp_Rmat, Xp_Rmat] = meshgrid(Xu,Xrb);
            [Ymp_Rmat, Yp_Rmat] = meshgrid(Yu,Yrb);
            D_ue_Rbs = sqrt((Xmp_Rmat-Xp_Rmat).^2 + (Ymp_Rmat-Yp_Rmat).^2);
            
            fadeRand = exprnd(1,NRF,Nu);
            SRF=gammaI.*fadeRand.*PR.*D_ue_Rbs.^(-alpha); %signal matrix for RF
            NP=(10)^-10;
            interf=repmat(sum(SRF,1),NRF,1)-SRF; %interference for RF
            RPrAllu1 = log2(1+SRF./(NP+interf)); %power from all base-stations to all users
            RPrAllu = log2(1+SRF./NP);
            
            
            
            fadeRand1 = exprnd(1,NTHzT,Nu);
            STHz=gammaII.*fadeRand1.*PT.*exp(-kf.*D_ue_Tbs)./(D_ue_Tbs.^2); %signal matrix for THZ
            interfT=repmat(sum(STHz,1),NTHzT,1)-STHz; %interference matrix for THz
            TPrAllu1 =log2(1+STHz./(NP+interfT));
            TPrAllu =log2(1+STHz./(NP));
            
            SINR_Matrix = [RPrAllu1;TPrAllu1]; 
            
            
            [SINR_max, I] = max(SINR_Matrix, [], 1);
            
            SINR_sum(i) = sum(SINR_max); %benchmark for max SINR
            
                       
                        
            zeroes_matrix = zeros(size(SINR_Matrix));
            
            user_throughput = zeros(size(SINR_Matrix));
            
            %zeroes_matrix() =  1
            
            for j = 1:Nu
               
                
                zeroes_matrix(I(j),j) = 1; 
                user_throughput(I(j),j) = SINR_Matrix(I(j),j);
                
            end
            
            zeroes_sum = sum(zeroes_matrix,2);
            
            users_sum = sum(user_throughput,2);
            
            SINR_div = users_sum./zeroes_sum;
            
            for k = 1:size(SINR_div)
            
                if (isnan(SINR_div(k)))
                    SINR_div(k) = 0;
                end
                
            end
            
            final_user_div_SINR(i) = sum(SINR_div); %benchmark #2 for max sinr with division
            
            
            
            Signal_Matrix = [SRF;STHz]/NP; 
            
            Interference_Matrix = [interf;interfT]; 
            
            Comb_reshaped_s(:,i) = ((reshape(Signal_Matrix,[], 1)));
            
            Comb_reshaped_i(:,i) = ((reshape(Interference_Matrix,[], 1)));
            
            Comb_reshaped(:,i) = ((reshape(SINR_Matrix,[], 1)));

            

   
            
            
            %RateR(i) = Wr*log(1 + RPrAllu);
            %Recevied power & SNIR from THz BSs
            


            
             
            
            
            
            %TDesiredIndex = TUserIndex(UoI);
            %TPrDesired = TPrAllu(TDesiredIndex,UoI);
            %CumPr = sum(TPrAllu(:,UoI));
            %TSINR = TPrDesired/(prob*(CumPr-TPrDesired)); %for terahertz
            %RateT(i) = Wt*log(1 + TPRAllu);
            
              
            
            
      end
        
      avg_max_sinr = mean(SINR_sum)
      
      avg_div_sinr = mean(final_user_div_SINR)
        
   
    


 




semilogx(pika,PCov1,'r-s','LineWidth',1.5,'DisplayName', 'RF')
hold on
semilogx(pika,PCov2,'b^','LineWidth',1.5,'DisplayName', 'THz')
 hold on
semilogx(pika,PCov,'ko--','LineWidth',1.5,'DisplayName', 'Opportunistic RF/THz 1')
 hold on
semilogx(pika,PCov3,'g:>','LineWidth',1.5,'DisplayName', 'Hybrid RF/THz')
hold on
semilogx(pika,PCov4,'m:>','LineWidth',1.5,'DisplayName', 'Biased Opportunistic')
legend


% figure
% plot(pika/(pi*R_max^2),PAss_T,'ro','LineWidth',1.5,'DisplayName', 'RSRP')
% hold on
% plot(pika/(pi*R_max^2),PAss_T1,'s-.k','LineWidth',1.5,'DisplayName', 'Biased')
% hold on
% 
% % %%%%%%%%%%%%%0.05
% % Theory=[0.518523, 0.717766, 0.892221, 0.950681, 0.968294, 0.976587, 0.988338]
% % Theory1=[0.99, 0.987, 0.8922, 0.1677, 0.014, 0.018, 0.0047297]
% % 
% % %%%%%%%%%%%%%0.2
% Theory=[0.098974, 0.185687, 0.388227, 0.598855, 0.720572, 0.795095, 0.915049];
% Theory1=[0.88, 0.949376, 0.988, 0.998, 0.999, 0.9999, 0.999];
% plot(pika/(pi*R_max^2),Theory,':ks','LineWidth',1.5,'DisplayName', 'Theory')
% hold on
% plot(pika/(pi*R_max^2),Theory1,':ks','LineWidth',1.5,'DisplayName', 'Theory1')
% legend

