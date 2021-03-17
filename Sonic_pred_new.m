% Dry Rock Approximation of DEM model (Keys and Xu, 2002)
close all; clc; clear variables;
%--------------------------------------------------------------------------
  v1=input('Enter the name of file in .xlsx:','s');
 %Definition of Constants
  k_m1 = 75*10^9; %Bulk Modulus of carbonate matrix
  k_m2 = 20*10^9; %Bulk Modulus of shale
  mu_m1 = 40*10^9; %Shear Modulus of carbonate matrix
  mu_m2 = 10*10^9; %Shear Modulus of shale
  k_fl = 2.25*10^9; %Bulk Modulus of Water
  v=input('Enter the value of poissons ratio:'); %Pore aspect r
  for i=1:2
   pr(i)=input('Enter the value of pore aspect ratio:');
  end   
%--------------------------------------------------------------------------
 %Initialising K-T parameters for dry rock
  A = -1; B=0;
  R = (1-2*v)/(2-2*v);
   for i=1:2
    p1(i) = (pr(i)/(1-(pr(i)).^2)^1.5)*(acos(pr(i))-(pr(i)*(1-(pr(i)^2))^0.5));
     g(i)=(pr(i)^2/(1-(pr(i)^2)))*(3*p1(i)-2);
   end
%    disp(p1);
%    disp(g);
   
  
  %Defining K-T Scalars (related to p and q)
  for i=1:2
      F1(i)= 1+A*(1.5*(g(i)+p1(i))-R*(1.5*g(i)+2.5*p1(i)-(4/3)));
      
      F2(i)= 1+A*(1+(1.5*(g(i)+p1(i)))-(R/2)*(3*g(i)+5*p1(i)))+B*(3-4*R)+(A/2)*(A+3*B)*(3-4*R)*(g(i)+p1(i)-R*(g(i)-p1(i)+2*(p1(i).^2)));
      
      F3(i)= 1+(A/2)*(R*(2-p1(i))+((1+pr(i).^2)/pr(i).^2)*g(i)*(R-1));
      
      F4(i)= 1+(A/4)*(3*p1(i)+g(i)-R*(g(i)-p1(i)));
      
      F5(i)= A*(R*(g(i)+p1(i)-(4/3))-g(i))+B*p1(i)*(3-4*R);
      
      F6(i)= 1+A*(1+g(i)-R*(g(i)+p1(i)))+B*(1-p1(i))*(3-4*R);
      
      F7(i)= 2+(A/4)*(9*p1(i)+3*g(i)-R*(5*p1(i)+3*g(i)))+B*p1(i)*(3-4*R);
      
      F8(i)= A*(1-2*R+(g(i)/2)*(R-1)+(p1(i)/2)*(5*R-3))+B*(1-p1(i))*(3-4*R);
      
      F9(i)= A*(g(i)*(R-1)-R*p1(i))+B*p1(i)*(3-4*R);
  end
%--------------------------------------------------------------------------
 %Defining set of coefficients p and q
    for i=1:2
        T(i)= 3*F1(i)/F2(i);
        F(i)= (2/F3(i))+(1/F4(i))+(F4(i)*F5(i)+F6(i)*F7(i)-F8(i)*F9(i))/(F2(i)*F4(i));
    end
%    disp(T);
%    disp(F);

    data=xlsread(v1);
  
    GR=data(:,2);
    v_sa=data(:,9);
    v_sh=data(:,3);
    roh=1000*data(:,7); %mineral density
    DTP_obs=data(:,4);
    DTS_obs=data(:,5);
    vp=(1000000/3.28)*(1./DTP_obs);
    vs=(1000000/3.28)*(1./DTS_obs);
    
    for i=1:length(data)
        p(i)=(1/3)*(v_sa(i)*T(1)+v_sh(i)*T(2));
        q(i)=(1/5)*(v_sa(i)*F(1)+v_sh(i)*F(2));
    end
%     disp(p);
%     disp(q);
  %------------------------------------------------------------------------
   %Calculating Bulk and Shear Moduli of rock frame
   phi=data(:,8);
   %roh=1000*data(:,7);
   roh=2710*ones(length(data),1);
   
   for i =1:length(data)
     Mv_k(i)=v_sa(i)*k_m1+v_sh(i)*k_m2;               %Voigt Average of bulk modulus
     Mv_mu(i)=v_sa(i)*mu_m1+v_sh(i)*mu_m2;
     
     Mr_k(i)= k_m1*k_m2/(v_sa(i)*k_m2+v_sh(i)*k_m1);  %Reuss Average of bulk modulus
     Mr_mu(i)= mu_m1*mu_m2/(v_sa(i)*mu_m2+v_sh(i)*mu_m1);
     
     Mh_k(i) =(Mv_k(i) + Mr_k(i))/2;           %Voigt-Reuss-Hill Averageof bulk modulus
     Mh_mu(i)=(Mv_mu(i)+Mr_mu(i))/2;
   end
%    disp(Mh_k);
%    disp(Mh_mu);
%    
   for i=1:length(data)
       k_dry(i)=Mh_k(i).*((1-phi(i)).^p(i));
       mu_dry(i)=Mh_mu(i).*((1-phi(i)).^q(i));
       
       %Gassmann Isotropic Model for calculating saturated rock moduli
       k_sat(i)= k_dry(i) + ((1-(k_dry(i)/Mh_k(i)))^2)/((phi(i)/k_fl)+(1-phi(i))/Mh_k(i)-(k_dry(i)/(Mh_k(i))^2));
       mu_sat=mu_dry;
       
       v_p(i)=3.28*sqrt((k_sat(i)+(4/3)*mu_sat(i))/roh(i));   %1m=3.28ft calibration
       v_s(i)=3.28*sqrt(mu_sat(i)/roh(i));                    %1m=3.28ft
       DTP_pred(i)=10^6*(1/v_p(i));
       DTS_pred(i)=10^6*(1/v_s(i));
   end
%    disp(k_sat);
%    disp(mu_sat);
   
   v_p = (3.28/10)*v_p; v_p=v_p';
   v_s = (3.28/10)*v_s; v_s=v_s';
   
%    DTP_error=(DTP_pred'-DTP_obs)*100/DTP_obs;
%    DTS_error=(DTS_pred'-DTS_obs)*100/DTS_obs;
DTP_error=(DTP_pred'-DTP_obs);
DTS_error=(DTS_pred'-DTS_obs);

   %-----------------------------------------------------------------------
    %Calcultaing relative RMS Error of predicted model
    m =length(data); 
    sum0=0;
    for i=1:m
        s(i) = ((v_p(i) - vp(i))/vp(i)).^2;
        sum0 = sum0+s(i);
        err(i) = (1/2)*sqrt(sum0/m);
%         diff_pred_obsP=(DTP_pred(i)-DTP_obs(i))^2;
%         diff_pred_obsS=(DTS_pred(i)-DTS_obs(i))^2;
%         
    end
%     RMSD_DTP_error=sqrt(sum(diff_pred_obsP)/m)
%     RMSD_DTS_error=sqrt(sum(diff_pred_obsS)/m)
    J = (1/2)*sqrt(sum0/m);
   
   
% m1= mean(DTP_pred,'omitnan');
% m2=mean(DTS_pred);
% %  temp1=DTP_pred'-m1;
% % temp2=DTS_pred'-m2;
%  RMSD_DTP_error= sqrt((sum()^2)/ m)  
%  RMSD_DTS_error= sqrt((sum(temp2)^2)/ m)
   
   
   depth=data(:,1);
   figure
   subplot(1,5,1)
   plot(GR,depth,'r-')
   axis([0 150 min(depth) max(depth)])
   title('GR log')
   xlabel('GR (GAPI)')
   ylabel('depth(ft)')
   set(gca,'Ydir','reverse')
   set(gca,'XAxisLocation','top')
   set(gca,'FontSize',15)
   
   subplot(1,5,2)
   plot(DTP_obs,depth,'r-',DTP_pred,depth,'b-')
   axis([40 240 min(depth) max(depth)])
   title('DEM model')
   xlabel('DTP (us/ft)')
   ylabel('depth(ft)')
   set(gca,'Ydir','reverse')
   set(gca,'Xdir','reverse')
   set(gca, 'XAxisLocation', 'top')
   set(gca,'FontSize',15)
   
   subplot(1,5,3)
   plot(DTP_error,depth,'r-')
   axis([-50 50 min(depth) max(depth)])
   title('DTP Error Log')
   xlabel('Error in (%)')
   ylabel('depth(ft)')
   set(gca,'Ydir','reverse')
   set(gca,'Xdir','reverse')
   set(gca, 'XAxisLocation', 'top')
   set(gca,'FontSize',15)
  
   subplot(1,5,4)
   plot(DTS_obs,depth,'r-',DTS_pred,depth,'b-')
   axis([80 480 min(depth) max(depth)])
   title('DEM model')
   xlabel('DTS (us/ft)')
   ylabel('depth(ft)')
   set(gca,'Ydir','reverse')
   set(gca,'Xdir','reverse')
   set(gca, 'XAxisLocation', 'top')
   set(gca,'FontSize',15)
   
   subplot(1,5,5)
   plot(DTS_error,depth,'r-')
   axis([-50 50 min(depth) max(depth)])
   title('DTS Error log')
   xlabel('Error in (%)')
   ylabel('depth(ft)')
   set(gca,'Ydir','reverse')
   set(gca,'Xdir','reverse')
   set(gca, 'XAxisLocation', 'top')
   set(gca,'FontSize',15)
  %--------------------------------------------------------------------------------
   %Comparison of time-Wylie method with DEM
%    DTP_wylie = data(:,15);
%    Vp_wylie = (1000000./(3.28*DTP_wylie)); 
%   
%    Vs_wylie = sqrt((Vp_wylie.^2)*(((2*0.3) - 1)/((2*0.3) - 2)));
%    DTS_wylie = (1000000./(3.28*Vs_wylie));
%    
%    DTP_error=(DTP_obs-DTP_wylie)*100/DTP_wylie;
%    DTS_error=(DTS_obs-DTS_wylie)*100/DTS_wylie;
%   
%    figure
%    subplot(1,3,1)
%    plot(GR,depth,'r-')
%    xlabel('GR(GAPI)')
%    ylabel('depth(ft)')
%    title('GAMMA RAY LOG')
%    axis([0 150 7245 7920])
%    set(gca, 'Ydir','reverse')
%    set(gca, 'XAxisLocation','top')
%    
%    subplot(1,3,2)
%    plot(DTP_obs,depth,'b-',DTP_wylie,depth,'k-')
%    xlabel('DTP(us/ft)')
%    ylabel('depth(ft)')
%    title('Time-Wylie Method')
%    axis([40 140 7245 7920])
%    set(gca, 'Ydir','reverse')
%    set(gca, 'Xdir', 'reverse')
%    set(gca, 'XAxisLocation','top')
%    
%     
%    
%    subplot(1,3,3)
%    plot(DTS_obs,depth,'b-',DTS_wylie,depth,'k-')
%    xlabel('DTS(us/ft)')
%    ylabel('depth(ft)')
%    title('Time-Wylie Method')
%    axis([80 280 7245 7920])
%    set(gca, 'Ydir','reverse')
%    set(gca, 'Xdir', 'reverse')
%    set(gca, 'XAxisLocation','top')
%    
%   

   
   
   
   
   %-----------------------------------------------------------------------
   