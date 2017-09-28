clear
clc
close all
addpath Functions
tic

%%

T=0.01 ;
Pu1 = blkdiag( 1,1,1,1) ;
 
runtime =5000;
effTimeArray = 1 : runtime ;
 
Xu1 =[5;1;90;-2.5];
 
ff = [ 1 T ; 0 1 ] ;
fy = blkdiag( ff , ff ) ;
 
Z = zeros(2,1);
%%
Pu2 = blkdiag( 1,1,1,1) ;
Xu2 =[5;1;90;-2.5];
 
%%
sigma1 = 2;
sigma2 = 2;
 
Qxyz1 = sigma1^2 * [T^3/3      T^2/2 ;
                                        T^2/2        T      ] ;
Qxyz2 = sigma2^2 * [T^3/3      T^2/2 ;
                                        T^2/2        T      ] ;
Q = blkdiag( Qxyz1 , Qxyz2 ) ;

%%
% for R
Rhoerror = 1.5 ;

for i= 1: runtime
    Xnoise =0.2*randn(1,1);    Ynoise =0.2*randn(1,1);
    if Xnoise>1
        Xref =2*i*T+1 ;
    elseif  Xnoise<-1      Xref =2*i*T-1 ;
    else
        Xref =2*i*T +Xnoise;
    end
    
     if Ynoise>1
           Yref =100-2*i*T+1 ;
    elseif  Ynoise<-1      Yref =100-2*i*T -1 ;
    else
         Yref =100-2*i*T +Ynoise;
    end
      
    horizion_ref1 = sqrt(Xref^2+Yref^2) ;    
    horizion_ref2 = sqrt((Xref-100)^2+(Yref-100)^2) ;
    
    Z1noise = 1*randn(1,1);    Z2noise = 1*randn(1,1);
    if Z1noise>1  Z(1,1) =horizion_ref1+1 ;
    elseif  Z1noise <-1        Z(1,1) =horizion_ref1-1 ;
    else
        Z(1,1) =horizion_ref1 + Z1noise;
    end
    
     if Z2noise>1  Z(2,1) =horizion_ref2+1 ;
    elseif  Z2noise <-1       Z(2,1) =horizion_ref2-1;
    else
        Z(2,1) =horizion_ref2 + Z2noise;
     end
  
      R =  eye( size( Z ,1 ) ) * Rhoerror^2;
    
    %%
    %KF
    
    Pp1 = fy * Pu1 * fy' + Q ;
    Xp1 = fy * Xu1 ;                  % fy為State transition matrix
    
    horizion11 = sqrt(Xp1(1)^2+Xp1(3)^2) ;
    horizion12 = sqrt( (Xp1(1)-100)^2+(Xp1(3)-100)^2) ;
    horizion_total1 = [horizion11 ; horizion12];
    
    H1=[Xp1(1)/horizion11 ,           0 ,     Xp1(3)/horizion11,                0  ;
        (Xp1(1)-100)/horizion12 , 0 ,      (Xp1(3)-100)/horizion12,       0  ];    
    
    %K=(inv(Pp)+H' * inv(R )* H)\H'*inv(R);
    K1 = Pp1 * H1' / (H1 * Pp1 * H1' + R) ;
    Xu1 = Xp1 +  K1 * ( Z -horizion_total1) ;        %   Z為修正過後的pseudorange
    
    I = eye( size(Xu1,1) ) ;
    Pu1 = (I - K1 * H1) * Pp1;
    
    Xref_run(i,:) = [Xref Yref ];
    Xxyz1(i,:)=Xu1(:,1);
    Zxyz(i,:)=Z(:,1)';
    
    %%
    %REFINEMENT
     Xp2 = fy * Xu2 ;                  % fy為State transition matrix
     Pp2 = fy * Pu2 * fy' + Q ;
    
    horizion21 = sqrt(Xp2(1)^2+Xp2(3)^2) ;
    horizion22 = sqrt((Xp2(1)-100)^2+(Xp2(3)-100)^2) ;
    horizion_total2 = [horizion21 ; horizion22];
    
    H2=[Xp2(1)/horizion21 ,           0 ,     Xp2(3)/horizion21,                0  ;
            (Xp2(1)-100)/horizion22 , 0 ,      (Xp2(3)-100)/horizion22,       0  ];
    
    K2 = Pp2 * H2' / (H2 * Pp2 * H2' + R) ;
    deltaX2 = K2 * ( Z - horizion_total2);  %   Z為修正過後的pseudorange
    Xu2 = Xp2 + deltaX2 ;        
    Vx2 = Xu2 - Xp2;
    
    I = eye( size(Xu2,1) ) ;
    Pu2 = (I - K2 * H2) * Pp2;
    
    [ Xu2 , Pu2 ] = KF_refine_pseudo( Xu2 , Z ,Pu2 , R , Vx2);
        
    Xxyz2(i,:)=Xu2(:,1);
    
end
toc

figure( 1 ) ;
plot( effTimeArray , Xxyz1( effTimeArray , 1 ) , 'k-' ) ;
grid on ;
hold on ;
plot( effTimeArray , Xref_run( effTimeArray ,1 ) , 'g--' ) ;
title( 'KF' ) ;
print( '-dpng',  'KF' , '-r600'  ) ;       %Change "-r600" to the required DPI

%{
figure( 2 ) ;
plot( effTimeArray , Zxyz( effTimeArray ,1 ) , 'k-' ) ;
grid on ;
%}
figure( 3 ) ;
plot( effTimeArray , Xxyz1( effTimeArray , 1 )-Xref_run( effTimeArray ,1 ) , 'k-' ) ;
grid on ;

figure( 11 ) ;
plot( effTimeArray , Xxyz2( effTimeArray , 1 ) , 'k-' ) ;
grid on ;
hold on ;
plot( effTimeArray , Xref_run( effTimeArray ,1 ) , 'g--' ) ;
title( 'RKF' ) ;
print( '-dpng',  'RKF' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 13 ) ;
plot( effTimeArray , Xxyz2( effTimeArray , 1 )-Xref_run( effTimeArray ,1 ) , 'k-' ) ;
grid on ;
%{
figure( 4 ) ;
plot( effTimeArray , Xxyz( effTimeArray , 3 ) , 'k-' ) ;
grid on ;
%hold on ;

%plot( effTimeArray , Xref_run( effTimeArray ,2 ) , 'g--' ) ;
%}