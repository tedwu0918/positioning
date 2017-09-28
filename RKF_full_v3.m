clear
clc
close all
% addpath Functions
tic

%%

T=0.01 ;

Xu1 =[45; pi-2 ; 5 ;1.5];
Pu1 = blkdiag( 5,2,5,3) ;
 
runtime =10000;
effTimeArray = 1 : runtime ;
  
ff = [ 1 T ; 0 1 ] ;
fy = blkdiag( ff , ff ) ;
 
Z = zeros(2,1);
%%
Xu2 =[45; pi-2 ; 5 ;1.5];
Pu2 = blkdiag( 5,2,5,3 ) ;

%%
sigma1 = 0.3;
sigma2 = 0.3;
 
Qxyz1 = sigma1^2 * [T^3/3      T^2/2 ;
                                        T^2/2        T      ] ;
Qxyz2 = sigma2^2 * [T^3/3      T^2/2 ;
                                        T^2/2        T      ] ;
Q = blkdiag( Qxyz1 , Qxyz2 ) ;

%%
% for R (measurement noise)
Rhoerror = 0.01;

for i= 1: runtime
    Xnoise =0.2*randn(1,1);    Ynoise =0.2*randn(1,1);
    if Xnoise>5
        Xref =50*sin(2*pi*i*T/100)+50+5 ;
    elseif  Xnoise<-5      Xref =50*sin(2*pi*i*T/100)+50 -5 ;
    else
        Xref =50*sin(2*pi*i*T/100)+50 +Xnoise;
    end
    
     if Ynoise>5    Yref = i*T +5 ;
    elseif  Ynoise<-5      Yref =i*T - 5 ;
    else
         Yref =i*T + Ynoise;
    end
      
    horizion_ref1 = sqrt(Xref^2+Yref^2) ;    
    horizion_ref2 = sqrt( (Xref-100)^2+(Yref-100)^2) ;
    
    Z1noise = 0.01*randn(1,1);    Z2noise = 0.01*randn(1,1);
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
    Xu1 = Xp1 +  K1 * ( Z - horizion_total1) ;        %   Z為修正過後的pseudorange
    
    I = eye( size(Xu1,1) ) ;
    Pu1 = (I - K1 * H1) * Pp1;
    
    Xref_run(i,:) = [Xref Yref ];
    Xxyz1(i,:)=Xu1(:,1);
    Zxyz(i,:)=Z(:,1)';
    SE1(i)=sqrt( (Xu1(1)-Xref)^2+ (Xu1(3)-Yref)^2);
    
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
    
    I = eye( size(Xu2,1) ) ;
    Pu2 = (I - K2 * H2) * Pp2;
    
     %initial setting for refinement
    loop = 0;
    
    Xp_refine = Xu2 ;
    Pp_refine = Pu2 ;
    Xu_refine = Xu2 ;
    Pu_refine = Pu2 ;
    norm_delta_p =  norm(deltaX2);
    
    horizion1 = sqrt(Xp_refine(1)^2+Xp_refine(3)^2) ;
    horizion2 = sqrt((Xp_refine(1)-100)^2+(Xp_refine(3)-100)^2) ;
    horizion_total = [horizion1 ; horizion2];
    
    probe = [ 0, 0, 0, 0, 0, 0, 1, 0,0,0;
                      0, 0, 0, 0, 0, 0, 0, 1,0,0;
                      0, 0, 0, 0, 0, 0, 0, 0,1,0;
                      0, 0, 0, 0, 0, 0, 0, 0,0,1 ];
    bigmat = [ Pp_refine(1,1), Pp_refine(1,2), Pp_refine(1,3), Pp_refine(1,4),0,0,1,0,0,0;
                        Pp_refine(2,1), Pp_refine(2,2), Pp_refine(2,3), Pp_refine(2,4),0,0,0,1,0,0;
                        Pp_refine(3,1), Pp_refine(3,2), Pp_refine(3,3), Pp_refine(3,4),0,0,0,0,1,0;
                        Pp_refine(4,1), Pp_refine(4,2), Pp_refine(4,3), Pp_refine(4,4),0,0,0,0,0,1;
                        0,0,0,0, R(1,1), R(1,2), Xp_refine(1)/horizion1 ,   0 ,   Xp_refine(3)/horizion1,   0;
                        0,0,0,0, R(2,1), R(2,2), (Xp_refine(1)-100)/horizion2 , 0 ,   (Xp_refine(3)-100)/horizion2,  0;
                        1,0,0,0,  Xp_refine(1)/horizion1, (Xp_refine(1)-100)/horizion2 ,0,0,0,0;
                        0,1,0,0,  0,                                        0,                                                 0,0,0,0;
                        0,0,1,0,  Xp_refine(3)/horizion1,  (Xp_refine(3)-100)/horizion2,0,0,0,0;
                        0,0,0,1,  0,                                        0,                                                 0,0,0,0 ];
    virtualmeas = [0;0;0;0; Z(1) - horizion_total(1); Z(2) - horizion_total(2); 0;0;0;0 ];
    
  %  deter = det(bigmat);
    pseudoinv = inv(bigmat);
    delta_X_refine =  probe * pseudoinv * virtualmeas;
    norm_delta_n = norm(delta_X_refine([1,3])) ;
    threshold = sqrt(trace(Pp_refine));
    
    while norm_delta_n <= norm_delta_p * 0.99 && norm_delta_n >= threshold*0.05  %&&  loop <=3
        Xu_refine = Xp_refine + delta_X_refine ;        % update
        Pu_refine = - probe * pseudoinv * probe';
        
        norm_delta_p = norm_delta_n;
        Xp_refine = Xu_refine;
        Pp_refine = Pu_refine;
        
        horizion1 = sqrt(Xu_refine(1)^2+Xu_refine(3)^2) ;
        horizion2 = sqrt((Xu_refine(1)-100)^2+(Xu_refine(3)-100)^2) ;
        horizion_total = [horizion1 ; horizion2];
      horizion_total 
        probe = [ 0, 0, 0, 0, 0, 0, 1, 0,0,0;
                            0, 0, 0, 0, 0, 0, 0, 1,0,0;
                            0, 0, 0, 0, 0, 0, 0, 0,1,0;
                            0, 0, 0, 0, 0, 0, 0, 0,0,1 ];
        bigmat = [ Pu_refine(1,1), Pu_refine(1,2), Pu_refine(1,3), Pu_refine(1,4),0,0,1,0,0,0;
                            Pu_refine(2,1), Pu_refine(2,2), Pu_refine(2,3), Pu_refine(2,4),0,0,0,1,0,0;
                            Pu_refine(3,1), Pu_refine(3,2), Pu_refine(3,3), Pu_refine(3,4),0,0,0,0,1,0;
                            Pu_refine(4,1), Pu_refine(4,2), Pu_refine(4,3), Pu_refine(4,4),0,0,0,0,0,1;
                            0,0,0,0, R(1,1), R(1,2), Xu_refine(1)/horizion1 ,   0 ,   Xu_refine(3)/horizion1,   0;
                            0,0,0,0, R(2,1), R(2,2), (Xu_refine(1)-100)/horizion2 , 0 ,   (Xu_refine(3)-100)/horizion2,  0;
                            1,0,0,0,  Xu_refine(1)/horizion1, (Xu_refine(1)-100)/horizion2 ,0,0,0,0;
                            0,1,0,0,  0,                                        0,                                                 0,0,0,0;
                            0,0,1,0,  Xu_refine(3)/horizion1,  (Xu_refine(3)-100)/horizion2,0,0,0,0;
                            0,0,0,1,  0,                                        0,                                                 0,0,0,0 ];
        virtualmeas = [0;0;0;0; Z(1) - horizion_total(1); Z(2) - horizion_total(2); 0;0;0;0 ];
        pseudoinv = inv (bigmat);
        delta_X_refine =  probe * pseudoinv * virtualmeas;
               
        norm_delta_n =  norm(delta_X_refine([1,3])) ;    %記錄前一個step的norm delta
        threshold = sqrt(trace(Pp_refine));
        loop=loop+1;
        
    end
   % loop
    
    Xu2 = Xu_refine;
    Pu2 = Pu_refine;
    
    Xxyz2(i,:)=Xu2(:,1);
    SE2(i)=sqrt( (Xu2(1)-Xref)^2+ (Xu2(3)-Yref)^2);
    
end
%% 
MSSE1 = sum(SE1)/runtime
MSSE2 = sum(SE2)/runtime
MSSEReduction = MSSE1 - MSSE2
ImprovedPercentage = MSSEReduction / MSSE1 * 100
toc

figure( 1 ) ;
plot( effTimeArray , Xxyz1( effTimeArray , 1 ) , 'k-' ) ;
grid on ;
hold on ;
plot( effTimeArray , Xref_run( effTimeArray ,1 ) , 'g--' ) ;
title( 'KF' ) ;
print( '-dpng',  'KF' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 2 ) ;
plot( effTimeArray , Xxyz1( effTimeArray , 3 ) , 'g-' ) ;
grid on ;
hold on ;
plot( effTimeArray , Xref_run( effTimeArray ,2 ) , 'r--' ) ;
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

figure( 12 ) ;
plot( effTimeArray , Xxyz2( effTimeArray , 3 ) , 'g-' ) ;
grid on ;
hold on ;
plot( effTimeArray , Xref_run( effTimeArray ,2 ) , 'r--' ) ;
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