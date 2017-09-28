clear
clc
close all
addpath Functions
tic

%%


Pu = blkdiag( 10,1,10,1,10,1) ;
runtime =1000;
effTimeArray = 1 : runtime ;
Xu =[200;0;2;0;-1;0];
dXu =[1;0;0;0;0;0];
ff = [ 1 1 ; 0 1 ] ;
fy = blkdiag( ff , ff , ff ) ;
T=1 ;
Z = zeros(1,1);
c=3*10^8;
%%

sigma = 0.001 ;
%Sb = (1.1e-19)*c^2 ;
%Sd = ((pi^2)*56e-21)*c^2 ;%Single-frequency, single-receiver terristrial and spaceborne point positioning(p.67)
Sd = 0.001 ;%Single-frequency, single-receiver terristrial and spaceborne point positioning(p.67)

Qc = [Sd*T*T*T/3     Sd*T*T/2 ;
    Sd*T*T/2                            Sd*T ] ;
Qxyz = sigma^2 * [T^3/3      T^2/2 ;
    T^2/2        T      ] ;
Q = blkdiag( Qxyz , Qxyz , Qxyz) ;

%%
% for R
Rhoerror = 1 ;

%%
Xs = [] ;
ref_pos = [] ;
Ansat = [] ;
for i= 1: runtime
    Xref=100+50*sin(i/200*pi) +5*randn(1,1);
    
    
    
    Pp =  fy * Pu * fy.' + Q ;
    Xp = fy * Xu ;                  % fy為State transition matrix
    
    Z(1,1) = 4*rand(1,1)+ Xref;%+50*cos(i/100*pi);
    H=[1,0,0,0,0,0];
    
    
    R =  eye( size( Z , 1 ) ) * Rhoerror;
    
    %K=(inv(Pp)+H' * inv(R )* H)\H'*inv(R);
    K = Pp * H' / (H * Pp * H.' + R) ;
    dXu =  K * ( Z - H*Xp) ;        %   Z為修正過後的pseudorange
    
    I = eye( size(dXu,1),size(dXu,1) ) ;
    Pu = (I - K * H) * Pp;
    Xu =Xp + dXu;
    
    [ Xu , Pu ] = KF_refine_test2( Xu , Z ,Pu , R  ,dXu,H);
    
    
    Xref_run(i,1) = Xref;
    Xxyz(i,:)=Xu([1,3,5],1);
    Zxyz(i,:)=Z(1,1);
    
    %}
end
toc
figure( 1 ) ;
plot( effTimeArray , Xxyz( effTimeArray , 1 ) , 'k-' ) ;
grid on ;
hold on ;

plot( effTimeArray , Xref_run( effTimeArray ,1 ) , 'g--' ) ;

figure( 2 ) ;
plot( effTimeArray , Zxyz( effTimeArray ,1 ) , 'k-' ) ;
grid on ;

figure( 3 ) ;
plot( effTimeArray , Xxyz( effTimeArray , 1 )-Xref_run( effTimeArray ,1 ) , 'k-' ) ;
grid on ;


%{
figure( 2 ) ;
plot( effTimeArray , Xxyz( effTimeArray , 2 ) , 'k-' ) ;
grid on ;
hold on ;

plot( effTimeArray , Zxyz( effTimeArray , 2 ) , 'g--' ) ;

figure( 3 ) ;
plot( effTimeArray , Xxyz( effTimeArray , 3 ) , 'k-' ) ;
grid on ;
hold on ;

plot( effTimeArray , Zxyz( effTimeArray , 3 ) , 'g--' ) ;
%}

