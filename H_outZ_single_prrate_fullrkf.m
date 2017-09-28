function [ H ,Z ,Val ] = H_outZ_single_prrate_fullrkf( Xs,  Xr ,  ARho0_G, nsat , Vs, Ds)
%UNTITLED5 Summary of this function goes here
X=Xr([1 2 3 ])';
b=Xr(4);
V=Xr([5 6 7])';
%   Detailed explanation goes here
dX = bsxfun(@minus , Xs , X ) ;                                   % Xs - X
Nor = sum(dX .^2, 2) .^0.5 ;                                        % || Xs - X ||
Unit_Mtrix = bsxfun(@rdivide , dX , Nor ) ;               % ( Xs - X ) / || Xs - X ||

Val_pr =Nor+b ;  %   Val = || Xs - X || + b

doppler = Unit_Mtrix.*Vs ;
V_rate = doppler(:,1)+ doppler(:,2)+doppler(:,3);
Val_prr = zeros(nsat,1)+V_rate;
Val = [Val_pr ; -Val_prr];
%% 
%for Ds
Zprr=Ds;
% for ps
Zmeasurement = ARho0_G ;
%for output Z
Z = [Zmeasurement ; Zprr];
%%
%for H matrix
g=[(-1./Nor+(Xs(:,1)-X(1)).^2)/(Nor.^3)*(V(1)-Vs(:,1))  +  (Xs(:,1)-X(1)).*(Xs(:,2)-X(2))/(Nor.^3)*(V(2)-Vs(:,2))  +  (Xs(:,1)-X(1)).*(Xs(:,3)-X(3))/(Nor.^3)*(V(3)-Vs(:,3)),...
       (Xs(:,1)-X(1)).*(Xs(:,2)-X(2))/(Nor.^3)*(V(1)-Vs(:,1))   -  (1./Nor+(Xs(:,2)-X(2)).^2)/(Nor.^3)*(V(2)-Vs(:,2))  +  (Xs(:,2)-X(2)).*(Xs(:,3)-X(3))/(Nor.^3)*(V(3)-Vs(:,3)),...
    (Xs(:,3)-X(3)).*(Xs(:,1)-X(1))/(Nor.^3)*(V(1)-Vs(:,1))  +  (Xs(:,2)-X(2)).*(Xs(:,3)-X(3))/(Nor.^3)*(V(2)-Vs(:,2))  -  (1./Nor+(Xs(:,3)-X(3)).^2)/(Nor.^3)*(V(3)-Vs(:,3))];

H =  zeros( size(Xs,1)*2 , 8 ) ;
H( 1:nsat , [1,2,3] ) = -Unit_Mtrix ;               % ( X - Xs ) / || Xs - X ||     for pr
H(1:nsat , [5,6,7] ) = 0 ;                 
H( 1:nsat , 4) = 1 ;                                          %bias
H(  1:nsat, 8) = 0 ;                                          %drift

H(nsat+1:nsat*2 , [1,2,3] ) = -g ;               % ( Xs - X ) / || Xs - X || for prr
H(nsat+1:nsat*2 , [5,6,7] ) =-Unit_Mtrix ;                 
H( nsat+1:nsat*2 , 4) =0 ;                            %bias
H(  nsat+1:nsat*2, 8) = 1 ;                           %drift



end

