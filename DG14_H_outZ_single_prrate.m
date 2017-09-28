function [ H ,Z  ] = DG14_H_outZ_single_prrate( Xs, Vs, X , ARho0_G ,   Ds  ,nsat,rec_bias )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
dX = bsxfun(@minus , Xs , X ) ;                                   % Xs - X
Nor = sum(dX .^2, 2) .^0.5 ;                                        % || Xs - X ||
Unit_Mtrix = bsxfun(@rdivide , dX , Nor ) ;               % ( Xs - X ) / || Xs - X ||

doppler = Unit_Mtrix.*Vs ;
V_rate = doppler(:,1)+ doppler(:,2)+doppler(:,3);
Val_prr = zeros(nsat,1)+V_rate;

Val_pr = sum(dX .^2, 2) .^0.5+rec_bias ;  %   Val = || Xs - X || + b
Val = [Val_pr ; -Val_prr];
%%
%for Ds

Zprr = Ds ;
%%
% for pr
Zmeasurement = ARho0_G  ;
%%
%for output Z

Z = [Zmeasurement ; Zprr];
%%
%for H matrix
H = zeros( size(Xs,1)*2 , 8 ) ;
H( 1:nsat , [1,3,5] ) = -Unit_Mtrix ;               % ( X - Xs ) / || Xs - X ||     for pr
H(nsat+1:nsat*2 , [2,4,6] ) = Unit_Mtrix ;               % ( Xs - X ) / || Xs - X || for prr

H( 1:nsat , 7) = 1 ;                            %bias
H( nsat+1:nsat*2 , 8) = -1 ;            %drift
%H( nsat*2+1 , 7) = 1 ;                      %clock bias to clock bias


end

