clear ;
clc;
%%

runtime = 10000;
%%
X_accurate=[2 1 1]';
H=[cos(30/180*pi) sin(30/180*pi) 1 ;
    cos(45/180*pi) sin(45/180*pi) 1 ;
    cos(60/180*pi) sin(60/180*pi) 1 ;
    cos(75/180*pi) sin(75/180*pi) 1 ];
Z_accurate=H*X_accurate;
DOP = inv(H'*H);
PDOP=sqrt( DOP(1,1)+DOP(2,2) )
%%
for i = 1 : runtime
    randerror = 0.5*randn(1,4);
    Z = Z_accurate + randerror';
    X = (H'*H)\H'*Z;
    Xmatrix = X ;
end
dX = bsxfun( @minus , Xmatrix , X_accurate ) ;
std = (sum( dX .^2 , 2 )/runtime) .^0.5 