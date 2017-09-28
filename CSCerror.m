clear ;
clc ;
close all ;
%%
CNRmatrix =[ ] ;
Rmatrix = [] ;
ELmatrix = [] ;
ERmatrix = [] ;

%%
%{
%0512
[CNR , R , EL ] =textread( 'origin0512error.csv' , '%d %f %f \n' , 'delimiter' , ','  ) ;
CNRmatrix = CNR ;
Rmatrix =R ;
ELmatrix = EL ;


clear CNR CNRerror EL ELerror ;

%0623
[CNR , R , EL ] =textread( 'origin0623error.csv' , '%d %f %f\n' , 'delimiter' , ','  ) ;
CNRmatrix =[ CNRmatrix ; CNR] ;
Rmatrix = [Rmatrix ; R] ;
ELmatrix = [ELmatrix ; EL] ;
clear CNR CNRerror EL ELerror ;


%%

CNRmatrix = round(CNRmatrix) ;
[CNRmatrix , ic] = sort(CNRmatrix) ;
CRhoc_error = Rmatrix(ic) ;

ELmatrix = round(ELmatrix) ;
[ELmatrix , id] = sort(ELmatrix) ;
ERhoc_error = Rmatrix(id) ;
%}

%%

%{
[CNR , R , EL , ER ] =textread( '071011hr_error.csv' , '%d %f %d %f\n' , 'delimiter' , ','  ) ;
CNRmatrix =[CNRmatrix ;CNR] ;
Rmatrix = [Rmatrix; R] ;
ELmatrix = [ELmatrix ;EL] ;
ERmatrix = [ERmatrix; ER] ;

clear  CNR  R  EL  ER 
%}
%

[CNR , R , EL , ER ] =textread( 'hksterror.csv' , '%d %f %f %f\n' , 'delimiter' , ','  ) ;
CNRmatrix =[CNRmatrix; CNR] ;
Rmatrix = [Rmatrix ;R] ;
ELmatrix = [ELmatrix ;EL] ;
ERmatrix = [ERmatrix ;ER] ;

clear  CNR  R  EL  ER 

%

[CNR , R , EL , ER ] =textread( 'hksserror.csv' , '%d %f %f %f\n' , 'delimiter' , ','  ) ;
CNRmatrix =[CNRmatrix ;CNR] ;
Rmatrix = [Rmatrix ;R] ;
ELmatrix = [ELmatrix; EL] ;
ERmatrix = [ERmatrix ;ER] ;

clear  CNR  R  EL  ER 

%

[CNR , R , EL , ER ] =textread( '0623error.csv' , '%d %f %f %f\n' , 'delimiter' , ','  ) ;
CNRmatrix =[CNRmatrix ;CNR] ;
Rmatrix = [Rmatrix ;R] ;
ELmatrix = [ELmatrix ;EL] ;
ERmatrix = [ERmatrix ;ER] ;

clear  CNR  R  EL  ER 


[CNR , R , EL , ER ] =textread( '0512error.csv' , '%d %f %d %f\n' , 'delimiter' , ','  ) ;
CNRmatrix =[CNRmatrix ;CNR] ;
Rmatrix = [Rmatrix ;R] ;
ELmatrix = [ELmatrix ;EL] ;
ERmatrix = [ERmatrix ;ER] ;

clear  CNR  R  EL  ER 

[CNR , R , EL , ER ] =textread( '0309error.csv' , '%d %f %d %f\n' , 'delimiter' , ','  ) ;
CNRmatrix =[CNRmatrix ;CNR] ;
Rmatrix = [Rmatrix ;R] ;
ELmatrix = [ELmatrix ;EL] ;
ERmatrix = [ERmatrix ;ER] ;

clear  CNR  R  EL  ER 
%}
%{
[CNR , R , EL ,  ] =textread( '0706errorBDS.csv' , '%d %f %f \n' , 'delimiter' , ','  ) ;
CNRmatrix =[CNRmatrix ;CNR] ;
Rmatrix = [Rmatrix ;R] ;
ELmatrix = [ELmatrix ;EL] ;


clear  CNR  R  EL  
%}
CNRmatrix = round(CNRmatrix) ;
[CNRmatrix , ic] = sort(CNRmatrix) ;
Rmatrix = Rmatrix(ic) ;

ELmatrix = round(ELmatrix) ;
[ELmatrix , id] = sort(ELmatrix) ;
ERmatrix = ERmatrix(id) ;

%%
%PLOT all
%{



figure( 11 ) ;
scatter(CNRmatrix , CNRerrormatrix );
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'covariance error(m)' ) ;
grid on;                        hold on;
title('CNR-CSC error');
print( '-dpng',  'CNR-CSC error' , '-r600' ) ;

figure( 12 ) ;
scatter(ELmatrix , ELerrormatrix );
grid on ;
xlabel( 'EL(Deg)' ) ;
ylabel( 'covariance error(m)' ) ;
grid on;                        hold on;
title('EL-CSC error');
print( '-dpng',  'EL-CSC error' , '-r600' ) ;
%}
%%



cref = CNRmatrix(1) ;
Eref = ELmatrix(1) ;

ccalculate = [] ;
Ecalculate = [] ;
n=1;x_matrix = [] ;mean_matrix = [] ;std_matrix = [];
ss = 1 ;Ex_matrix = [] ;Emean_matrix = [] ;Estd_matrix = [];
for i = 1 : length(CNRmatrix)
    if CNRmatrix(i)~= cref
        temp_mean = mean(ccalculate) ;
        temp_std = std(ccalculate);
        x_matrix = [x_matrix ; cref];
        mean_matrix = [mean_matrix ; temp_mean];
        std_matrix = [std_matrix ; temp_std];
        cref =  CNRmatrix(i);
        ccalculate = [];
        clear temp_mean temp_std ;
        n=1;
    end
    if ELmatrix(i)~= Eref
        temp_mean = mean(Ecalculate) ;
        temp_std = std(Ecalculate);
        Ex_matrix = [Ex_matrix ; Eref];
        Emean_matrix = [Emean_matrix ; temp_mean];
        Estd_matrix = [Estd_matrix ; temp_std];
        Eref =  ELmatrix(i);
        Ecalculate = [];
        clear temp_mean temp_std ;
        ss=1;
    end
    ccalculate(n ,1) = Rmatrix(i) ;
    Ecalculate(ss,1) = ERmatrix(i) ;
    n = n+1;
    ss = ss + 1 ;
    
    
    if i==length(CNRmatrix)
        temp_mean = mean(ccalculate) ;
        temp_std = std(ccalculate);
        x_matrix = [x_matrix ; cref];
        mean_matrix = [mean_matrix ; temp_mean];
        std_matrix = [std_matrix ; temp_std];
        ccalculate = [];
        clear temp_mean temp_std ;
        n=1;
        
        temp_mean = mean(Ecalculate) ;
        temp_std = std(Ecalculate);
        Ex_matrix = [Ex_matrix ; Eref];
        Emean_matrix = [Emean_matrix ; temp_mean];
        Estd_matrix = [Estd_matrix ; temp_std];
        Ecalculate = [];
        clear temp_mean temp_std ;
        ss=1;
    end
end

%{
figure( 11 ) ;
errorbar(x_matrix , mean_matrix , std_matrix );
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'covariance error(m)' ) ;
grid on;                        hold on;
title('CNR-CSC error');
print( '-dpng',  'CNR-CSC error' , '-r600' ) ;

figure( 12 ) ;
errorbar(Ex_matrix , Emean_matrix , Estd_matrix );
grid on ;
xlabel( 'EL(Deg)' ) ;
ylabel( 'covariance error(m)' ) ;
grid on;                        hold on;
title('EL-CSC error');
print( '-dpng',  'EL-CSC error' , '-r600' ) ;
%}

%%
%STD plot
 figure( 101 ) ;
plot(x_matrix,std_matrix);
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'variance(m)' ) ;
grid on;                        hold on;
title('CNR noise variance');
print( '-dpng',  'CNR noise variance' , '-r600' ) ;


figure( 102 ) ;
plot(Ex_matrix,Estd_matrix);
grid on ;
xlabel( 'EL(Deg)' ) ;
ylabel( 'variance(m)' ) ;
grid on;                        hold on;
title('EL noise variance');
print( '-dpng',  'EL noise variance' , '-r600' ) ;


%mean Plot

 figure( 103 ) ;
bar(x_matrix,abs(mean_matrix));
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error(m)' ) ;
grid on;                        hold on;
title('CNR mean noise');
print( '-dpng',  'CNR mean noise' , '-r600' ) ;


figure( 104 ) ;
bar(Ex_matrix,abs(Emean_matrix));
grid on ;
xlabel( 'EL(Deg)' ) ;
ylabel( 'error(m)' ) ;
grid on;                        hold on;
title('EL mean noise');
print( '-dpng',  'EL mean noise' , '-r600' ) ;