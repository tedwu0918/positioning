%---------------------顯示參數設定---------------------
[ Data_start_end_obs ] = { Data sampling_time_end obs  }
sum0 = 0;
%--------------------將xyz轉成enu------------------------
rec_pos_act = mean( Xxyz(effTimeArray,:),1 )       %計算STD
Xenu = zeros( runtime , 3 ) ;
for i = 1 : sampling_time_end
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
    sum0 = sum0 + Xenu( i , : ).^2 ;
end
%Xenu = Xenu - solid_tide_error ;

%--------------------計算STD & RMS------------------------
std = sqrt( sum0/runtime ) ;
%std = sqrt( var( Xenu(effTimeArray,:))  ) ;

%mean_errorE = mean(Xenu(effTimeArray,1))
%mean_errorN = mean(Xenu(effTimeArray,2))
%mean_errorU = mean(Xenu(effTimeArray,3))
%mean_errorH = sqrt(mean_errorE^2 + mean_errorN^2)
%rms = bsxfun( @power , mean( Xenu(effTimeArray , :).^2,1 ) , 0.5 ) ;

%-----------------四捨五入(round)取到小數點第X位---------------------
std_enu = round( std .* 1e4 ) ./ 1e4 ;
std_en = round( sqrt( sum( std(1:2).^2 ) ) * 1e4 ) / 1e4 ;
STD = [ std_enu' ; std_en' ]

%rms_enu = round( rms .* 1e4 ) ./ 1e4 ;
%rms_en = round( sqrt( sum( rms(1:2).^2 ) ) * 1e4 ) / 1e4 ;
%RMS = [ rms_enu' ; rms_en' ] ;

%error = [ std' , rms' ]
%Error = [ STD , RMS ]

figure( 1 ) ;
plot( effTimeArray , PDOP( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'PDOP Value' ) ;
title( 'PDOP (GPS)' )
print( '-dpng',  'PDOP_GPS' , '-r600'  ) ;

%--------------------Number of visible satellites----------------------
figure( 2 ) ;
plot( effTimeArray , nsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
ylim( [ min( nsat_Mtrix(effTimeArray) )-0.5  max( nsat_Mtrix(effTimeArray) )+0.5 ]  ) ;
title( 'Number of Visible Satellites (GPS)' )
print( '-dpng',  'Number of visible satellites_GPS' , '-r600'  ) ;
%}
%--------------------2D plot----------------------
figure( 3 ) ;
plot( Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East Error (m)' ) ;
ylabel( 'North Error (m)') ;
axis equal ;
%axis( [-3 , 3 ,-3 , 3 ] ) ;
title( 'Scatter Plot (GPS)' ) ;
print( '-dpng',  'Scatter plot_GPS' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 4 );
subplot(3,1,1);
plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS)');
%axis( [xlim,-2,2] ) ;
grid on;                        hold on;

subplot(3,1,2);
plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-3,3] ) ;
grid on;                        hold on;

subplot(3,1,3);
plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
%axis( [xlim,-5,5] ) ;
grid on;                        hold on;
print( '-dpng',  'ENU_GPS_LS' , '-r600' ) ;       %Change "-r600" to the required DPI

%%

%-------------------------------MSE----------------------------------------------------------
sum0 = 0;
%rec_pos_act = [-2424425.1013  5377188.1768  2418617.7454 ];
%rec_pos_act = [-2417142.8718  5382345.4730  2415036.9459  ];
%rec_pos_act = lla2ecef([25.149442  121.777873  73.6], 'WGS84' ) ;
rec_pos_act = lla2ecef([22.4259617417365 , 114.211355684335 ,53.2231195336208], 'WGS84' ) ;  %華大測站參考位置
Xenu = zeros( runtime , 3 ) ;
for i = 1 : sampling_time_end
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
    sum0 = sum0 + Xenu( i , : ).^2 ;
end
%Xenu = Xenu - solid_tide_error ;

%--------------------計算STD & RMS------------------------
RMSE = sqrt( sum0/runtime ) ;
%std = sqrt( var( Xenu(effTimeArray,:))  ) ;

%mean_errorE = mean(Xenu(effTimeArray,1))
%mean_errorN = mean(Xenu(effTimeArray,2))
%mean_errorU = mean(Xenu(effTimeArray,3))
%mean_errorH = sqrt(mean_errorE^2 + mean_errorN^2)
%rms = bsxfun( @power , mean( Xenu(effTimeArray , :).^2,1 ) , 0.5 ) ;

%-----------------四捨五入(round)取到小數點第X位---------------------
RMSE_enu = round( RMSE .* 1e4 ) ./ 1e4 ;
RMSE_en = round( sqrt( sum( RMSE(1:2).^2 ) ) * 1e4 ) / 1e4 ;
RMSE = [ RMSE_enu' ; RMSE_en' ]
%{
figure( 5 ) ;
plot( Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East Error (m)' ) ;
ylabel( 'North Error (m)') ;
axis equal ;
%axis( [-1.5 , 1.5 ,-1.5 , 1.5 ] ) ;
title( 'Scatter Plot (GPS)' ) ;
print( '-dpng',  'Scatter plot_GPS_MSE' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 6 );
subplot(3,1,1);
plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS)');
%axis( [xlim,-2,1] ) ;
grid on;                        hold on;

subplot(3,1,2);
plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-2,2] ) ;
grid on;                        hold on;

subplot(3,1,3);
plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
%axis( [xlim,-5,5] ) ;
grid on;                        hold on;
print( '-dpng',  'ENU_GPS_LS_MSE' , '-r600' ) ;       %Change "-r600" to the required DPI
%}
%%
%{
figure( 11 ) ;
plot( effTimeArray , run_delta_B(6, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
%ylim( [ -0.5 1.5 ]  ) ;
title( 'BD6' )
print( '-dpng',  'BD6' , '-r600'  ) ;
figure( 12 ) ;
plot( effTimeArray , run_delta_B(7, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
%ylim( [ -0.5 1.5 ]  ) ;
title( 'BD7' )
print( '-dpng',  'BD7' , '-r600'  ) ;
figure( 13 ) ;
plot( effTimeArray , run_delta_B(8, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
%ylim( [ -0.5 1.5 ]  ) ;
title( 'BD8' )
print( '-dpng',  'BD8' , '-r600'  ) ;
figure( 14 ) ;
plot( effTimeArray , run_delta_B(9, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
%ylim( [ -0.5 1.5 ]  ) ;
title( 'BD9' )
print( '-dpng',  'BD9' , '-r600'  ) ;
figure( 15 ) ;
plot( effTimeArray , run_delta_B(10, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
%ylim( [ -0.5 1.5 ]  ) ;
title( 'BD10' )
print( '-dpng',  'BD10' , '-r600'  ) ;
figure( 16 ) ;
plot( effTimeArray , run_delta_B(13, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
%ylim( [ -0.5 1.5 ]  ) ;
title( 'BD13' )
print( '-dpng',  'BD13' , '-r600'  ) ;

figure( 7 ) ;
plot( effTimeArray , cycle_slip_time(27, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'cycle-slip' ) ;
ylim( [ -0.5 1.5 ]  ) ;
title( 'GP27' )
print( '-dpng',  'GP27' , '-r600'  ) ;
figure( 8 ) ;
plot( effTimeArray , cycle_slip_time(31, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'cycle-slip' ) ;
ylim( [ -0.5 1.5 ]  ) ;
title( 'GP31' )
print( '-dpng',  'GP31' , '-r600'  ) ;
figure( 9 ) ;
plot( effTimeArray , cycle_slip_time(21, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'cycle-slip' ) ;
ylim( [ -0.5 1.5 ]  ) ;
title( 'GP21' )
print( '-dpng',  'GP21' , '-r600'  ) ;
figure( 10 ) ;
plot( effTimeArray , cycle_slip_time(14, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'cycle-slip' ) ;
ylim( [ -0.5 1.5 ]  ) ;
title( 'GP14' )
print( '-dpng',  'GP14' , '-r600'  ) ;
figure( 11 ) ;
plot( effTimeArray , cycle_slip_time(22, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'cycle-slip' ) ;
ylim( [ -0.5 1.5 ]  ) ;
title( 'GP22' )
print( '-dpng',  'GP22' , '-r600'  ) ;
figure( 12 ) ;
plot( effTimeArray , cycle_slip_time(29, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'cycle-slip' ) ;
ylim( [ -0.5 1.5 ]  ) ;
title( 'GP29' )
print( '-dpng',  'GP29' , '-r600'  ) ;
figure( 13 ) ;
plot( effTimeArray , cycle_slip_time(32, effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'cycle-slip' ) ;
ylim( [ -0.5 1.5 ]  ) ;
title( 'GP32' )
print( '-dpng',  'GP32' , '-r600'  ) ;
%}
%-------------sat CP-----------------------------------------
%{
figure( 11 ) ;
plot( effTimeArray , run_delta( 3 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP3' )
print( '-dpng',  'GP3-error' , '-r600'  ) ;

figure( 12 ) ;
plot( effTimeArray , run_delta( 14 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP14' )
print( '-dpng',  'GP14-error' , '-r600'  ) ;

figure( 13 ) ;
plot( effTimeArray , run_delta( 16 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP16' )
print( '-dpng',  'GP16-error' , '-r600'  ) ;

figure( 14 ) ;
plot( effTimeArray , run_delta( 21 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP21' )
print( '-dpng',  'GP21-error' , '-r600'  ) ;

figure( 15 ) ;
plot( effTimeArray , run_delta( 22 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP22' )
print( '-dpng',  'GP22-error' , '-r600'  ) ;

figure( 16 ) ;
plot( effTimeArray , run_delta( 26 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP26' )
print( '-dpng',  'GP26-error' , '-r600'  ) ;

figure( 17 ) ;
plot( effTimeArray , run_delta( 27 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP27' )
print( '-dpng',  'GP27-error' , '-r600'  ) ;

figure( 18 ) ;
plot( effTimeArray , run_delta( 31 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP31' )
print( '-dpng',  'GP31-error' , '-r600'  ) ;

figure( 19 ) ;
plot( effTimeArray , run_delta( 23 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP23' )
print( '-dpng',  'GP23-error' , '-r600'  ) ;

figure( 20 ) ;
plot( effTimeArray , run_delta( 29 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP29' )
print( '-dpng',  'GP29-error' , '-r600'  ) ;

figure( 21 ) ;
plot( effTimeArray , run_delta( 32 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'G32' )
print( '-dpng',  'GP32-error' , '-r600'  ) ;

figure( 22 ) ;
plot( effTimeArray , run_delta( 29 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'GP29' )
print( '-dpng',  'GP29-error' , '-r600'  ) ;

figure( 23 ) ;
plot( effTimeArray , run_delta_B( 32 ,effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'G32' )
print( '-dpng',  'GP32-error' , '-r600'  ) ;



figure( 10 ) ;
plot( effTimeArray , run_recbias( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'error(m)' ) ;
title( 'receiver bias' )
print( '-dpng',  'receiver bias' , '-r600'  ) ;
%}
%{
figure( 101 ) ;
plot( EL_error(2,:) , Rhoc_error(2,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP2');
print( '-dpng',  'GP2' , '-r600' ) ;
figure( 102 ) ;
plot( EL_error(5,:) , Rhoc_error(5,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP5');
print( '-dpng',  'GP5' , '-r600' ) ;
figure( 103 ) ;
plot( EL_error(6,:) , Rhoc_error(6,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP6');
print( '-dpng',  'GP6' , '-r600' ) ;
figure( 104 ) ;
plot( EL_error(9,:) , Rhoc_error(9,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP9');
print( '-dpng',  'GP9' , '-r600' ) ;
figure( 105 ) ;
plot( EL_error(12,:) , Rhoc_error(12,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP12');
print( '-dpng',  'GP12' , '-r600' ) ;
figure( 106 ) ;
plot( EL_error(13,:) , Rhoc_error(13,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP13');
print( '-dpng',  'GP13' , '-r600' ) ;
figure( 107 ) ;
plot( EL_error(15,:) , Rhoc_error(15,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP15');
print( '-dpng',  'GP15' , '-r600' ) ;
figure( 108 ) ;
plot( EL_error(17,:) , Rhoc_error(17,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP17');
print( '-dpng',  'GP17' , '-r600' ) ;
figure( 109 ) ;
plot( EL_error(19,:) , Rhoc_error(19,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP19');
print( '-dpng',  'GP19' , '-r600' ) ;
figure( 110 ) ;
plot( EL_error(20,:) , Rhoc_error(20,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP20');
print( '-dpng',  'GP20' , '-r600' ) ;
figure( 111 ) ;
plot( EL_error(25,:) , Rhoc_error(25,:) ) ;
grid on ;
xlabel( 'elevation' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP25');
print( '-dpng',  'GP25' , '-r600' ) ;

figure( 151 ) ;
scatter( CNR_error(2,:) , Rhoc_error(2,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP2');
print( '-dpng',  'CNR-GP2' , '-r600' ) ;

figure( 152 ) ;
scatter( CNR_error(5,:) , Rhoc_error(5,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP5');
print( '-dpng',  'CNR-GP5' , '-r600' ) ;

figure( 153 ) ;
scatter( CNR_error(6,:) , Rhoc_error(6,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP6');
print( '-dpng',  'CNR-GP6' , '-r600' ) ;

figure( 154 ) ;
scatter( CNR_error(9,:) , Rhoc_error(9,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP9');
print( '-dpng',  'CNR-GP9' , '-r600' ) ;

figure( 155 ) ;
scatter( CNR_error(12,:) , Rhoc_error(12,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP12');
print( '-dpng',  'CNR-GP12' , '-r600' ) ;

figure( 156 ) ;
scatter( CNR_error(13,:) , Rhoc_error(13,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP13');
print( '-dpng',  'CNR-GP13' , '-r600' ) ;

figure( 157 ) ;
scatter( CNR_error(15,:) , Rhoc_error(15,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP15');
print( '-dpng',  'CNR-GP15' , '-r600' ) ;

figure( 158 ) ;
scatter( CNR_error(17,:) , Rhoc_error(17,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP17');
print( '-dpng',  'CNR-GP17' , '-r600' ) ;

figure( 159 ) ;
scatter( CNR_error(19,:) , Rhoc_error(19,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP19');
print( '-dpng',  'CNR-GP19' , '-r600' ) ;

figure( 160 ) ;
scatter( CNR_error(20,:) , Rhoc_error(20,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP20');
print( '-dpng',  'CNR-GP20' , '-r600' ) ;

figure( 161 ) ;
scatter( CNR_error(25,:) , Rhoc_error(25,:) ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error' ) ;
grid on;                        hold on;
title('GP25');
print( '-dpng',  'CNR-GP25' , '-r600' ) ;

figure( 91 ) ;
plot( TimeArray , Rcovariance( 3,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP3');
print( '-dpng',  'GP3' , '-r600' ) ;

figure( 92 ) ;
plot( TimeArray , Rcovariance( 14,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP14');
print( '-dpng',  'GP14' , '-r600' ) ;

figure( 93 ) ;
plot( TimeArray , Rcovariance( 16,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP16');
print( '-dpng',  'GP16' , '-r600' ) ;


figure( 94 ) ;
plot( TimeArray , Rcovariance( 22,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP22');
print( '-dpng',  'GP22' , '-r600' ) ;


figure( 95 ) ;
plot( TimeArray , Rcovariance( 23,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP23');
print( '-dpng',  'GP23' , '-r600' ) ;


figure( 96 ) ;
plot( TimeArray , Rcovariance( 26,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP26');
print( '-dpng',  'GP26' , '-r600' ) ;


figure( 97 ) ;
plot( TimeArray , Rcovariance( 31,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP31');
print( '-dpng',  'GP31' , '-r600' ) ;

figure( 98 ) ;
plot( TimeArray , Rcovariance( 32,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP32');
print( '-dpng',  'GP32' , '-r600' ) ;

figure( 99 ) ;
plot( TimeArray , Rcovariance(29,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP29');
print( '-dpng',  'GP29' , '-r600' ) ;

figure( 100 ) ;
plot( TimeArray , Rcovariance( 27,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP27');
print( '-dpng',  'GP27' , '-r600' ) ;

figure( 101 ) ;
plot( TimeArray , Rcovariance( 21,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP21');
print( '-dpng',  'GP21' , '-r600' ) ;

figure( 102 ) ;
plot( TimeArray , Rcovariance( 8,TimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'covariance R' ) ;
grid on;                        hold on;
title('GP8');
print( '-dpng',  'GP8' , '-r600' ) ;
%}
%{
figure( 99 ) ;
plot( 1:length(allVk) , allVk( : ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
%axis( [xlim,-0.5,1.5] ) ;
%}