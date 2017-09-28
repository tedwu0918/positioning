%---------------------顯示參數設定----------------------
[ Data_time_obs_mask_converge ] = { Data effTimeArray(1) effTimeArray(end) obs smaller_than_elevation converge_time Remark }

%--------------------將xyz轉成enu------------------------
rec_pos_act = mean( Xxyz(effTimeArray,:),1 ) ;      %計算STD
Xenu = zeros( runtime , 3 ) ;

for i = 1 : runtime 
    Xenu( i , : ) = xyz2enu( Xxyz( i , : ) , rec_pos_act )' ;
end
Xenu = Xenu - solid_tide_error ;

%--------------------計算STD & RMS------------------------
std = sqrt( var( Xenu(effTimeArray,:))  ) ;
%rms = bsxfun( @power , mean( Xenu(effTimeArray , :).^2,1 ) , 0.5 ) ;

%-----------------四捨五入(round)取到小數點第X位---------------------
std_enu = round( std .* 1e4 ) ./ 1e4 ;
std_en = round( sqrt( sum( std(1:2).^2 ) ) * 1e4 ) / 1e4 ; 
STD = [ std_enu' ; std_en' ] 

%rms_enu = round( rms .* 1e4 ) ./ 1e4 ;
%rms_en = round( sqrt( sum( rms(1:2).^2 ) ) * 1e4 ) / 1e4 ; 
%RMS = [ rms_enu' ; rms_en' ] ;

%error = [ std' , rms' ]
%Range = [ min( Xenu(effTimeArray,:) )' , max( Xenu(effTimeArray,:) )' ] 
%Error = [ STD , RMS ] 

%--------------------Number of visible satellites----------------------
figure( 1 ) ;
plot( effTimeArray , nsat_Mtrix( effTimeArray ) ,  'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
ylim( [ min( nsat_Mtrix(effTimeArray) )-0.5  max( nsat_Mtrix(effTimeArray) )+0.5 ]  ) ;
title( 'Number of Visible Satellites (GPS/BDS)' )
print( '-dpng',  'Number of Dual' , '-r600'  ) ;

%{
figure( 7 ) ;
plot( effTimeArray , nGsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
title( 'Number of Visible Satellites (GPS)' ) ;
print( '-dpng',  'Number of GPS' , '-r600'  ) ;

figure( 8 ) ;
plot( effTimeArray , nBsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'time(s)' ) ;
ylabel( 'number' ) ;
title( 'Number of Visible Satellites (BDS)' ) ;
print( '-dpng',  'Number of BDS' , '-r600'  ) ;
%}

figure( 2 ) ;
plot( effTimeArray , PDOP( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'PDOP Value' ) ;
title( 'PDOP (GPS/BDS)' )
print( '-dpng',  'PDOP_Dual' , '-r600'  ) ;

%--------------------2D plot----------------------
figure( 3 ) ;
plot( Xenu( effTimeArray , 1 ) , Xenu( effTimeArray , 2 ) , 'b.' ) ;
grid on ;
xlabel( 'East (m)' ) ;
ylabel( 'North (m)') ;
axis equal ;
%axis( [-6 , 6 ,-5 , 5 ] ) ;
title( 'Scatter Plot (GPS/BDS)' ) ;
print( '-dpng',  'Scatter plot_Dual_LS' , '-r600'  ) ;       %Change "-r600" to the required DPI

figure( 4 );
subplot(3,1,1);
plot( effTimeArray , Xenu( effTimeArray , 1 ) , 'r-' ) ;
xlabel('Time (s)');
ylabel('East (m)');
title('ENU位置與時間關係 (GPS/BDS)');
%axis( [xlim,-2,2] ) ;
grid on;                        hold on;

subplot(3,1,2);
plot(effTimeArray, Xenu( effTimeArray , 2 ) ,'b-');
xlabel('Time (s)');
ylabel('North (m)');
%axis( [xlim,-2,2] ) ;
grid on;                        hold on;

subplot(3,1,3);
plot(effTimeArray, Xenu( effTimeArray , 3 ) , 'g-' ) ;
xlabel('Time (s)');
ylabel('Up (m)');
%axis( [xlim,-5,5] ) ;
grid on;                        hold on;
print( '-dpng',  'ENU_Dual_LS' , '-r600' ) ;       %Change "-r600" to the required DPI

