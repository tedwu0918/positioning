format long g
clc ;
clear ;
close all ;

addpath DG14_0622
obs_filename = '0622.csv' ;         % GPS的數據檔案

%-----------讀取接收機資訊------------------------
[ time , prn , SV_x , SV_y , SV_z , ADR , pr , pr_rate , SV_x_dot , SV_y_dot , SV_z_dot ] = ...
    textread( obs_filename , '%f %f %f %f %f %f %f %f %f %f %f ' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:跳過第一行

%-----------initialize--------------------------
runtime = round( max( time ) ) - round ( min( time ) ) ;
runtime = 9000 ;
count_times = ( min( time ) ) ;               % 資料的時刻(不是時間,只表示此為接收機第幾筆資料)
count = 1 ;                                                % 第幾行資料(讀取row data)

for i = 1 : runtime

    id = [] ;
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %每一秒都把各顆GPS的資料依序抓進來
    
    while judge == 1
        
        if time( count ) == count_times
            if prn( count ) ~= 193
                
                id( j , 1) = prn( count ) ;
                j = j+1 ;
                
            end
            count = count + 1 ;       
        
        else        
            count_times = time( count ) ;          
            judge = 0 ;                
        end
            
         if count>length( time )
            count = count - 1 ;
                break
        end    
              
    end
    
    %----------------每次(每一秒)定位衛星的顆數---------------
    nsat = length( id ) ;
    nsat_Mtrix( i ) = nsat ;    
    
end

%---------------------資料取樣時間----------------------
sampling_time_end = runtime ;
effTimeArray = 10 : sampling_time_end ;

%--------------------Number of visible satellites----------------------
figure( 2 ) ;
plot( effTimeArray , nsat_Mtrix( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'Time (s)' ) ;
ylabel( 'Number' ) ;
ylim( [ min( nsat_Mtrix(effTimeArray) )-0.5  max( nsat_Mtrix(effTimeArray) )+0.5 ]  ) ;
title( 'Number of Visible Satellites (GPS)' )
print( '-dpng',  'Number of visible satellites_GPS' , '-r600'  ) ;

