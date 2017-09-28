format long g
clc ;
clear ;
close all ;

GMAX = 32 ;
BMAX = 14 ;

addpath Data_0322
obs_filename = '03283.csv' ;     % BDS�PGPS���ƾ��ɮ�

smaller_than_elevation = 15 ;
sampling_time_start = 10 ;                     %��ƨ��ˮɶ�(�}�l�I)
%sampling_time_end = 3935 ;                  %��ƨ��ˮɶ�(�����I)  

[ week_num , hwtime , useless , sat_sys , prn , cnr0 , TOW0 , Elapse_Epoch , Elapse_Code , pr , pr_rate , ADR , epr , sat_EL , sat_AZ , time ] = ...
    textread( obs_filename , '%f %f %c %c %f %f %f %f %f %f %f %f %f %f %f %f' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:���L�Ĥ@��

runtime = round( max( time ) ) - round ( min( time ) ) ;
%runtime = 2000 ;
count_times = ( min( time ) ) ;               % ��ƪ��ɨ�(���O�ɶ�,�u��ܦ����������ĴX�����)
count = 1 ;                                                % �ĴX����(�Ψ�Ū��raw data)
id_set = zeros( runtime , GMAX  ) ;
id_set_B = zeros( runtime , BMAX  ) ;

for i = 1 : runtime
    i
    sat_id = [] ; sat_id_B = [] ;
    
    judge = 1 ;
    j = 1 ;                            %�C�@����U��GPS����ƨ̧ǧ�i��
    k = 1 ;                           %�C�@����U��BDS����ƨ̧ǧ�i��
    
    while judge == 1
        if sat_EL( count ) < smaller_than_elevation
            count = count + 1 ;
        elseif time( count ) == count_times
            if char( sat_sys( count ) ) == 80 &&  prn( count ) ~= 193
                sat_id( j ) = prn( count ) ;
                id_set( i , prn( count ) ) = 1 ;
                                
                 j = j+1 ;
                 
            elseif char( sat_sys( count ) ) == 66  %&&  prn( count ) ~= 5
                sat_id_B(k)=prn(count);
                id_set_B( i , prn( count ) ) = 1 ;
                
                k = k+1;
              
            end
            
            count=count+1;       
        
        else        
            count_times=time(count);          
            judge=0;                
        end
            
         if count>length(time)
            count=count-1;
                break
        end    
              
    end
    
    nsat = length( sat_id ) ;
    sat_num( i ) = nsat ;
    nsat_B = length( sat_id_B ) ;
    sat_num_B( i ) = nsat_B ;
    nsat_total = nsat + nsat_B ;
    sat_num_total( i ) = nsat_total ;
    
end

sampling_time_end = i 

%------------------��ƨ��ˮɶ�----------------------
effTimeArray = sampling_time_start : sampling_time_end ;

%--------------------number of satellites----------------------
figure( 1 ) ;
plot( effTimeArray , sat_num_total( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'time(s)' ) ;
ylabel( 'number' ) ;
legend( ' Number of Satellites ' ) ;
title( 'GPS/BDS satellite' )
%print( '-dpng',  'Number of Daul' , '-r600'  ) ;

figure( 2 ) ;
plot( effTimeArray , sat_num( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'time(s)' ) ;
ylabel( 'number' ) ;
legend( 'Number of satellites' ) ;
title( 'GPS satellte' ) ;
%print( '-dpng',  'Number of GPS' , '-r600'  ) ;

figure( 3 ) ;
plot( effTimeArray , sat_num_B( effTimeArray ) , 'b-' ) ;
grid on ;
xlabel( 'time(s)' ) ;
ylabel( 'number' ) ;
legend( 'Number of satellites' ) ;
title( 'BDS satellite' ) ;
%print( '-dpng',  'Number of BDS' , '-r600'  ) ;
%{
%--------------------number of  single satellite----------------------
id = 4
figure( 4 ) ;
plot( effTimeArray , id_set( effTimeArray , id ) , 'b-' ) ;
grid on ;
xlabel( 'time(s)' ) ;
ylabel( 'Appear/Disappear' ) ;
title( 'GP' ) ; 

id_B = 10 
figure( 5 ) ;
plot( effTimeArray , id_set_B( effTimeArray , id_B ) , 'b-' ) ;
grid on ;
xlabel( 'time(s)' ) ;
ylabel( 'Appear/Disappear' ) ;
title( 'GB' ) ; 
%}