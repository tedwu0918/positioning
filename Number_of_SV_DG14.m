format long g
clc ;
clear ;
close all ;

addpath DG14_0622
obs_filename = '0622.csv' ;         % GPS���ƾ��ɮ�

%-----------Ū����������T------------------------
[ time , prn , SV_x , SV_y , SV_z , ADR , pr , pr_rate , SV_x_dot , SV_y_dot , SV_z_dot ] = ...
    textread( obs_filename , '%f %f %f %f %f %f %f %f %f %f %f ' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:���L�Ĥ@��

%-----------initialize--------------------------
runtime = round( max( time ) ) - round ( min( time ) ) ;
runtime = 9000 ;
count_times = ( min( time ) ) ;               % ��ƪ��ɨ�(���O�ɶ�,�u��ܦ����������ĴX�����)
count = 1 ;                                                % �ĴX����(Ū��row data)

for i = 1 : runtime

    id = [] ;
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %�C�@����U��GPS����ƨ̧ǧ�i��
    
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
    
    %----------------�C��(�C�@��)�w��ìP������---------------
    nsat = length( id ) ;
    nsat_Mtrix( i ) = nsat ;    
    
end

%---------------------��ƨ��ˮɶ�----------------------
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

