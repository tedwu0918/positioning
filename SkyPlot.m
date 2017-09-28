tic
format long g
clc ;
clear ;
close all ;

addpath Functions

GMAX = 32 ;     BMAX = 15 ;    

Data = 0322 ;
smaller_than_elevation = 10 ;                   %仰角小於elevation即濾掉

switch Data        
    case { 1110 }
        addpath Data_1110 ; 
        obs_filename = 'Obs_2015-11-10_10.01.24.csv' ; 
    case{ 11151 }
        addpath Data_1115_1 ;
        obs_filename = 'NTOU_1115_1.csv' ;
    case{ 11152 }
        addpath Data_1115_2 ;
        obs_filename = 'NTOU_1115_2.csv' ;
    case{ 11153 }
        addpath Data_1115_3 ;
        obs_filename = 'NTOU_1115_3.csv' ;
    case{ 1116 }
        addpath Data_1116 ;
        obs_filename = 'NTU_1116.csv' ;
    case{ 1210322 }
        addpath Data_ntou0322 ;
        obs_filename = 'ntou0322_10000.csv' ;
    case{ 0512 }
        addpath Data_0512 ;
        obs_filename = 'ntou0512.csv' ;
    case{ 0322 }
        addpath Data_ntou0322 ;
        obs_filename = 'ntou0322_16000.csv' ;
        case { 0913 }
        addpath Data_0913 ;
        obs_filename = 'HD9100_2HR_RCV1.16o' ;            
        
end

%-----------讀取接收機資訊------------------------
[ week_num , hwtime , useless , sat_sys , prn , cnr0 , TOW0 , Elapse_Epoch , Elapse_Code , pr , pr_rate , ADR , epr , sat_EL , sat_AZ , time ] = ...
    textread( obs_filename , '%f %f %c %c %f %f %f %f %f %f %f %f %f %f %f %f' , 'delimiter' , ',' , 'headerlines' , 1 ) ;%headerlines:跳過第一行

%-----------initialize--------------------------
runtime = round( max( time ) ) - round ( min( time ) ) ;
%runtime = 6093 ;
sat_num = zeros( 1 , runtime ) ;      sat_num_B = sat_num ;       sat_num_total = sat_num ;     rec_bias_total = sat_num ;
count_times = ( min( time ) ) ;               % 資料的時刻
count = 1 ;                                                % 第幾行資料(讀取row data)
id_set_G = zeros( runtime , GMAX );     id_set_B = zeros( runtime , BMAX );
appeared_id = zeros( 1 , GMAX ) ;
appeared_id_B = zeros( 1 , BMAX ) ;

beautify = 1 ;

%-----0322(within 9500s),0322(within 7200s)----------
GP_sat = [ 12,14,18,25,32 ] ;
BD_sat = [ 1,2,3,4,6,7,9,10 ] ;

for i = 1 : runtime
    
    if rem( i , 500 ) == 0
        i
    end
    
    sat_id = [] ;
    CNR = [] ;
    sat_elevation = [] ;
    sat_azimuth = [] ;
    
    sat_id_B = [] ;
    CNR_B = [] ;
    sat_elevation_B = [] ;
    sat_azimuth_B = [] ;
    
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %每一秒都把各顆GPS的資料依序抓進來
    k = 1 ;                           %每一秒都把各顆BDS的資料依序抓進來
    
    while judge == 1
        
        while sat_EL( count ) < smaller_than_elevation
            count = count + 1 ;
        end
        
        if time( count ) == count_times                                                       % all( )==1 , 即AX矩陣內元素都不為零
            
            if char( sat_sys( count ) ) == 80 &&  prn( count ) ~= 193 %&&  ismember( prn(count) ,GP_sat ) == 1
                sat_id( j ) = prn( count ) ;
                id_set_G( i , prn(count) ) = 111 ;
                CNR( j ) = cnr0( count ) ;
                sat_elevation( j ) = sat_EL( count ) ;
                sat_azimuth( j ) = sat_AZ( count ) ;
                j = j+1 ;
            elseif char( sat_sys( count ) ) == 66 %&& ismember( prn(count) ,BD_sat ) == 1
                sat_id_B( k ) = prn( count ) ;
                id_set_B( i , prn(count) ) = 111 ;
                CNR_B( k ) = cnr0( count ) ;
                sat_elevation_B( k ) = sat_EL( count ) ;
                sat_azimuth_B( k ) = sat_AZ( count ) ;
                k = k + 1 ;
            end
            count = count + 1 ;
            
        else
            count_times = time(count);
            judge=0;
        end
        
        if count>length(time)
            count=count-1;
            break
        end
        
    end

    %----------------每次(每一秒)定位衛星的顆數---------------
    nsat = length( sat_id ) ;
    nsat_B = length( sat_id_B ) ;
  
    %----------------------SkyPlot ( GPS )--------------------------------
    svx = zeros( 1 , nsat ) ;               svy = zeros( 1 , nsat ) ;
    svx_B = zeros( 1 , nsat_B ) ;      svy_B = zeros( 1 , nsat_B ) ;
    az_radians = sat_azimuth * ( pi/180) ;
    az_radians_B = sat_azimuth_B * ( pi/180) ;
    zenith = 90 - sat_elevation ;                                        %Convert elevation angle to zenith
    zenith_B = 90 - sat_elevation_B ;                                        %Convert elevation angle to zenith

    for s = 1 : nsat
        svx( s ) = zenith( s ) * cos( az_radians( s ) ) ; 
        svy( s ) = zenith( s ) * sin( az_radians( s ) ) ; %Calculate polar co-ordinates
    end
    
    for s = 1 : nsat_B
        svx_B( s ) = zenith_B( s ) * cos( az_radians_B( s ) ) ; 
        svy_B( s ) = zenith_B( s ) * sin( az_radians_B( s ) ) ; %Calculate polar co-ordinates
    end
    
    az_radians = [] ;   zenith = [] ;   az_radians_B = [] ;     zenith_B = [] ;
    figure(1)
    polarhg( [ 15 30  45 60 75 ] ) 
    hold on
    
    % appeared_id中為出現過的衛星編號,故第一秒時為零向量,則第一秒[tf,index]=ismember(sat_id , appeared_id )的
    % tf亦為零向量(0代表沒出現過),反操作後(~)變成1向量 → 變成沒出現過的output會是1
    % any:只要向量內有一個元素為非零(true),則傳回1
    if any(~ismember( sat_id , appeared_id )) == 1       % any內只要有一個衛星沒出現過(tf為1),則進到for loop進行text的編號動作
        for s = 1 : nsat
            if ismember( sat_id( s ) , appeared_id ) == 0     %一個一個取,沒出現過的才要編號 ; 出現過的衛星則不編號,掉進else後直接畫點
                plot( svx(s) , svy(s) , 'r.' , 'markersize', 20 )
                id_string = num2str( sat_id( s ) ) ;
                sv_id = [ 'G' id_string ] ;
                text( svx( s )+6 , svy( s ) , sv_id , 'FontSize' , 12 ,'color' , 'r' ) ;
                appeared_id( sat_id( s ) ) = sat_id( s ) ;     
                axis('square')
                grid on;  
            else
                plot( svx(s) , svy(s) , 'r.' , 'markersize', 1 )
            end
        end
        
    else
        plot( svx , svy , 'r.' , 'markersize', 1 )        % 如果都是出現過的衛星,直接向量化畫點,節省運算速度
    end
    
 
    if any(~ismember( sat_id_B , appeared_id_B )) == 1          
        for s = 1 : nsat_B
            if ismember( sat_id_B( s ) , appeared_id_B ) == 0
                plot( svx_B(s) , svy_B(s) , 'b*' , 'markersize', 10 )
                id_string_B = num2str( sat_id_B( s ) ) ;
                sv_id_B = [ 'B' id_string_B ] ;
                text( svx_B( s )+6 , svy_B( s ) , sv_id_B , 'FontSize' , 12 ,'color' , ' b' ) ;
                appeared_id_B( sat_id_B( s ) ) = sat_id_B( s ) ;     
                axis('square')
                grid on;  
            else
                plot( svx_B(s) , svy_B(s) , 'b.' , 'markersize', 1 )
            end
        end
        
    else
        plot( svx_B , svy_B , 'b.' , 'markersize', 1 )
    end
    
    if beautify == 1
        set(gcf, 'Color', 'w');                        %Change background of figure from grey to white
        ti = get(gca,'TightInset')   ;                 %Remove extra spacing around figure
        set(gca, 'LooseInset', [0,0,0,0.01]);          %Depending on the figure, you may need to add extra spacing [left bottom width height])        
        beautify = 0 ;
    end
    
    if i == runtime
        i
        Data
        print( '-dpng',  'Skyplot', '-r600') ;       %Change "-r600" to the required DPI
    end
    
    
end
toc