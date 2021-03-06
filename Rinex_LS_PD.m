tic
format long g
clc ;
clear ;
close all ;
addpath Functions

Data_List_rinex                  %   從Data_List_rinex.m內讀取資料

Preparation               %   前置作業
Data

%%
% error reference
Rhoc_error = [];
EL_error = [];
CNR_error = [];
%ref_pos_0=[-2419080.98069294     5379843.08493094    2418103.28729912];                 %mean HED
%ref_pos_0=[-2419080.39060001          5379846.01229999          2418100.60679998];                 %07102mean HED
%ref_pos_0=[ -2417142.72255327          5382343.78217749          2415036.61995793  ];      %mean hkst
%ref_pos_0=[ -2424425.1013  5377188.1768  2418617.7454    ];                                %hkss
%ref_pos_0=[ -3042356.01766021          4911053.45651187   2694094.56271902    ];                                %NTOU mean
%ref_pos_0 = lla2ecef([22.4259617417365 , 114.211355684335 ,53.2231195336208], 'WGS84' ) %ref


%%
first_P = 0 ;       %   用來記錄第一次成功定位的時刻數
bar = waitbar(0,'Please wait...');
%%
for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %每900秒,重新更新lagrange內插節點
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
            SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX為(prn*data_num)的矩陣,即 第一列為衛星prn1在13個時刻的x位置 ,第二列為衛星prn2在13個時刻的x位置
    hel_count=0;
    wl_count = 0 ;                                      %計算GPS/BDS權重疊代次數
    GP_id = [] ;
    Rho0_G = [] ;
    Rho0_Gpr = [] ; Ds =[];
    CNR_G = [] ;
    EL_G = [] ;
    AZ_G = [] ;
    Rg = zeros( 32,1 ) ;
    Rrg( :,1 ) = Rrg( :,2 ) ;
    Rrg( :,2 ) = 0 ;
    Carrier_G( :,1 ) = Carrier_G( :,2 ) ;
    Carrier_G( :,2 ) = 0 ;
    SVxyz = zeros( 1,3 ) ;
    AEL_G = [] ;            AAZ_G = [] ;
    ARr_G = [] ;            ARho0_G = [] ;
    ASCBg = [] ;            Asat_G = [] ;
    AGP_id = [] ;            ACNR_G = [] ;
    Atime_tx = [] ;
    ACNR_G = [] ;
    
    stEL_G = [] ;
    m = 1 ;
    l = 1 ;
    s = 1 ;
    first_in = 0 ;
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %每一秒都把各顆GPS的資料依序抓進來
    k = 1 ;                           %每一秒都把各顆BDS的資料依序抓進來
    
    while judge == 1
        
        
        
        if time( count ) == count_times                                                       % all( )==1 , 即AX矩陣內元素都不為零
            if prn( count ) < 2000 && prn( count ) ~= 1193 && all( AX( prn(count)-1000,: ) ) == 1 && AS( prn(count)-1000,1 ) ~= 999999.999999...
                    &&( prn( count ) == 1002  || prn( count ) == 1005  || prn( count ) == 1006 || prn( count ) == 1012|| prn( count ) == 1013 || prn( count ) == 1019|| prn( count ) == 1020 || prn( count ) == 1025|| prn( count ) == 1029);
                prn( count ) = prn( count ) -1000 ;
                
                Rg( prn(count) , 1 ) = pr( count ) ;
                Rrg( prn(count) , 2 ) = pr_rate( count )*lambda_G ;
                
                if rem(ADR( count ),0.001) ==0
                    Carrier_G( prn(count) , 2 ) = ADR( count )*lambda_G;
                    if  obs_flag == 3 && i <= initial_time || obs_flag <= 2
                        GP_id( j , 1) = prn( count ) ;
                        Rho0_Gpr( j , 1 ) = pr( count ) ;
                        Ds( j , 1 )  = pr_rate( count ) *lambda_G ;
                        CNR_G( j , 1) = cnr0( count ) ;
                        j = j+1 ;
                    elseif obs_flag == 3 && i > initial_time && count_G( prn(count),1 ) > converge_time && Carrier_G(prn(count),2) ~= 0 && Carrier_G(prn(count),1) ~= 0
                        GP_id( j , 1) = prn( count ) ;
                        CNR_G( j , 1) = cnr0( count ) ;
                        
                        j = j+1 ;
                    end
               
                end
                
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
    nGsat = length( GP_id ) ;
    nsat = nGsat ;
    
    %-----------------------------------------------------------------------------------------------
    if nsat ~= 0
        
        %-----------得到接收時間---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        
        %----------------選擇觀測量------------------
        if obs_flag == 1
            Rho0_G = Rho0_Gpr ;
        elseif obs_flag == 2
            [ Rho0_G , S_G , count_G ] = Doppler_Smoothed_Code_Rinex( Rrg , Rg , S_G , GP_id , count_G , GMax ) ;
        elseif obs_flag == 3
            [ Rho0_G , S_G , count_G ] = Carrier_Smoothed_Code( Carrier_G , Rg , S_G , GP_id , count_G , GMax ) ;
        end
        
        Rr_G = Rho0_G ;
        
        Tg = zeros( nGsat , 1 ) ;
        ig = zeros( nGsat , 1 ) ;
        
        %--------------------得到transmission time ( time_tx )----------------------
        [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
            satellite_positions( r_gpst , Rr_G , GP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
        
        %------------使用 time_tx 計算精密星歷位置-----------------
        %--------------initialize----------------
        sat_G = zeros( nGsat , 3 ) ;
        SCBg = zeros( nGsat , 1 ) ;
        
        %--------------------計算精密星歷衛星位置及鐘差---------------------------
        for k = 1 : nGsat
            sat_G( k , 1 ) = LagrangeInter( Data_time , AX( GP_id(k),: ) , time_tx( k ) ) ;
            sat_G( k , 2 ) = LagrangeInter( Data_time , AY( GP_id(k),: ) , time_tx( k ) ) ;
            sat_G( k , 3 ) = LagrangeInter( Data_time , AZ( GP_id(k),: ) , time_tx( k ) ) ;
            for t = 1 : length( AX( GP_id(k),: ) )
                if Data_time( t ) - time_tx( k ) > 0
                    Interval = t-1 : t ;
                    break
                end
            end
            ephemeris_SCB = AS( GP_id(k),: ) ; %SCB=satellite clock bias
            SCBg( k , 1 ) = LagrangeInter( Data_time( Interval ) , ephemeris_SCB( Interval ) , time_tx( k ) ) * c * 1e-6 ;
        end
        
        %-------------------計算仰角----------------------------------
        for a= 1 : nGsat
            [ EL_G( a , 1 ) , AZ_G( a,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , sat_G( a , : ) ) ;
            %if EL_G( a , 1 ) > smaller_than_elevation
            AEL_G( l , 1) = EL_G( a , 1 ) ;
            AAZ_G( l , 1) = AZ_G( a , 1 ) ;
            ARr_G( l , 1) = Rr_G( a , 1 ) ;
            ARho0_G( l , 1) = ARr_G( l , 1) ;
            ASCBg( l , 1) = SCBg( a , 1 ) ;
            Asat_G( l , :) = sat_G(a , :) ;
            AGP_id( l , 1) = GP_id(a , 1) ;
            ACNR_G( l , 1) = CNR_G(a , 1) ;
            Atime_tx( l , 1) = time_tx( a , 1) ;
            runtime_GEL(AGP_id(l,1),i)=AEL_G( l , 1);
            l = l + 1 ;
            %else
            %stEL_G(s,1)= GP_id( a , 1 ) ;
            %s = s + 1 ;
            %end
        end
        
        AnGsat = length( AGP_id ) ;
        Ansat = AnGsat  ;
        
        %-------------------計算covariance----------------------------
        EL =  AEL_G  ;
        %CNR = zeros(Ansat,1) ;
        CNR = [ CNR_G  ] ;
        [ sat_var_elevation_all , sat_var_CNR_all , sat_var_SIGMA_all , sat_var_CandE_all ] = weightingfunc( EL , CNR ) ;
        
        %-------------------選擇權重矩陣------------------------------
        W = diag( ( 1./( sat_var_SIGMA_all ) ).^2 ) ;
        %W = eye( Ansat ) ;

        if AnGsat ~=0
            
            Delta_Mtrix( 1:4 ) = 1 ;
            
            while norm( Delta_Mtrix(1:3) ) > 0.01
                
                llh_0 = xyz2llh( rec_pos_0 ) ;
                lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
                lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
                height_rec = llh_0( 3 ) ;                  % unit : meter
                
                Tg = zeros( AnGsat , 1 ) ;
                ig = zeros( AnGsat , 1 ) ;
                
                %--------------------得到transmission time ( time_tx )----------------------
                [ XS , dtS , XS_tx , VS_tx , Atime_tx , no_eph , sys_idx ] = ...                 %dtS已包含相對論誤差
                    satellite_positions( r_gpst , ARr_G , AGP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
                
                Vs =  VS_tx ;
                
                %------------利用 UNB3 模型計算對流層誤差-------------------------
                for n = 1 : AnGsat
                    [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , AEL_G(n) * pi/180 ) ;
                    Tg( n , 1 ) = RTROP ;
                end
                
                %------------使用GIM計算電離層誤差--------------------------
                for n = 1 : AnGsat
                    [ iono , VTEC ]=...
                        GIM_corr( AAZ_G(n) , AEL_G(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                    ig( n , 1 ) = iono ;
                end
                
                %------------使用 time_tx 計算精密星歷位置-----------------
                %--------------initialize----------------
                sat_G = zeros( AnGsat , 3 ) ;
                SCBg = zeros( AnGsat , 1 ) ;
                
                %--------------------計算精密星歷衛星位置及鐘差---------------------------
                for k = 1 : AnGsat
                    sat_G( k , 1 ) = LagrangeInter( Data_time , AX( AGP_id(k),: ) , Atime_tx( k ) ) ;
                    sat_G( k , 2 ) = LagrangeInter( Data_time , AY( AGP_id(k),: ) , Atime_tx( k ) ) ;
                    sat_G( k , 3 ) = LagrangeInter( Data_time , AZ( AGP_id(k),: ) , Atime_tx( k ) ) ;
                    for t = 1 : length( AX( AGP_id(k),: ) )
                        if Data_time( t ) - Atime_tx( k ) > 0
                            Interval = t-1 : t ;
                            break
                        end
                    end
                    ephemeris_SCB = AS( AGP_id(k),: ) ; %SCB=satellite clock bias
                    SCBg( k , 1 ) = LagrangeInter( Data_time( Interval ) , ephemeris_SCB( Interval ) , Atime_tx( k ) ) * c * 1e-6 ;
                end
                
                %------------------計算GPS相對論誤差-------------------------------
                relative_IGS = zeros( AnGsat , 1 ) ;
                tgd = 0 ;
                
                for n = 1 : AnGsat
                    icol = find_eph( Eph_tatol , AGP_id( n ) , r_gpst ) ;
                    Eph = Eph_tatol( : , icol ) ;
                    
                    [ satp , satv ] = satellite_orbits( Atime_tx( n ) , Eph , icol , [] ) ;
                    
                    Sp3_corr = -2*dot( satp , satv )/c ;
                    relative_IGS( n , 1 ) = Sp3_corr - tgd*c ;
                end
                
                %--------------------修正虛擬距離(各項單位皆為公尺)--------------------------
                % 北斗的SCBb包含sat_bias,relative與group_delay
                Rhoc_G = ARho0_G - ig - Tg + SCBg + relative_IGS ;
                Rhoc = [ Rhoc_G ];
                
                %---------------------------修正地球自轉--------------------------------
                [ Xs_G , ARr_G, Vs ] = Fix_Earth_Rotation_Ds( sat_G , rec_pos_0, Vs ) ;
                Xs = [ Xs_G ] ;
                
                %-----------------------最小平方法定位-------------------------
                if Ansat > 3
                    
                    [ Delta_Mtrix , Delta_Rho ] = LeastSquare_PD( Xs , Vs  , rec_pos_0 , Rhoc,Ds , W , Ansat , rec_bias) ;
                    
                    rec_pos = rec_pos_0 + Delta_Mtrix(1:3)' ;
                    rec_bias = rec_bias + Delta_Mtrix(4) ;
                    rec_pos_0 = rec_pos ;
                else
                    break
                end
                if hel_count >= 100
                    i
                    norm( Delta_Mtrix(1:3) )
                    AGP_id
                    break
                end
                                %}
            end
            
            if Ansat < 4
                Xxyz( i , : ) = rec_pos_0 ;
            else
                Xxyz( i , : ) = rec_pos ;
            end
            %---------------count DOP--------------------
            [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
            rec_pos_0 = Xxyz( i , : ) ;            %   將所得之位置更新為下一時刻之初始位置

            %--------------將低仰角衛星之CSC資料清除----------
            %{
            if obs_flag > 1
                for q = 1 : size(stEL_G)
                    S_G(stEL_G(q),1) = 0 ;
                    count_G(stEL_G(q),1) = 0 ;
                end
                for q = 1 : size(stEL_B)
                    S_G(stEL_B(q),1) = 0 ;
                    count_B(stEL_B(q),1) = 0 ;
                end
            end
            %}
        end
%{
        %------------紀錄carrier phase&true range-----------------
        if i > initial_time
            [ sss , ARr_G ] = Fix_Earth_Rotation( sat_G , ref_pos_0 ) ;
            Rhoc_error = [ Rhoc_error ; (Rhoc - ARr_G-rec_bias) ];
            EL_error =[EL_error ; EL] ;
            CNR_error = [CNR_error ;CNR] ;
        end
        %}
    end
    
    %----------------時間+1-------------------------
    time_week_last_epoch = time_week ;
    time_week = time( count ) ;
    
    time_day_UTC = time_day_UTC + 1 ;
    local_time = local_time + 1 ;
    
    %---------------------計算固體潮汐誤差(solid tide)----------------
    if Data ~= no_solid_tide
        e_solid_error = interp1( time_solid , e_solid , time_day_UTC ) ;
        n_solid_error = interp1( time_solid , n_solid , time_day_UTC ) ;
        u_solid_error = interp1( time_solid , u_solid , time_day_UTC ) ;
        solid_tide_error_temp = [ e_solid_error , n_solid_error , u_solid_error ] ;
        solid_tide_error( i , : ) = solid_tide_error_temp ;
    end
    
    if i > 1 && Ansat ~= nsat_prior
        nsat_change(chg) = i ;
        chg = chg + 1 ;
    end
    
    nsat_prior = nsat ; %   紀錄上一秒衛星顆數
    nGsat_Mtrix( i ) = AnGsat ;
    nsat_Mtrix( i ) = Ansat ;
    if first_P == 0 && all( Xxyz( i,: ) ) == 1  %第一次定位成功後,紀錄當下時刻
        first_P = i ;
    end
    
end

close( bar ) ;
switch obs_flag
    case { 0 }
        obs = 'Dual_EPR_LS' ;
    case { 1 }
        obs = 'Dual_PR_LS' ;
    case { 2 }
        obs = 'Dual_DSC_LS' ;
    case { 3 }
        obs = 'Dual_CSC_LS' ;
end

%---------------------資料取樣時間----------------------
if obs_flag ==3
    CutTimeArray = 1 : cut_time ;
    Xxyz(CutTimeArray , :) = [] ;
    nsat_Mtrix(CutTimeArray) = [] ;
    PDOP(CutTimeArray) = [] ;
    sampling_time_end = runtime - cut_time ;
    TimeArray = first_P : sampling_time_end ;
    %effTimeArray_error = round(length(TimeArray)/100) -1 + first_P : sampling_time_end ;
    effTimeArray = TimeArray;
elseif obs_flag == 1
    sampling_time_end = runtime ;
    effTimeArray = 1 : sampling_time_end ;
end
%Rhoc_error = abs(Rhoc_error);
%{
figure( 151 ) ;
scatter( CNR_error , Rhoc_error,'.' ) ;
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'error(m)' ) ;
grid on;                        hold on;
title('CNR-CSC error');
print( '-dpng',  'CNR-CSC error' , '-r600' ) ;

figure( 152 ) ;
scatter( EL_error , Rhoc_error,'.' ) ;
grid on ;
xlabel( 'EL(deg)' ) ;
ylabel( 'error(m)' ) ;
grid on;                        hold on;
title('EL-CSC error');
print( '-dpng',  'EL-CSC error' , '-r600' ) ;
%}

%{
CNR_error = round(CNR_error) ;
[CNR_error , ic] = sort(CNR_error) ;
CRhoc_error = Rhoc_error(ic) ;

EL_error = round(EL_error) ;
[EL_error , id] = sort(EL_error) ;
ERhoc_error = Rhoc_error(id) ;

% csvwrite('0710error.csv',[CNR_error  Rhoc_error  EL_error   ]);
csvwrite('071011hr_error.csv',[CNR_error  Rhoc_error  EL_error   ]);
 %Plotting_Rinex
 %%

cref = CNR_error(1) ;
Eref = EL_error(1) ;

ccalculate = [] ;
Ecalculate = [] ;
n=1;x_matrix = [] ;mean_matrix = [] ;std_matrix = [];
ss = 1 ;Ex_matrix = [] ;Emean_matrix = [] ;Estd_matrix = [];
for i = 1 : length(CNR_error)
    if CNR_error(i)~= cref
        temp_mean = mean(ccalculate) ;
        temp_std = std(ccalculate);
        x_matrix = [x_matrix ; cref];
        mean_matrix = [mean_matrix ; temp_mean];
        std_matrix = [std_matrix ; temp_std];
        cref =  CNR_error(i);
        ccalculate = [];
        clear temp_mean temp_std ;
        n=1;
    end
      if EL_error(i)~= Eref
        temp_mean = mean(Ecalculate) ;
        temp_std = std(Ecalculate);
        Ex_matrix = [Ex_matrix ; Eref];
        Emean_matrix = [Emean_matrix ; temp_mean];
        Estd_matrix = [Estd_matrix ; temp_std];
        Eref =  EL_error(i);
        Ecalculate = [];
        clear temp_mean temp_std ;
        ss=1;
    end
    ccalculate(n ,1) = CRhoc_error(i) ;
    Ecalculate(ss,1) = ERhoc_error(i) ;
    n = n+1;
    ss = ss + 1 ;
    
    
        if i==length(CNR_error)
        temp_mean = mean(ccalculate) ;
        temp_std = std(ccalculate);
        x_matrix = [x_matrix ; cref];
        mean_matrix = [mean_matrix ; temp_mean];
        std_matrix = [std_matrix ; temp_std];
        cref =  CNR_error(i);
        ccalculate = [];
        clear temp_mean temp_std ;
        n=1;
    end
      if i==length(CNR_error)
        temp_mean = mean(Ecalculate) ;
        temp_std = std(Ecalculate);
        Ex_matrix = [Ex_matrix ; Eref];
        Emean_matrix = [Emean_matrix ; temp_mean];
        Estd_matrix = [Estd_matrix ; temp_std];
        Eref =  EL_error(i);
        Ecalculate = [];
        clear temp_mean temp_std ;
        ss=1;
    end
end

figure( 11 ) ;
scatter(CNR_error , CRhoc_error );
grid on ;
xlabel( 'CNR' ) ;
ylabel( 'covariance error(m)' ) ;
grid on;                        hold on;
title('CNR-CSC error');
print( '-dpng',  'CNR-CSC error' , '-r600' ) ;

figure( 12 ) ;
scatter(EL_error , ERhoc_error );
grid on ;
xlabel( 'EL(Deg)' ) ;
ylabel( 'covariance error(m)' ) ;
grid on;                        hold on;
title('EL-CSC error');
print( '-dpng',  'EL-CSC error' , '-r600' ) ;

%}



toc