tic
format long g
clc ;
clear ;
close all ;
addpath Functions

Data_List_rinex                  %   �qData_List_rinex.m��Ū�����

Preparation               %   �e�m�@�~
Data


%----------------KF setting--------------------------
T = 1; % positioning interval
% State vector is as [x Vx y Vy z Vz b d].', i.e. the coordinate (x,y,z),
% the clock bias b, and their derivatives.

% Kalman Parameters initial setting
sigma= 0.1 ;      Rhoerror = 3 ;       Pvalue = 10 ;
%Sb = (1.1e-19)*c^2 ;        Sd = (4.3e-20)*c^2 ;
%Sb = (4e-19)*c^2 ;        Sd = ((pi^2)*16e-20)*c^2 ;%The study of GPS Time Transfer Based on Extended Kalman Filter
Sb = (1.1e-19)*c^2 ;          Sd = ((pi^2)*56e-21)*c^2 ;%Single-frequency, single-receiver terristrial and spaceborne point positioning(p.67)
%Sb = (4e-19)*c^2 ;        Sd = (1.58e-18)*c^2 ;

Xp = zeros( 8 , 1 ) ;
Pp = eye( 8 )*Pvalue ;
Xu = zeros( 8 , 1 ) ;
Pu_prior = [];



Qc = [Sb*T+Sd*T*T*T/3     Sd*T*T/2 ;
    Sd*T*T/2                            Sd*T ] ;
Qxyz = sigma^2 * [T^3/3      T^2/2 ;
    T^2/2        T      ] ;
Q = blkdiag( Qxyz , Qxyz , Qxyz , Qc ) ;

ff = [ 1 T ; 0 1 ] ;
fy = blkdiag( ff , ff , ff , ff  ) ;

Rcovariance = zeros( 32,runtime - cut_time ) ;
first_P = 0 ;
KF = 0 ;
IDcount = zeros(32 , 1);
Vac = [] ;
IDac = [] ;
NSATac = [] ;
Vxac = [] ;
Pu = [] ;
id_comn = [] ;
Vc_matrix=zeros(runtime,1);
%[ CNRref ,CNRmean ,CNRstd ] =textread( '0710CSCref.csv' , '%d %f %f ' , 'delimiter' , ',' ) ;%headerlines:���L�Ĥ@��
%%
%KALman static setting
Ppstatic= eye( 4 )*Pvalue ;
Xpstatic = zeros(4,1);
Xustatic = zeros(4,1);
Qstatic=blkdiag(Sb*T+Sd*T*T*T/3 ,Sb*T+Sd*T*T*T/3,Sb*T+Sd*T*T*T/3,Sb*T+Sd*T*T*T/3);
fystatic = blkdiag( 1 , 1 , 1 , 1  ) ;

%%
bar = waitbar(0,'Please wait...');
for i = 1 : runtime
    
    str=['Positinoing ',num2str(i),'s'];
    waitbar( i/runtime , bar , str )   % computation here %
    
    if rem( i , 900 ) == 0  %�C900��,���s��slagrange�����`�I
        [ SP3_data_num , SP3_data , Data_time , DesiretimeS , DesiretimeE , AX , AY , AZ , AS ]  = ...
            SP3_read( year , month , day , SP3_filename , week , time_week , interpolation_order ) ;
    end   %AX��(prn*data_num)���x�},�Y �Ĥ@�C���ìPprn1�b13�Ӯɨ誺x��m ,�ĤG�C���ìPprn2�b13�Ӯɨ誺x��m
    hel_count=0;
    wl_count = 0 ;                                      %�p��GPS/BDS�v���|�N����
    GP_id = [] ;
    Rho0_G = [] ;
    Rho0_Gpr = [] ;
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
    
    
    stEL_G = [] ;
    m = 1 ;
    l = 1 ;
    s = 1 ;
    %-----------------------------------------------
    judge = 1 ;
    j = 1 ;                            %�C�@����U��GPS����ƨ̧ǧ�i��
    k = 1 ;                           %�C�@����U��BDS����ƨ̧ǧ�i��
    
    while judge == 1
        
        if time( count ) == count_times                                                       % all( )==1 , �YAX�x�}�������������s
            if prn( count ) < 2000 && prn( count ) ~= 1193 && all( AX( prn(count)-1000,: ) ) == 1 && AS( prn(count)-1000,1 ) ~= 999999.999999...
                    &&( prn( count ) == 1003  || prn( count ) == 1005  || prn( count ) == 1006 || prn( count ) == 1012|| prn( count ) == 1013 || prn( count ) == 1019|| prn( count ) == 1020 || prn( count ) == 1025|| prn( count ) == 1029);
                
                prn( count ) = prn( count ) -1000 ;
                
                Rg( prn(count) , 1 ) = pr( count ) ;
                Rrg( prn(count) , 2 ) = pr_rate( count )*lambda_G ;
                if rem(ADR( count ),0.001) ==0
                    Carrier_G( prn(count) , 2 ) = ADR( count )*lambda_G;
                    checkmatrix( prn(count),1 ) = 1 ;
                    checkmatrix( prn(count),3 ) = checkmatrix( prn(count),3 ) + 1 ;
                    if  obs_flag == 3 && i <= initial_time || obs_flag <= 2
                        GP_id( j , 1) = prn( count ) ;
                        Rho0_Gpr( j , 1 ) = pr( count ) ;
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
    
    %----------------�C��(�C�@��)�w��ìP������---------------
    nGsat = length( GP_id ) ;
    nsat = nGsat ;
    
    %-----------------------------------------------------------------------------------------------
    if nsat ~= 0
        
        %-----------�o�챵���ɶ�---------------------
        r_gpst = weektow2time ( week , time_week , 'G' ) ;  %GPST
        
        %----------------����[���q------------------
        if obs_flag == 1
            Rho0_G = Rho0_Gpr ;
        elseif obs_flag == 2
            [ Rho0_G , S_G , count_G ] = Doppler_Smoothed_Code( Rrg , Rg , S_G , GP_id , count_G , GMax ) ;
        elseif obs_flag == 3
            %[ Carrier_G  , RMSE , mean_delta , ht ] = DCDRM( Carrier_G , Rrg , lambda ,  RMSE , mean_delta , ht , scale_detect_cycle_slip) ;  %check cycle slip
            [ Rho0_G , S_G , count_G ] = Carrier_Smoothed_Code( Carrier_G , Rg , S_G , GP_id , count_G , GMax ) ;
        end
        
        %-------------------���ƶ���----------------------------------
        [ GP_id , id_index ]=sortrows( GP_id ) ;
        Rho0_G = Rho0_G( id_index ) ;
        Rr_G = Rho0_G ;
        CNR_G = CNR_G( id_index ) ;
        
        Tg = zeros( nGsat , 1 ) ;
        ig = zeros( nGsat , 1 ) ;
        
        %--------------------�o��transmission time ( time_tx )----------------------
        [ XS , dtS , XS_tx , VS_tx , time_tx , no_eph , sys_idx ] = ...                 %dtS�w�]�t�۹�׻~�t
            satellite_positions( r_gpst , Rr_G , GP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
        %------------�ϥ� time_tx �p���K�P����m-----------------
        %--------------initialize----------------
        sat_G = zeros( nGsat , 3 ) ;
        
        %--------------------�p���K�P���ìP��m�����t---------------------------
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
        end
        
        %-------------------�p�����----------------------------------
        for a= 1 : nGsat
            [ EL_G( a , 1 ) , AZ_G( a,1) ] = Calc_Azimuth_Elevation( rec_pos_0 , sat_G( a , : ) ) ;
            %if EL_G( a , 1 ) > smaller_than_elevation
            AEL_G( l , 1) = EL_G( a , 1 ) ;
            AAZ_G( l , 1) = AZ_G( a , 1 ) ;
            ARr_G( l , 1) = Rr_G( a , 1 ) ;
            ARho0_G( l , 1) = ARr_G( l , 1) ;
            Asat_G( l , :) = sat_G(a , :) ;
            AGP_id( l , 1) = GP_id(a , 1) ;
            Atime_tx( l , 1) = time_tx( a , 1) ;
            ACNR_G( l , 1) = CNR_G( a , 1 ) ;
            l = l + 1 ;
            %else
            %stEL_G(s,1)= GP_id( a , 1 ) ;
            %s = s + 1 ;
            %end
        end
        
        AnGsat = length( AGP_id ) ;
        Ansat = AnGsat  ;
        %-------------------�p��covariance----------------------------
        EL =  AEL_G  ;
        %CNR = zeros(Ansat,1) ;
        CNR = [ ACNR_G  ] ;
        [ sat_var_elevation_all , sat_var_CNR_all , sat_var_SIGMA_all , sat_var_CandE_all ] = weightingfunc( EL , CNR ) ;
        
        %-------------------����v���x�}------------------------------
        W = diag( ( 1./( sat_var_SIGMA_all ) ).^2 ) ;
        %W = eye( Ansat ) ;
        if obs_flag == 3 && chenckcount ~=0 && i > initial_time
            [W , checkmatrix ] = adjweightingfunc( W , checkmatrix , AnGsat , chenckcount ) ;
        end
        if AnGsat ~=0
            
            Delta_Mtrix( 1:4 ) = 1 ;
            
            while norm( Delta_Mtrix(1:3) ) > 0.01
                
                llh_0 = xyz2llh( rec_pos_0 ) ;
                lat_rec = llh_0( 1 ) * 180/pi ;          % unit : degree
                lon_rec = llh_0( 2 ) * 180/pi ;         % unit : degree
                height_rec = llh_0( 3 ) ;                  % unit : meter
                
                Tg = zeros( AnGsat , 1 ) ;
                ig = zeros( AnGsat , 1 ) ;
                
                %--------------------�o��transmission time ( time_tx )----------------------
                [ XS , dtS , XS_tx , VS_tx , Atime_tx , no_eph , sys_idx ] = ...                 %dtS�w�]�t�۹�׻~�t
                    satellite_positions( r_gpst , ARr_G , AGP_id , Eph_tatol , [] , [] , Tg , ig , rec_bias_satpos ) ;
                
                %------------�Q�� UNB3 �ҫ��p���y�h�~�t-------------------------
                for n = 1 : AnGsat
                    [ RTROP ,HZD , HMF , WZD , WMF ] = UNB3M( lat_rec * pi/180 , height_rec , day_of_year , AEL_G(n) * pi/180 ) ;
                    Tg( n , 1 ) = RTROP ;
                end
                
                %------------�ϥ�GIM�p��q���h�~�t--------------------------
                for n = 1 : AnGsat
                    [ iono , VTEC ]=...
                        GIM_corr( AAZ_G(n) , AEL_G(n) , time_week , lat_rec , lon_rec , VTEC_MAP , time_total , lat_total , lon_total ) ;
                    ig( n , 1 ) = iono ;
                end
                
                %------------�ϥ� time_tx �p���K�P����m-----------------
                %--------------initialize----------------
                sat_G = zeros( AnGsat , 3 ) ;
                SCBg = zeros( AnGsat , 1 ) ;
                
                %--------------------�p���K�P���ìP��m�����t---------------------------
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
                
                %------------------�p��GPS�۹�׻~�t-------------------------------
                relative_IGS = zeros( AnGsat , 1 ) ;
                tgd = 0 ;
                
                for n = 1 : AnGsat
                    icol = find_eph( Eph_tatol , AGP_id( n ) , r_gpst ) ;
                    Eph = Eph_tatol( : , icol ) ;
                    
                    [ satp , satv ] = satellite_orbits( Atime_tx( n ) , Eph , icol , [] ) ;
                    
                    Sp3_corr = -2*dot( satp , satv )/c ;
                    relative_IGS( n , 1 ) = Sp3_corr - tgd*c ;
                end
                
                %--------------------�ץ������Z��(�U�����Ҭ�����)--------------------------
                % �_�檺SCBb�]�tsat_bias,relative�Pgroup_delay
                Rhoc_G = ARho0_G - ig - Tg + SCBg + relative_IGS ;
                Rhoc = [ Rhoc_G ] ;
                
                %---------------------------�ץ��a�y����--------------------------------
                [ Xs_G , ARr_G ] = Fix_Earth_Rotation( sat_G , rec_pos_0 ) ;
                Xs = [ Xs_G ] ;
                
                %-----------------------�̤p����k�w��-------------------------
                if Ansat > 3 %&& i <= initial_time
                    
                    [ Delta_Mtrix , Delta_Rho ] = LeastSquare( Xs , rec_pos_0 , Rhoc , W , Ansat , rec_bias ) ;
                    
                    rec_pos_0 = rec_pos_0 + Delta_Mtrix(1:3)' ;
                    rec_bias = rec_bias + Delta_Mtrix(4) ;
                else
                    break
                end
                if hel_count >= 100
                    i
                    norm( Delta_Mtrix(1:3) )
                    AGP_id
                    
                    break
                end
            end
            
            if Ansat < 4
                Xxyz( i , : ) = rec_pos_0 ;
            elseif Ansat>4 %&& i <= initial_time
                Xxyz( i , : ) = rec_pos_0 ;
                Xp([1,3,5]) = Delta_Mtrix(1:3).' ;
                Xp(7) = Delta_Mtrix(4) ;                        
                if KF ==0 && i==initial_time
                    ref_pos = rec_pos_0 ;
                   ref_bias =  rec_bias;
                    KF = 1 ;
                end

            end
            if KF == 1 && i > initial_time % && any( Xxyz(i-1,:) ) == 1
                
                %--------------------Kalman Filter------------------------------
                [ Val , H ] = H_Mtrix_Single( ref_pos , Xs , Ansat  , ref_bias ) ;
                Z = Rhoc - Val ;    % Z= z - zp
                R =  eye( size( Xs , 1 ) ) * Rhoerror;
                %{
                if ~isempty(id_comn) && ~isequal(id_comn , AGP_id)
                    [common , ing , ic ] = intersect(AGP_id,id_comn);
                    R_diag = [] ;
                    R_diag(1:Ansat) = Rhoerror ;
                    R_diag(ing) = R0(ic) ;
                    R = diag( R_diag ) ;
                elseif ~isempty(id_comn) && isequal(id_comn , AGP_id)
                    R = diag(R0) ;
                else
                    R =  eye( size( Xs , 1 ) ) * Rhoerror;
                end
                %}
                
                Vc_matrix(i)=trace(Vc*Vc');
                Vc_delta(i)=trace(Vc*Vc')-Vc_matrix(i-1);
                %{
                if Vc_matrix(i-1)~=0 && Vc_matrix(i-2)~=0
                    if  abs(Vc_delta(i))>=13
                        innovation_window = 5;
                        
                    elseif abs(Vc_delta(i))< 13 && abs(Vc_delta(i))>=10
                        innovation_window = 15;
                        
                    elseif abs(Vc_delta(i)) <10 && abs(Vc_delta(i))>=7
                        innovation_window = 25;
                        
                    elseif abs(Vc_delta(i)) <7 && abs(Vc_delta(i))>=4
                        innovation_window = 30;
                        
                    elseif abs(Vc_delta(i)) <4 && abs(Vc_delta(i))>=1
                        innovation_window = 35;
                    
                    elseif abs(Vc_delta(i)) <=1
                        innovation_window = 40;                       
                    end
                end
                %}
                K=(inv(Pp)+H' * inv(R )* H)\H'*inv(R);
                %K = Pp * H' / (H * Pp * H.' + R) ;
                Xu = Xp + K * ( Z - (H*Xp) ) ;        %   Z���ץ��L�᪺pseudorange
                Vc = Z - H*Xp;       %innovation sequence                
                %[  R , IDcount , IDac , NSATac , Vac , id_comn ] =IAE_change( AGP_id , Vc , H , Pp , Vac , IDac , NSATac , innovation_window , IDcount,R ) ;
                
                Vx = Xu - Xp ;        % state correction matrix
                
                
                I = eye( size(Xu,1),size(Xu,1) ) ;
                Pu = (I - K * H) * Pp;
                
                if i > cut_time
                    Rcovariance(AGP_id , jj) = diag(R);
                    jj=jj+1 ;
                end
                [ Q , Vxac ,Pu_prior] = Adaptive_state_co_V2v( Vx , Vxac , innovation_window , Q, Pu, Pu_prior);%adaptive Q
                
                Pp =  fy * Pu * fy.' + Q ;
                Xp = fy * Xu ;                  % fy��State transition matrix
                
                Xpstatic(1:4,1)=Xu( [1,3,5,7] ) ;
                Xustatic=Xpstatic;
                
                %%
                %{
                %static Kalman filter loop
                Delta_Xu( 1:4 ) = 1 ;
                while norm( Delta_Xu(1:3) ) > 0.001
                    
                    %--------------------Kalman Filter------------------------------
                    [ Val , H ] = H_Mtrix_Single_new( ref_pos , Xs , Ansat  , ref_bias ) ;
                    
                    
                    
                    K = Ppstatic * H' / (H * Ppstatic * H.' + R) ;
                    Xustatic = Xpstatic + K * ( Z - (H*Xpstatic) ) ;        %   Z���ץ��L�᪺pseudorange
                    
                    I = eye( size(Xustatic,1),size(Xustatic,1) ) ;
                    Pustatic = (I - K * H) * Ppstatic;
                    
                    Ppstatic =  fystatic * Pustatic * fystatic.' + Qstatic ;
                    
                    Delta_Xu(1:4)= Xustatic-Xpstatic;
                    
                    Xpstatic = fystatic * Xustatic ;                  % fy��State transition matrix
                    
                end
                %}
                Xxyz( i , : ) = ref_pos + Xu( [1:3] )' ;
                %Xp([1,3,5,7]) = Xustatic( [1:4] )' ;
            end
            
            %---------------count DOP--------------------
            [ EDOP(i) , NDOP(i) , VDOP(i) , HDOP(i) , PDOP(i) , GDOP(i) ] = Calculate_DOP( Xs , Xxyz( i , : ) ) ;
            rec_pos_0 = Xxyz( i , : ) ;            %   �N�ұo����m��s���U�@�ɨ褧��l��m
            %rec_pos_01 = rec_pos1;          %   �N�ұo����m��s���U�@�ɨ褧��l��m
            %--------------�N�C�����ìP��CSC��ƲM��----------
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
    end
    %{
    %------------����carrier phase&true range-----------------
    for r=1:AnGsat
        run_ADR(GP_id(r),i) = Carrier_G( GP_id(r) , 2 ) + ig(r , 1) - Tg(r , 1) + SCBg(r , 1) + relative_IGS(r , 1) - rec_bias ;
        run_trueRange(GP_id(r),i) = ARr_G(r) ;
    end
    run_recbias(i)= rec_bias;
    %}
    %----------------�ɶ�+1-------------------------
    time_week_last_epoch = time_week ;
    time_week = time( count ) ;
    
    time_day_UTC = time_day_UTC + 1 ;
    local_time = local_time + 1 ;
    
    %---------------------�p��T�����~�t(solid tide)----------------
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
    
    nsat_prior = nsat ; %   �����W�@��ìP����
    nGsat_Mtrix( i ) = AnGsat ;
    nsat_Mtrix( i ) = Ansat ;
    if first_P == 0 && all( Xxyz( i,: ) ) == 1  %�Ĥ@���w�즨�\��,������U�ɨ�
        first_P = i ;
        jj= 1 ;
    end
    
end

close( bar ) ;
%run_delta = run_ADR-run_trueRange;
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

%---------------------��ƨ��ˮɶ�----------------------
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


figure( 13 ) ;
plot( 70:length(Vc_delta) , Vc_delta(70:length(Vc_delta)) ) ;
grid on ;

Plotting_Rinex

toc