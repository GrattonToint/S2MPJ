function varargout = PORTFL3(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PORTFL3
%    *********
% 
%    Establish the sensitivity of certain obligations combinations in portfolio
%    analysis, based on the methodology of W. Sharpe.
% 
%    The asset allocation (the allocation of an investor's portfolio across a 
%    number of major asset classes) accounts for a large part of the variability
%    in the return on a typical investor's portfolio. It's important to
%     determine 
%    the exposures of each component of an investor's overall portofolio to
%    movements in their returns. The purpose of this problem is to determine the
%    conscious or inconscious style developped by the manager in it's fund and 
%    thereafter establish if the manager has been able to add value to his style
%    by security selection or active country allocation process. In our case the 
%    factor 1 to 12 represent the performance of specific asset classes
%     influencing 
%    the return of the fund: 12 JP Morgan country factors.
% 
%    Realized fund returns to establish the typical exposures of
%    the fund to the following asset classes: BBL International, BBL HY, 
%    CG MULTI International, CG MULTI HY, G BOND RENTINPLUS, PANELFUND HY.
% 
%    The factor's coefficients are required to sum to 1 with the constraint
%    of positive coefficients
% 
%    Source: 
%    DATASTREAM   Period: 15.1.91 to 15.3.96 (62 observations),
%    collected by D. Baudoux, July 1996.
% 
%    SIF input: Ph. Toint, July 1996.
% 
%    classification = 'C-CSLR2-MN-12-1'
% 
%    number of sensitivities
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 21 VI 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PORTFL3';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('NS') = 12;
        v_('NR') = 62;
        v_('F1,1') = 0.02802144;
        v_('F1,2') = -0.000237;
        v_('F1,3') = 0.01635846;
        v_('F1,4') = 0.01119664;
        v_('F1,5') = -0.0063052;
        v_('F1,6') = 0.00162501;
        v_('F1,7') = 0.01506489;
        v_('F1,8') = 0.01248192;
        v_('F1,9') = 0.00819364;
        v_('F1,10') = 0.00469729;
        v_('F1,11') = 0.0141744;
        v_('F1,12') = 0.02114737;
        v_('F1,13') = 0.00566105;
        v_('F1,14') = -0.0045603;
        v_('F1,15') = 0.01052255;
        v_('F1,16') = 0.00701282;
        v_('F1,17') = 0.00499437;
        v_('F1,18') = 0.0;
        v_('F1,19') = 0.00524953;
        v_('F1,20') = 0.02764239;
        v_('F1,21') = 0.0328613;
        v_('F1,22') = 0.01121753;
        v_('F1,23') = 0.01998054;
        v_('F1,24') = 0.01507346;
        v_('F1,25') = 0.00545113;
        v_('F1,26') = 0.01265034;
        v_('F1,27') = 0.00086154;
        v_('F1,28') = 0.01629365;
        v_('F1,29') = 0.02674088;
        v_('F1,30') = -0.0028283;
        v_('F1,31') = 0.00608639;
        v_('F1,32') = -0.0120404;
        v_('F1,33') = 0.0242554;
        v_('F1,34') = 0.02048871;
        v_('F1,35') = 0.02713002;
        v_('F1,36') = 0.00382081;
        v_('F1,37') = -0.0174316;
        v_('F1,38') = -0.0057826;
        v_('F1,39') = -0.0075103;
        v_('F1,40') = -0.0137687;
        v_('F1,41') = 0.00490366;
        v_('F1,42') = 0.00430564;
        v_('F1,43') = -0.015834;
        v_('F1,44') = -0.0020909;
        v_('F1,45') = 0.01344508;
        v_('F1,46') = 0.01659775;
        v_('F1,47') = -0.0014688;
        v_('F1,48') = 0.0115983;
        v_('F1,49') = 0.01107383;
        v_('F1,50') = 0.01687134;
        v_('F1,51') = 0.01485068;
        v_('F1,52') = 0.03430532;
        v_('F1,53') = -0.0043532;
        v_('F1,54') = 0.0195711;
        v_('F1,55') = 0.01179293;
        v_('F1,56') = 0.00847671;
        v_('F1,57') = 0.0146095;
        v_('F1,58') = 0.01750579;
        v_('F1,59') = 0.00833576;
        v_('F1,60') = 0.02196482;
        v_('F1,61') = -0.0148144;
        v_('F1,62') = 0.00907008;
        v_('F2,1') = 0.06839808;
        v_('F2,2') = 0.11464542;
        v_('F2,3') = 0.02391108;
        v_('F2,4') = 0.01730916;
        v_('F2,5') = 0.03985465;
        v_('F2,6') = -0.0284541;
        v_('F2,7') = 0.03014387;
        v_('F2,8') = -0.0043781;
        v_('F2,9') = 0.04352908;
        v_('F2,10') = -0.0313116;
        v_('F2,11') = -0.0604625;
        v_('F2,12') = 0.03727582;
        v_('F2,13') = 0.01361736;
        v_('F2,14') = -0.0138246;
        v_('F2,15') = 0.00565259;
        v_('F2,16') = -0.004328;
        v_('F2,17') = -0.0246133;
        v_('F2,18') = 0.02714435;
        v_('F2,19') = -0.0513326;
        v_('F2,20') = -0.0578522;
        v_('F2,21') = 0.11776573;
        v_('F2,22') = -0.0294416;
        v_('F2,23') = 0.04846583;
        v_('F2,24') = 0.0021062;
        v_('F2,25') = 0.07112832;
        v_('F2,26') = -0.0261799;
        v_('F2,27') = -0.0217402;
        v_('F2,28') = 0.01046127;
        v_('F2,29') = 0.08218002;
        v_('F2,30') = 0.06404283;
        v_('F2,31') = -0.0309792;
        v_('F2,32') = -0.03096;
        v_('F2,33') = 0.06890906;
        v_('F2,34') = -0.0099322;
        v_('F2,35') = 0.02470467;
        v_('F2,36') = 0.01070497;
        v_('F2,37') = -0.0648169;
        v_('F2,38') = -0.0867818;
        v_('F2,39') = -0.0118189;
        v_('F2,40') = -0.0292303;
        v_('F2,41') = -0.0474585;
        v_('F2,42') = 0.01299229;
        v_('F2,43') = 0.02902949;
        v_('F2,44') = 0.00934787;
        v_('F2,45') = -0.0394707;
        v_('F2,46') = 0.02680211;
        v_('F2,47') = -0.0302386;
        v_('F2,48') = -0.0227665;
        v_('F2,49') = 0.02188145;
        v_('F2,50') = -0.0726653;
        v_('F2,51') = 0.06472895;
        v_('F2,52') = 0.04080201;
        v_('F2,53') = -0.0176916;
        v_('F2,54') = -0.0174385;
        v_('F2,55') = 0.11707885;
        v_('F2,56') = -0.0196384;
        v_('F2,57') = 0.00977683;
        v_('F2,58') = 0.04120185;
        v_('F2,59') = -0.0051549;
        v_('F2,60') = 0.04551689;
        v_('F2,61') = -0.0252174;
        v_('F2,62') = 0.01001894;
        v_('F3,1') = 0.03009985;
        v_('F3,2') = 0.00113202;
        v_('F3,3') = 0.02014134;
        v_('F3,4') = 0.00609629;
        v_('F3,5') = 0.00041314;
        v_('F3,6') = -0.0022713;
        v_('F3,7') = 0.01628035;
        v_('F3,8') = 0.01568015;
        v_('F3,9') = 0.00915592;
        v_('F3,10') = -0.0006622;
        v_('F3,11') = 0.01722995;
        v_('F3,12') = 0.01700326;
        v_('F3,13') = 0.00634168;
        v_('F3,14') = -0.0091024;
        v_('F3,15') = 0.0141967;
        v_('F3,16') = 0.01805168;
        v_('F3,17') = -0.0075903;
        v_('F3,18') = -0.0169895;
        v_('F3,19') = -0.0010204;
        v_('F3,20') = 0.0119382;
        v_('F3,21') = 0.04706328;
        v_('F3,22') = -0.0065674;
        v_('F3,23') = 0.00533721;
        v_('F3,24') = 0.0316723;
        v_('F3,25') = -0.0015788;
        v_('F3,26') = 0.04234509;
        v_('F3,27') = 0.02123953;
        v_('F3,28') = 0.02167813;
        v_('F3,29') = 0.02929614;
        v_('F3,30') = -0.0015172;
        v_('F3,31') = -0.0125235;
        v_('F3,32') = 0.03889626;
        v_('F3,33') = 0.02053325;
        v_('F3,34') = -0.0012012;
        v_('F3,35') = 0.01257767;
        v_('F3,36') = 0.00841293;
        v_('F3,37') = -0.0258624;
        v_('F3,38') = -0.0102267;
        v_('F3,39') = -0.004479;
        v_('F3,40') = -0.0310854;
        v_('F3,41') = -0.0043269;
        v_('F3,42') = 0.0192379;
        v_('F3,43') = -0.030522;
        v_('F3,44') = 0.00311075;
        v_('F3,45') = 0.0101588;
        v_('F3,46') = 0.02429471;
        v_('F3,47') = -0.0176725;
        v_('F3,48') = 0.01294056;
        v_('F3,49') = 0.00420648;
        v_('F3,50') = 0.00449915;
        v_('F3,51') = 0.0;
        v_('F3,52') = 0.06800865;
        v_('F3,53') = -0.0112316;
        v_('F3,54') = 0.02515601;
        v_('F3,55') = 0.02306449;
        v_('F3,56') = 0.0073444;
        v_('F3,57') = 0.01564303;
        v_('F3,58') = 0.03084961;
        v_('F3,59') = 0.01560227;
        v_('F3,60') = 0.01436445;
        v_('F3,61') = -0.0139043;
        v_('F3,62') = 0.01262528;
        v_('F4,1') = 0.02765259;
        v_('F4,2') = 0.00825095;
        v_('F4,3') = 0.02300206;
        v_('F4,4') = -0.0109541;
        v_('F4,5') = 0.0075051;
        v_('F4,6') = -0.0013018;
        v_('F4,7') = 0.02042146;
        v_('F4,8') = 0.01036122;
        v_('F4,9') = 0.00786683;
        v_('F4,10') = -0.0013241;
        v_('F4,11') = 0.02267969;
        v_('F4,12') = 0.01487547;
        v_('F4,13') = 0.00531164;
        v_('F4,14') = -0.0015382;
        v_('F4,15') = 0.01232501;
        v_('F4,16') = 0.021637;
        v_('F4,17') = -0.0068652;
        v_('F4,18') = -0.0156514;
        v_('F4,19') = -0.0011925;
        v_('F4,20') = 0.03608384;
        v_('F4,21') = 0.03060179;
        v_('F4,22') = -0.0001242;
        v_('F4,23') = 0.01261183;
        v_('F4,24') = 0.01822198;
        v_('F4,25') = 0.0208484;
        v_('F4,26') = 0.02921733;
        v_('F4,27') = 0.01559901;
        v_('F4,28') = 0.00265402;
        v_('F4,29') = 0.03300293;
        v_('F4,30') = 0.01390252;
        v_('F4,31') = 0.02618702;
        v_('F4,32') = 0.02305596;
        v_('F4,33') = 0.00676091;
        v_('F4,34') = -0.0050875;
        v_('F4,35') = 0.02495398;
        v_('F4,36') = -0.0082318;
        v_('F4,37') = -0.0199205;
        v_('F4,38') = -0.0122157;
        v_('F4,39') = -0.0138737;
        v_('F4,40') = -0.0156497;
        v_('F4,41') = -0.0069589;
        v_('F4,42') = 0.01854347;
        v_('F4,43') = -0.0196348;
        v_('F4,44') = -0.0026452;
        v_('F4,45') = -0.0041136;
        v_('F4,46') = 0.02315343;
        v_('F4,47') = -0.0198672;
        v_('F4,48') = 0.01392878;
        v_('F4,49') = -0.0044366;
        v_('F4,50') = 0.00977181;
        v_('F4,51') = 0.00212687;
        v_('F4,52') = 0.0311986;
        v_('F4,53') = 0.00308721;
        v_('F4,54') = 0.03221339;
        v_('F4,55') = 0.01610098;
        v_('F4,56') = -0.0060155;
        v_('F4,57') = 0.01156268;
        v_('F4,58') = 0.03117856;
        v_('F4,59') = 0.02853774;
        v_('F4,60') = 0.01889475;
        v_('F4,61') = -0.0084619;
        v_('F4,62') = 0.01389078;
        v_('F5,1') = 0.01934344;
        v_('F5,2') = -0.0001824;
        v_('F5,3') = 0.01368738;
        v_('F5,4') = -0.0011702;
        v_('F5,5') = 0.0063987;
        v_('F5,6') = -0.0012536;
        v_('F5,7') = 0.01631848;
        v_('F5,8') = 0.00855757;
        v_('F5,9') = 0.00629811;
        v_('F5,10') = 0.00756259;
        v_('F5,11') = 0.01708222;
        v_('F5,12') = 0.01891594;
        v_('F5,13') = 0.00524476;
        v_('F5,14') = -0.0021532;
        v_('F5,15') = 0.00431571;
        v_('F5,16') = 0.0098339;
        v_('F5,17') = 0.00196399;
        v_('F5,18') = 0.00236851;
        v_('F5,19') = 0.01417746;
        v_('F5,20') = 0.02771752;
        v_('F5,21') = 0.01993433;
        v_('F5,22') = 0.00774124;
        v_('F5,23') = 0.01201704;
        v_('F5,24') = 0.01525628;
        v_('F5,25') = 0.02701902;
        v_('F5,26') = 0.00641488;
        v_('F5,27') = -0.0007161;
        v_('F5,28') = 0.00121838;
        v_('F5,29') = 0.01281317;
        v_('F5,30') = 0.03569157;
        v_('F5,31') = 0.03316501;
        v_('F5,32') = 0.02476882;
        v_('F5,33') = 0.01933613;
        v_('F5,34') = -0.0151122;
        v_('F5,35') = -0.0032742;
        v_('F5,36') = -0.0117874;
        v_('F5,37') = -0.0197497;
        v_('F5,38') = 0.00079793;
        v_('F5,39') = -0.0090359;
        v_('F5,40') = -0.0077774;
        v_('F5,41') = 0.0058112;
        v_('F5,42') = 0.0069869;
        v_('F5,43') = -0.0068717;
        v_('F5,44') = -0.0115544;
        v_('F5,45') = 0.00686421;
        v_('F5,46') = 0.0161998;
        v_('F5,47') = -0.009432;
        v_('F5,48') = 0.02105545;
        v_('F5,49') = 0.00866881;
        v_('F5,50') = 0.0170584;
        v_('F5,51') = 0.01658024;
        v_('F5,52') = 0.02052897;
        v_('F5,53') = -0.0102431;
        v_('F5,54') = 0.01558603;
        v_('F5,55') = 0.01086556;
        v_('F5,56') = 0.01153823;
        v_('F5,57') = 0.01518881;
        v_('F5,58') = 0.02004731;
        v_('F5,59') = 0.01107311;
        v_('F5,60') = 0.01026376;
        v_('F5,61') = -0.0163459;
        v_('F5,62') = 0.00490451;
        v_('F6,1') = 0.02026545;
        v_('F6,2') = 0.02412925;
        v_('F6,3') = 0.02670218;
        v_('F6,4') = 0.00685114;
        v_('F6,5') = 0.00436018;
        v_('F6,6') = -0.0003288;
        v_('F6,7') = 0.01644953;
        v_('F6,8') = 0.01831952;
        v_('F6,9') = 0.01716356;
        v_('F6,10') = -0.0061246;
        v_('F6,11') = 0.01024964;
        v_('F6,12') = 0.02626665;
        v_('F6,13') = 0.00958273;
        v_('F6,14') = -0.001622;
        v_('F6,15') = 0.0159456;
        v_('F6,16') = 0.00527126;
        v_('F6,17') = -0.0166735;
        v_('F6,18') = -0.025704;
        v_('F6,19') = 0.00830207;
        v_('F6,20') = -0.1286899;
        v_('F6,21') = 0.06852863;
        v_('F6,22') = -0.0098919;
        v_('F6,23') = -0.0326849;
        v_('F6,24') = 0.01682627;
        v_('F6,25') = -0.0275797;
        v_('F6,26') = -0.0260791;
        v_('F6,27') = 0.08026138;
        v_('F6,28') = 0.04569663;
        v_('F6,29') = 0.06344316;
        v_('F6,30') = 0.02979956;
        v_('F6,31') = 0.04064994;
        v_('F6,32') = 0.01009655;
        v_('F6,33') = 0.01584007;
        v_('F6,34') = -0.0668351;
        v_('F6,35') = 0.0445981;
        v_('F6,36') = 0.00965304;
        v_('F6,37') = -0.0387893;
        v_('F6,38') = 0.04251449;
        v_('F6,39') = 0.00283502;
        v_('F6,40') = -0.0290312;
        v_('F6,41') = -0.0475923;
        v_('F6,42') = -0.0070546;
        v_('F6,43') = -0.0179988;
        v_('F6,44') = 0.01344507;
        v_('F6,45') = -0.0137426;
        v_('F6,46') = 0.00554952;
        v_('F6,47') = -0.0219556;
        v_('F6,48') = 0.00705348;
        v_('F6,49') = -0.0832572;
        v_('F6,50') = -0.0902205;
        v_('F6,51') = 0.04980283;
        v_('F6,52') = 0.07888147;
        v_('F6,53') = -0.0188265;
        v_('F6,54') = 0.05618347;
        v_('F6,55') = 0.06146955;
        v_('F6,56') = -0.0185803;
        v_('F6,57') = 0.01385571;
        v_('F6,58') = 0.04288407;
        v_('F6,59') = 0.0279033;
        v_('F6,60') = 0.07149137;
        v_('F6,61') = 0.00892353;
        v_('F6,62') = -0.0082346;
        v_('F7,1') = 0.02575342;
        v_('F7,2') = 0.04711538;
        v_('F7,3') = 0.04213856;
        v_('F7,4') = 0.00156648;
        v_('F7,5') = 0.04701857;
        v_('F7,6') = -0.0128839;
        v_('F7,7') = 0.02269933;
        v_('F7,8') = 0.00268196;
        v_('F7,9') = 0.02979155;
        v_('F7,10') = -0.0151365;
        v_('F7,11') = 0.00081848;
        v_('F7,12') = 0.06069968;
        v_('F7,13') = -0.0102801;
        v_('F7,14') = -0.0141954;
        v_('F7,15') = -0.0088682;
        v_('F7,16') = 0.02879164;
        v_('F7,17') = -0.0216998;
        v_('F7,18') = -0.0183082;
        v_('F7,19') = -0.0087868;
        v_('F7,20') = 0.04911805;
        v_('F7,21') = 0.06639076;
        v_('F7,22') = 0.03177555;
        v_('F7,23') = 0.02006112;
        v_('F7,24') = 0.00975647;
        v_('F7,25') = 0.10423007;
        v_('F7,26') = -0.010197;
        v_('F7,27') = 0.01517472;
        v_('F7,28') = 0.040661;
        v_('F7,29') = 0.09204718;
        v_('F7,30') = 0.075359;
        v_('F7,31') = -0.0047691;
        v_('F7,32') = -0.0071033;
        v_('F7,33') = 0.02237111;
        v_('F7,34') = 0.01449517;
        v_('F7,35') = -0.0150544;
        v_('F7,36') = -0.0181191;
        v_('F7,37') = 0.01805729;
        v_('F7,38') = -0.0015568;
        v_('F7,39') = 0.0053461;
        v_('F7,40') = -0.0277516;
        v_('F7,41') = 0.00586828;
        v_('F7,42') = -0.0188615;
        v_('F7,43') = -0.0192241;
        v_('F7,44') = 0.009359;
        v_('F7,45') = -0.0105551;
        v_('F7,46') = 0.02764189;
        v_('F7,47') = -0.0143381;
        v_('F7,48') = -0.0082625;
        v_('F7,49') = 0.00481108;
        v_('F7,50') = 0.08875394;
        v_('F7,51') = 0.05749222;
        v_('F7,52') = 0.05203367;
        v_('F7,53') = -0.0166313;
        v_('F7,54') = -0.0419138;
        v_('F7,55') = -0.0565902;
        v_('F7,56') = -0.0079726;
        v_('F7,57') = -0.0460882;
        v_('F7,58') = 0.04327144;
        v_('F7,59') = -0.0362028;
        v_('F7,60') = -0.0009689;
        v_('F7,61') = -0.0056484;
        v_('F7,62') = 0.00556576;
        v_('F8,1') = 0.01944241;
        v_('F8,2') = -0.0001799;
        v_('F8,3') = 0.01340651;
        v_('F8,4') = -0.0015981;
        v_('F8,5') = 0.0084482;
        v_('F8,6') = -0.0003527;
        v_('F8,7') = 0.01270289;
        v_('F8,8') = 0.00749129;
        v_('F8,9') = 0.00717621;
        v_('F8,10') = 0.00429221;
        v_('F8,11') = 0.01931789;
        v_('F8,12') = 0.0208805;
        v_('F8,13') = 0.0122392;
        v_('F8,14') = -0.0060861;
        v_('F8,15') = 0.01232854;
        v_('F8,16') = 0.00701669;
        v_('F8,17') = 0.00448502;
        v_('F8,18') = 0.00039866;
        v_('F8,19') = 0.01235355;
        v_('F8,20') = 0.0371595;
        v_('F8,21') = 0.0224685;
        v_('F8,22') = 0.00853749;
        v_('F8,23') = 0.02164152;
        v_('F8,24') = 0.01433821;
        v_('F8,25') = 0.03395369;
        v_('F8,26') = 0.0005496;
        v_('F8,27') = 0.00205987;
        v_('F8,28') = 0.00356311;
        v_('F8,29') = 0.02239519;
        v_('F8,30') = 0.02758114;
        v_('F8,31') = 0.04185351;
        v_('F8,32') = 0.02545069;
        v_('F8,33') = 0.0194659;
        v_('F8,34') = -0.014977;
        v_('F8,35') = 0.00411921;
        v_('F8,36') = -0.0118846;
        v_('F8,37') = -0.0260699;
        v_('F8,38') = -0.0097793;
        v_('F8,39') = -0.0127247;
        v_('F8,40') = -0.0138505;
        v_('F8,41') = 0.00429157;
        v_('F8,42') = 0.01786986;
        v_('F8,43') = -0.0171108;
        v_('F8,44') = -0.0081542;
        v_('F8,45') = 0.00593762;
        v_('F8,46') = 0.01381592;
        v_('F8,47') = -0.0106845;
        v_('F8,48') = 0.01985384;
        v_('F8,49') = 0.01020926;
        v_('F8,50') = 0.01952169;
        v_('F8,51') = 0.02056397;
        v_('F8,52') = 0.02473456;
        v_('F8,53') = -0.0103614;
        v_('F8,54') = 0.01671624;
        v_('F8,55') = 0.01351589;
        v_('F8,56') = 0.00889043;
        v_('F8,57') = 0.01487755;
        v_('F8,58') = 0.02683807;
        v_('F8,59') = 0.01081704;
        v_('F8,60') = 0.01298278;
        v_('F8,61') = -0.022308;
        v_('F8,62') = 0.00526547;
        v_('F9,1') = 0.03257355;
        v_('F9,2') = 0.11882583;
        v_('F9,3') = 0.0200098;
        v_('F9,4') = 0.00994581;
        v_('F9,5') = 0.04788101;
        v_('F9,6') = -0.023203;
        v_('F9,7') = 0.02328976;
        v_('F9,8') = -0.0285306;
        v_('F9,9') = 0.01515152;
        v_('F9,10') = -0.0187389;
        v_('F9,11') = -0.0351782;
        v_('F9,12') = 0.04382249;
        v_('F9,13') = 0.01976048;
        v_('F9,14') = 0.0000652;
        v_('F9,15') = 0.00978601;
        v_('F9,16') = -0.0084636;
        v_('F9,17') = -0.0374014;
        v_('F9,18') = -0.0066337;
        v_('F9,19') = -0.0414991;
        v_('F9,20') = 0.02509598;
        v_('F9,21') = 0.06727235;
        v_('F9,22') = 0.03385535;
        v_('F9,23') = 0.0318039;
        v_('F9,24') = 0.01803119;
        v_('F9,25') = 0.04092867;
        v_('F9,26') = -0.0186249;
        v_('F9,27') = -0.0076148;
        v_('F9,28') = -0.0016527;
        v_('F9,29') = 0.09317725;
        v_('F9,30') = 0.04970254;
        v_('F9,31') = -0.0021639;
        v_('F9,32') = -0.0077451;
        v_('F9,33') = 0.03382422;
        v_('F9,34') = -0.0108723;
        v_('F9,35') = -0.0003053;
        v_('F9,36') = 0.00463222;
        v_('F9,37') = -0.0396736;
        v_('F9,38') = -0.0400464;
        v_('F9,39') = -0.016379;
        v_('F9,40') = -0.0098904;
        v_('F9,41') = -0.0321688;
        v_('F9,42') = 0.01481136;
        v_('F9,43') = -0.005114;
        v_('F9,44') = -0.0292249;
        v_('F9,45') = -0.0286768;
        v_('F9,46') = 0.03834375;
        v_('F9,47') = -0.0092614;
        v_('F9,48') = -0.0008335;
        v_('F9,49') = -0.014421;
        v_('F9,50') = -0.0623375;
        v_('F9,51') = 0.02637348;
        v_('F9,52') = 0.05629201;
        v_('F9,53') = -0.0117171;
        v_('F9,54') = -0.0052961;
        v_('F9,55') = 0.07478219;
        v_('F9,56') = -0.0229678;
        v_('F9,57') = 0.00530076;
        v_('F9,58') = 0.03977533;
        v_('F9,59') = 0.00440966;
        v_('F9,60') = 0.04653715;
        v_('F9,61') = -0.0302045;
        v_('F9,62') = -0.0035146;
        v_('F10,1') = 0.02936247;
        v_('F10,2') = 0.02728753;
        v_('F10,3') = 0.01895524;
        v_('F10,4') = 0.02426429;
        v_('F10,5') = -0.0079572;
        v_('F10,6') = -0.0007959;
        v_('F10,7') = 0.03278387;
        v_('F10,8') = -0.0042126;
        v_('F10,9') = 0.01608771;
        v_('F10,10') = -0.0134287;
        v_('F10,11') = 0.01224441;
        v_('F10,12') = 0.0277158;
        v_('F10,13') = 0.01154154;
        v_('F10,14') = -0.0000565;
        v_('F10,15') = 0.01739818;
        v_('F10,16') = 0.01421354;
        v_('F10,17') = -0.0241966;
        v_('F10,18') = -0.0285554;
        v_('F10,19') = -0.0282975;
        v_('F10,20') = -0.083264;
        v_('F10,21') = 0.03474878;
        v_('F10,22') = -0.0091473;
        v_('F10,23') = 0.02921277;
        v_('F10,24') = 0.03329852;
        v_('F10,25') = 0.00214044;
        v_('F10,26') = 0.02699496;
        v_('F10,27') = -0.021606;
        v_('F10,28') = -0.0312352;
        v_('F10,29') = 0.07112818;
        v_('F10,30') = -0.0417093;
        v_('F10,31') = 0.09619381;
        v_('F10,32') = 0.00996696;
        v_('F10,33') = 0.04349692;
        v_('F10,34') = -0.0438425;
        v_('F10,35') = 0.00419287;
        v_('F10,36') = 0.02622986;
        v_('F10,37') = -0.0362526;
        v_('F10,38') = -0.003085;
        v_('F10,39') = -0.0045062;
        v_('F10,40') = -0.0205606;
        v_('F10,41') = -0.0288991;
        v_('F10,42') = 0.02540138;
        v_('F10,43') = -0.0243807;
        v_('F10,44') = 0.0041841;
        v_('F10,45') = 0.00393836;
        v_('F10,46') = 0.01091591;
        v_('F10,47') = -0.0323941;
        v_('F10,48') = -0.0103458;
        v_('F10,49') = 0.00099841;
        v_('F10,50') = -0.0562074;
        v_('F10,51') = 0.0586224;
        v_('F10,52') = 0.03905103;
        v_('F10,53') = -0.0012433;
        v_('F10,54') = 0.04589181;
        v_('F10,55') = 0.02250717;
        v_('F10,56') = -0.0007407;
        v_('F10,57') = 0.0114906;
        v_('F10,58') = 0.0482672;
        v_('F10,59') = 0.02911506;
        v_('F10,60') = 0.02470034;
        v_('F10,61') = 0.00525668;
        v_('F10,62') = 0.00579451;
        v_('F11,1') = 0.0142632;
        v_('F11,2') = 0.02913527;
        v_('F11,3') = 0.00045297;
        v_('F11,4') = -0.0110172;
        v_('F11,5') = 0.00610407;
        v_('F11,6') = 0.02745336;
        v_('F11,7') = 0.01660762;
        v_('F11,8') = 0.01974878;
        v_('F11,9') = 0.0011392;
        v_('F11,10') = -0.0152905;
        v_('F11,11') = 0.00137224;
        v_('F11,12') = 0.04046159;
        v_('F11,13') = 0.01122972;
        v_('F11,14') = -0.0300932;
        v_('F11,15') = 0.07110043;
        v_('F11,16') = 0.02401848;
        v_('F11,17') = -0.0161092;
        v_('F11,18') = -0.0216779;
        v_('F11,19') = -0.0265765;
        v_('F11,20') = -0.0605873;
        v_('F11,21') = 0.00856515;
        v_('F11,22') = -0.0087101;
        v_('F11,23') = 0.03675771;
        v_('F11,24') = -0.0057913;
        v_('F11,25') = 0.00028415;
        v_('F11,26') = 0.04552234;
        v_('F11,27') = 0.0116832;
        v_('F11,28') = -0.002417;
        v_('F11,29') = 0.06218872;
        v_('F11,30') = 0.0644405;
        v_('F11,31') = 0.01547711;
        v_('F11,32') = -0.0056861;
        v_('F11,33') = 0.03572692;
        v_('F11,34') = 0.01895492;
        v_('F11,35') = 0.02731691;
        v_('F11,36') = 0.00391517;
        v_('F11,37') = -0.0623984;
        v_('F11,38') = -0.0530907;
        v_('F11,39') = 0.00024404;
        v_('F11,40') = -0.0496492;
        v_('F11,41') = -0.0034015;
        v_('F11,42') = 0.00772798;
        v_('F11,43') = 0.00134202;
        v_('F11,44') = -0.0003191;
        v_('F11,45') = 0.01589632;
        v_('F11,46') = 0.02193175;
        v_('F11,47') = -0.0211536;
        v_('F11,48') = 0.00986305;
        v_('F11,49') = -0.0356454;
        v_('F11,50') = -0.0262546;
        v_('F11,51') = 0.01305068;
        v_('F11,52') = 0.04276746;
        v_('F11,53') = -0.04001;
        v_('F11,54') = 0.02436634;
        v_('F11,55') = 0.04495887;
        v_('F11,56') = -0.0062248;
        v_('F11,57') = 0.00411447;
        v_('F11,58') = 0.02287322;
        v_('F11,59') = 0.0187145;
        v_('F11,60') = 0.01414485;
        v_('F11,61') = -0.0208345;
        v_('F11,62') = -0.0014776;
        v_('F12,1') = 0.04961332;
        v_('F12,2') = 0.10977209;
        v_('F12,3') = 0.04399932;
        v_('F12,4') = -0.0011856;
        v_('F12,5') = 0.04759078;
        v_('F12,6') = -0.0030388;
        v_('F12,7') = 0.02769167;
        v_('F12,8') = -0.0027146;
        v_('F12,9') = 0.02207884;
        v_('F12,10') = -0.0142039;
        v_('F12,11') = -0.0637882;
        v_('F12,12') = 0.01581788;
        v_('F12,13') = 0.0276185;
        v_('F12,14') = 0.03916249;
        v_('F12,15') = 0.01596138;
        v_('F12,16') = -0.006546;
        v_('F12,17') = -0.0384127;
        v_('F12,18') = 0.00045683;
        v_('F12,19') = -0.121106;
        v_('F12,20') = 0.01876118;
        v_('F12,21') = 0.06550317;
        v_('F12,22') = 0.01441183;
        v_('F12,23') = 0.04188729;
        v_('F12,24') = 0.00322029;
        v_('F12,25') = 0.08396028;
        v_('F12,26') = 0.00564501;
        v_('F12,27') = 0.00358885;
        v_('F12,28') = -0.041766;
        v_('F12,29') = 0.07382422;
        v_('F12,30') = 0.10782392;
        v_('F12,31') = -0.0317326;
        v_('F12,32') = -0.055244;
        v_('F12,33') = 0.08599692;
        v_('F12,34') = -0.0217804;
        v_('F12,35') = 0.03927492;
        v_('F12,36') = 0.05519274;
        v_('F12,37') = -0.0400785;
        v_('F12,38') = -0.0666771;
        v_('F12,39') = -0.005139;
        v_('F12,40') = 0.0202388;
        v_('F12,41') = -0.0702606;
        v_('F12,42') = 0.02115788;
        v_('F12,43') = 0.01612974;
        v_('F12,44') = -0.0542889;
        v_('F12,45') = -0.0285207;
        v_('F12,46') = 0.08278316;
        v_('F12,47') = 0.01565405;
        v_('F12,48') = -0.0458125;
        v_('F12,49') = -0.0263263;
        v_('F12,50') = -0.0636084;
        v_('F12,51') = 0.02148485;
        v_('F12,52') = 0.05500192;
        v_('F12,53') = -0.0387375;
        v_('F12,54') = 0.03604668;
        v_('F12,55') = 0.1128135;
        v_('F12,56') = -0.0066382;
        v_('F12,57') = -0.0038363;
        v_('F12,58') = 0.03565365;
        v_('F12,59') = -0.0019992;
        v_('F12,60') = 0.05372596;
        v_('F12,61') = -0.0050188;
        v_('F12,62') = 0.02036761;
        v_('R1') = 0.00369588;
        v_('R2') = 0.03930262;
        v_('R3') = 0.03651482;
        v_('R4') = 0.01053366;
        v_('R5') = 0.03285931;
        v_('R6') = 0.00728512;
        v_('R7') = -0.00696702;
        v_('R8') = 0.01035681;
        v_('R9') = 0.02017062;
        v_('R10') = -0.00907559;
        v_('R11') = -0.01570064;
        v_('R12') = -0.01934069;
        v_('R13') = 0.03483565;
        v_('R14') = 0.01598009;
        v_('R15') = 0.01650229;
        v_('R16') = 0.00145856;
        v_('R17') = 0.00512918;
        v_('R18') = -0.02450702;
        v_('R19') = -0.00238956;
        v_('R20') = -0.01948598;
        v_('R21') = -0.03710551;
        v_('R22') = 0.04642074;
        v_('R23') = -0.01552978;
        v_('R24') = 0.00938498;
        v_('R25') = 0.03560831;
        v_('R26') = 0.03686724;
        v_('R27') = -0.00970278;
        v_('R28') = 0.01240233;
        v_('R29') = 0.00667647;
        v_('R30') = 0.05287496;
        v_('R31') = 0.03461627;
        v_('R32') = -0.02759314;
        v_('R33') = 0.03078867;
        v_('R34') = 0.01025355;
        v_('R35') = 0.0054057;
        v_('R36') = 0.04416525;
        v_('R37') = -0.01082388;
        v_('R38') = -0.04515032;
        v_('R39') = 0.00072319;
        v_('R40') = -0.04013564;
        v_('R41') = -0.01279898;
        v_('R42') = -0.04006805;
        v_('R43') = 0.01625619;
        v_('R44') = -0.02273137;
        v_('R45') = 0.01550674;
        v_('R46') = 0.00830152;
        v_('R47') = 0.00955529;
        v_('R48') = -0.02291803;
        v_('R49') = 0.00371634;
        v_('R50') = -0.04345979;
        v_('R51') = 0.00475919;
        v_('R52') = 0.04919793;
        v_('R53') = -0.01516884;
        v_('R54') = 0.01552472;
        v_('R55') = 0.03075534;
        v_('R56') = 0.05074156;
        v_('R57') = -0.02622951;
        v_('R58') = 0.0047937;
        v_('R59') = 0.02527404;
        v_('R60') = 0.02016397;
        v_('R61') = 0.02068853;
        v_('R62') = -0.007714;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('NS')
            [iv,ix_] = s2mpjlib('ii',['S',int2str(I)],ix_);
            pb.xnames{iv} = ['S',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('NR')
            for J=v_('1'):v_('NS')
                [ig,ig_] = s2mpjlib('ii',['A',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['S',int2str(J)]);
                valA(end+1) = v_(['F',int2str(J),',',int2str(I)]);
            end
        end
        for I=v_('1'):v_('NS')
            [ig,ig_] = s2mpjlib('ii','SUM',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'SUM';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = 1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        legrps = find(strcmp(gtype,'<='));
        eqgrps = find(strcmp(gtype,'=='));
        gegrps = find(strcmp(gtype,'>='));
        pb.nle = length(legrps);
        pb.neq = length(eqgrps);
        pb.nge = length(gegrps);
        pb.m   = pb.nle+pb.neq+pb.nge;
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('NR')
            pbm.gconst(ig_(['A',int2str(I)])) = v_(['R',int2str(I)]);
        end
        pbm.gconst(ig_('SUM')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.0*ones(pb.n,1);
        pb.xupper = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('RNS') = v_('NS');
        v_('SINI') = 1.0/v_('RNS');
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = v_('SINI')*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('NR')
            ig = ig_(['A',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN                3.27497293D-2
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CSLR2-MN-12-1';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

