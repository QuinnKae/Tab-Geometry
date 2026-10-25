clear;

P  = askNum('Enter axial in-plane force P (lbf): ');
V  = askNum('Enter transverse in-plane force V (lbf): ');
Fz = askNum('Enter force perpendicular to tab plane Fz (lbf): ');
 
D  = askNum('Enter mounting HOLE diameter (in): ');
Db = askNum('Enter nominal BOLT diameter (in): ');
OD = askNum('Enter chassis tube outer diameter (in): ');
t  = askNum('Enter tab thickness (in): ');
A  = askNum(['Enter tab angle A (0 to 90 deg): ']);
 
SyTab = askNum('Enter tab material yield strength (psi): ');
FuTab = askNum('Enter tab material ultimate strength (psi): ');
Fexx  = askNum('Enter weld filler tensile strength Fexx (psi): ');
Wf    = askNum('Enter fillet weld LEG size (in): ');
N     = askNum('Enter tab material factor of safety (>=1): ');
 
Hmin = askNum('Desired minimum hole-edge to tube clearance (in) [0.25]: ', 0.25);
  
if any([D,Db,OD,t,SyTab,FuTab,Fexx,Wf] <= 0) || Db > D
    error('Dimensions/strengths must be positive and bolt diameter <= hole diameter.');
end
 
if N < 1 || Hmin < 0 || A < 0 || A > 90
    error('Require N >= 1, Hmin >= 0, and 0 <= A <= 90 degrees.');
end
 
if FuTab < SyTab
    error('Entered ultimate strength must be >= entered yield strength.');
end
 
if P == 0 && V == 0 && Fz == 0
    error('Zero load: no positive stress-controlled minimum weld length exists.');
end
  
E = 29000000;  % Young's modulus (psi)
 
SyContinuous = min(SyTab, 44200);
FuContinuous = min(FuTab, 52900);
 
SyJoint = min(SyTab, 26000);
FuJoint = min(FuTab, 43500);
 
Fr = hypot(P,V);
 
TabAllow = SyJoint/N;
 
TauTabAllow = TabAllow/sqrt(3);
 
Q = Fr/(2*t*TauTabAllow);
 
W = Fr/(2*t*TabAllow);
 
H = max([Q,Hmin]);
  
OmegaHole = 2.0;
Fb = Fr;
 
RaBearing = 2.4*Db*t*FuJoint/OmegaHole;
BearingUtil = Fb/RaBearing;
tBearingRequired = OmegaHole*Fb/(2.4*Db*FuJoint);
BearingPass = BearingUtil <= 1;
 
LcTrial = H;
RaTearTrial = 1.2*LcTrial*t*FuJoint/OmegaHole;
 
if Fb == 0
    TearUtilTrial = 0;
elseif RaTearTrial == 0
    TearUtilTrial = Inf;
else
    TearUtilTrial = Fb/RaTearTrial;
end
 
TearTrialPass = TearUtilTrial <= 1;
 
LcRequired = OmegaHole*Fb/(1.2*t*FuJoint);
LcDesign = max([LcTrial,LcRequired]);
HoleCenterToFreeEdge = D/2 + LcDesign;
 
WeldAllow = 0.30*Fexx;   
te = Wf/sqrt(2);         
 
r = OD/2;   
 
if A == 0
    a = Inf;
else
    a = r/sind(A);
end
 
weldDemand = @(L) curvedWeldDemand(L,P,V,Fz,r,a,te,H,D,t);
  
Llow = 1e-10*OD;
nonMonotonicWarning = false;
fAtLmax = NaN;
 
if isinf(a)
 
    Lhigh = max([OD,D,H,Wf,0.1]);
 
    while weldDemand(Lhigh) > WeldAllow
        Lhigh = 2*Lhigh;
        if Lhigh > 1e6*max(OD,1)
            error('Could not bracket a passing straight-weld length.');
        end
    end
 
    if weldDemand(Llow) <= WeldAllow
        error('Lower search bound already passes: reduce Llow.');
    end
 
    L_solution = fzero(@(L) weldDemand(L)-WeldAllow,[Llow,Lhigh]);
 
else

    Lmax = 2*a;
 
    Lgrid = logspace(log10(Llow),log10(Lmax),1500);
    Lgrid(end) = Lmax;   
 
    stressError = arrayfun(weldDemand,Lgrid)-WeldAllow;
    fAtLmax = stressError(end)+WeldAllow;
 
    idx = find(stressError(1:end-1)>0 & stressError(2:end)<=0,1,'first');
 
    if isempty(idx)
        error(['No passing weld span on this ellipse branch. ', ...
               'Change design inputs.']);
    end
 
    if any(stressError(idx+1:end) > 0)
        nonMonotonicWarning = true;
    end
 
    L_solution = fzero(@(L) weldDemand(L)-WeldAllow,[Lgrid(idx),Lgrid(idx+1)]);
 
end
 
[fmax,S,e,Mx,Mz] = curvedWeldDemand(L_solution,P,V,Fz,r,a,te,H,D,t);
 
WeldUtilization = fmax/WeldAllow;
WeldNominalRatio = 0.60*Fexx/fmax;
 
TabWidth = max(L_solution, D + 2*W);
 
fprintf('\n==========================================\n');
fprintf(' FSAE PRELIMINARY TAB / WELD SIZING\n');
fprintf('==========================================\n');
 
fprintf('\nAPPLIED LOADS\n');
fprintf('Axial in-plane force P: %.2f lbf\n',P);
fprintf('Transverse in-plane force V: %.2f lbf\n',V);
fprintf('Out-of-plane force Fz: %.2f lbf\n',Fz);
fprintf('In-plane resultant Fr: %.2f lbf\n',Fr);
 
fprintf('\nMATERIAL PROPERTIES USED\n');
fprintf('Continuous steel yield: %.0f psi\n',SyContinuous);
fprintf('Continuous steel ultimate: %.0f psi\n',FuContinuous);
fprintf('Welded joint yield: %.0f psi\n',SyJoint);
fprintf('Welded joint ultimate: %.0f psi\n',FuJoint);
fprintf('Weld filler Fexx: %.0f psi\n',Fexx);
 
fprintf('\nPRELIMINARY TAB GEOMETRY\n');
fprintf('Shear tear-out ligament Q: %.5f in\n',Q);
fprintf('Net-section side ligament W: %.5f in\n',W);
fprintf('Minimum tab width at weld: %.5f in\n',TabWidth);
 
if BearingPass
    fprintf('Bearing: PASS\n');
else
    fprintf('Bearing: FAIL - change tab thickness/bolt/material\n');
end
  
if TearTrialPass
    fprintf('Provisional tear-out check: PASS\n');
else
    fprintf('Provisional tear-out check: FAIL\n');
end
  
fprintf('\nSTRESS-CONTROLLED WELD LENGTH\n');
fprintf('Minimum projected span: %.6f in\n',L_solution);
 
if ~isnan(fAtLmax)
    fprintf('Weld demand at full ellipse span (Lmax): %.2f psi\n',fAtLmax);
end
 
if nonMonotonicWarning
    fprintf(['WARNING: weld stress exceeds allowable again at some longer ', ...
             'span.\n         Do not lengthen the weld beyond the ', ...
             'reported span without rechecking.\n']);
end
 
input('Press Enter to finish...','s');
 
function x = askNum(msg,defaultVal)
   
    while true
        try
            x = input(msg);
        catch
            x = [];
            fprintf('Invalid entry. Enter a single number.\n');
            continue
        end
        if isempty(x) && nargin >= 2
            x = defaultVal;
            return
        end
        if isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x)
            return
        end
        fprintf('Enter a single finite real number.\n');
    end
end
 
function [fmax,S,e,Mx,Mz] = curvedWeldDemand(L,P,V,Fz,r,a,te,H,D,t)
 
    n = 1201;  
 
    if ~isscalar(L) || ~isfinite(L) || L<=0 || ...
            (~isinf(a) && L>2*a*(1+1e-12))
        error('Invalid projected weld span L.');
    end
 
    if isinf(a)
  
        q = linspace(-L/2,L/2,n);
        x = q;
        y = r*ones(size(q));
        dsdq = ones(size(q));
    else
        
        qmax = asin(min(1,L/(2*a)));
        q = linspace(-qmax,qmax,n);
        x = a*sin(q);
        y = r*cos(q);
        dsdq = hypot(a*cos(q),r*sin(q));
    end
 
    
 
    S = trapz(q,dsdq);     
    Aw = 2*te*S;           
 
    ybar = trapz(q,y.*dsdq)/S;
    zw = t/2;
 
    Jzz = 2*te*trapz(q,(x.^2+(y-ybar).^2).*dsdq);
    Ixx = 2*te*trapz(q,((y-ybar).^2+zw^2).*dsdq);
 
    yHole = r+H+D/2;
    e = yHole-ybar;
 
    Mz = -V*e;
    Mx = Fz*e;
 
    tauXDirect = V/Aw;
    tauYDirect = P/Aw;
    tauZDirect = Fz/Aw;
 
    tauXMoment = -Mz*(y-ybar)/Jzz;
    tauYMoment = Mz*x/Jzz;
    tauZMoment = Mx*(y-ybar)/Ixx;
 
    tauYMomentPlus = -Mx*zw/Ixx;
    tauYMomentMinus = Mx*zw/Ixx;
 
    taux = tauXDirect+tauXMoment;
    tauyPlus = tauYDirect+tauYMoment+tauYMomentPlus;
    tauyMinus = tauYDirect+tauYMoment+tauYMomentMinus;
    tauz = tauZDirect+tauZMoment;
 
    fPlus = hypot(hypot(taux,tauyPlus),tauz);
    fMinus = hypot(hypot(taux,tauyMinus),tauz);
 
    fmax = max([fPlus,fMinus]);
 
end