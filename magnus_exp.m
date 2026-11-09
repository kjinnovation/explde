function Omega = magnus_exp(fA,t0,t1,ord)
%MAGNUS_EXP Magnus exponent computation
%   This function calculates the Magnus exponent, as described in the
%   following publications (cited herein as "Doc0", "Doc1" and "Doc2"):
%   * Doc0 (tutorial introduction): K. C. Johnson, "Quadrature
%     approximation to the Magnus exponent for numerical solution of
%     non-autonomous coupled linear differential equations," Version
%     09-Nov-2026 (explde.pdf, https://vixra.org/abs/2604.0065)
%   * Doc1: K. C. Johnson, "The centered Magnus exponent and Magnus
%     quadrature formulas of orders 2, 4, 6 and 8," Version 09-Nov-2026
%     (magnus_order_2-8.pdf, https://ai.vixra.org/abs/2611.0005)
%   * Doc2: K. C. Johnson, "Magnus quadrature formulas of order 10 on a
%     centered interval," Version 09-Nov-2026 (magnus_order_10.pdf,
%     https://ai.vixra.org/abs/2611.0004)
%   * Code and documentation repositories:
%     https://github.com/kjinnovation/explde
%     https://www.mathworks.com/matlabcentral/fileexchange/183695-explde
%
% The Magnus function Omega(t) is the solution of the differential equation
%   (d/dt)exp(Omega(t)) = fA(t)*exp(Omega(t)), with Omega(t0) = zeros(...)
% The function call magnus_exp(fA,t0,t1,ord) returns an approximation to
% Omega(t1), based on a quadrature approximation (using Gauss-Legendre
% quadrature).
%
% Syntax:
%   Omega = magnus_exp(fA,t0,t1,ord);
%
% Inputs:
%
%   fA: function handle. fA accepts a scalar argument t (or more generally,
%   a size-[1,1,...] array t for paged operation) and returns an array A =
%   fA(t) with size(A,1)==size(A,2), and size(A,n)==size(t,n) for n>2.
%
%   t0, t1: integration limits, scalar (or size-[1,1,...] and
%   size-compatible for paged operation).
%
%   ord: quadrature approximation order, either 2, 4, 6, 8, or 10.
%
% Output:
%
%   Omega: Magnus exponent at t = t1, size-matched to A.

% version 09-Nov-2026, Ken Johnson, ken@kjinnovation.com

persistent x
% Argument validation
narginchk(4,4)
if ~isa(fA,'function_handle')
    error('explde:Validation', ...
        ['In magnus_exp(fA,t0,t1,ord): ' ...
        'Invalid fA (must be a function handle).']);
end
if ~(isscalar(t0) && isscalar(t1))
    size_t0 = size(t0);
    size_t1 = size(t1);
    size_t0(end+1:length(size_t1)) = 1;
    size_t1(end+1:length(size_t0)) = 1;
    if ~(all(size_t0(1:2)==1) && all(size_t1(1:2)==1) && ...
            all(size_t0==size_t1 | size_t0==1 | size_t1==1))
        error('explde:Validation', ...
            ['In magnus_exp(fA,t0,t1,ord): ' ...
            'Invalid t0, t1 (must be size-compatible and ' ...
            'singleton in dimensions 1 and 2).']);
    end
end

switch ord
    case 2
        % Doc1 Eq 12
        Omega = (t1-t0).*fA(0.5*(t0+t1));
    case 4
        % Doc1 Eqs 13, 10 and Table 2; Doc0 Eq 56 (P_k)
        t_ctr = 0.5*(t0+t1);
        dt = 0.5*(t1-t0);
        u2 = 1/sqrt(3);
        u2_dt = u2*dt;
        A1 = fA(t_ctr-u2_dt);
        A2 = fA(t_ctr+u2_dt);
        g0 = dt.*(A1+A2);
        g1 = dt.*(u2*(A2-A1));
        Omega = g0-0.5*cmt(g0,g1);
    case 6
        % Doc1 Eqs 15, 10 and Table 2; Doc0 Eqs 80, 56 (P_k)
        t_ctr = 0.5*(t0+t1);
        dt = 0.5*(t1-t0);
        u3 = sqrt(3/5);
        u3_dt = u3*dt;
        A1 = fA(t_ctr-u3_dt);
        A2 = fA(t_ctr);
        A3 = fA(t_ctr+u3_dt);
        g0 = dt.*(8/9*A2+5/9*(A1+A3));
        g1 = dt.*(u3*5/9*(A3-A1));
        g2 = dt.*(-4/9*A2+2/9*(A1+A3));
        C1 = cmt(g0,g1);
        C2 = cmt(g0,g2+1/10*C1);
        C3 = cmt(g0-g2-3/10*C1,g1-1/6*C2);
        Omega = g0-1/2*C3;
    case 8
        % Doc1 Eqs 16, 10 and Table 2; Doc0 Eq 56 (P_k)
        t_ctr = 0.5*(t0+t1);
        dt = 0.5*(t1-t0);
        u3 = 3/7;
        u4 = 2/7*sqrt(6/5);
        [u3,u4] = deal(sqrt(u3-u4),sqrt(u3+u4));
        u3_dt = u3*dt;
        u4_dt = u4*dt;
        A1 = fA(t_ctr-u4_dt);
        A2 = fA(t_ctr-u3_dt);
        A3 = fA(t_ctr+u3_dt);
        A4 = fA(t_ctr+u4_dt);
        w3_2 = (18+sqrt(30))/36; % w3*2
        w4_2 = (18-sqrt(30))/36; % w4*2
        g0 = dt.*(w4_2*(A1+A4)+w3_2*(A2+A3));
        g1 = dt.*(w4_2*u4*(A4-A1)+w3_2*u3*(A3-A2));
        P2_3 = (3*u3^2-1)/2;
        P2_4 = (3*u4^2-1)/2;
        g2 = dt.*(w4_2*P2_4*(A1+A4)+w3_2*P2_3*(A2+A3));
        P3_3 = 3\(5*u3*P2_3-2*u3);
        P3_4 = 3\(5*u4*P2_4-2*u4);
        g3 = dt.*(w4_2*P3_4*(A4-A1)+w3_2*P3_3*(A3-A2));
        C1 = cmt(-1/2*g0+5/7*g2,g1-g3);
        C2 = cmt(g0-5/2*g2,5*g2-C1);
        C3 = cmt(1/2*g0-5/7*g2+3/14*C1,-2*g1+2*g3+1/21*C2);
        C4 = cmt(3/14*C1,-g1+7/2*g3);
        C5 = cmt(2*g0-5*g2+15/14*C1,10*g2-C3+C4);
        C6 = cmt(g0-C1+C3,3*g1+7*g3+1/2*C2-1/6*C5);
        Omega = g0+7/20*C3-1/20*C6;
    case 10
        % Doc2 Eq 3 and Table 4; Doc1 Eq 10; Doc0 Eq 56 (P_k)
        if isempty(x)
            % copied from magnus10coef.m
            x = [
                5/3
                0.0901217511054286179200006432152
                0.132144147153367034682119532902
                -1/5
                -0.472685390399902712406002183602
                3.23401783247220824419907849435
                1
                0.161670971415242270251705225531
                -0.525919202398768195923304595384
                -1.53079589966726923775289156382
                -0.910924396469683553350121373571
                -19/20
                -0.455192455511641898707317211627
                -7/160
                1/20
                0.999065729598492021595447046029
                -0.169685143102619456375901958438
                3.92051929889526396734420883485
                1
                0.555796983434714384253066760376
                -0.607334116992801989897737814155
                1
                -8.07420917786156075183923102445
                -10.0495069461621075880489787299
                80/3
                19/2
                12.1017461973869000839204686308
                -31.4531547349911685094163487849
                1/20
                0.275188457138167181875630064491
                -0.954326847042167934105954068326
                -0.784727265144758812578823778889
                1
                -1/4
                1.61045764376873307258367256828
                2.94807728420384234993029710640
                -11
                1
                -0.246517496492115672033127579776
                -2.90152105321322178449420125787
                -0.720231951037809406297947770664
                1
                -1/8
                -0.0834850346446151127595195646048
                0.596648858778570464309441892471
                -0.0848726298881828778243828476824
                0.189439293642050895990550879245
                -1
                1/11
                -2/9
                -1.22504203400172707389155574478
                2.93982839425824894217792481282
                1.89056056479215080737398254598
                1
                -9/4
                0.493726863995418585916367639513
                0.143371545985261839634589032722
                -1.10171588444014784576499503191
                -1.33521791511463352749639989683
                3
                -0.300277475990164551753003485261
                0.0422226321169889153042100461307
                -0.707315165639034053327374104654
                -1.00852021385304266703414054617
                0.522319462460385203547544946442
                -0.0395250285511817450646611919490
                1/2
                0.124884180986072604975597118288
                -0.501432555510506461635553953416
                0.335346238212012539516596859030
                0.359077663794944462844741747862
                0.837847882188203746637730745011
                1/5
                -0.0277515431965363360406864485263
                -0.0751536499542879518128684277851
                -0.594739381596017298987763177053
                0.345400850217906386578044663611
                -0.119041629440661798962689039945
                0.195879050217731400879397915960
                -0.219699011141684586146119923134
                0.263055372971027419203395240270
                ];
        end
        t_ctr = 0.5*(t0+t1);
        dt = 0.5*(t1-t0);
        u4 = 2*sqrt(10/7);
        u5 = 3\sqrt(5+u4);
        u4 = 3\sqrt(5-u4);
        u4_dt = u4*dt;
        u5_dt = u5*dt;
        A1 = fA(t_ctr-u5_dt);
        A2 = fA(t_ctr-u4_dt);
        A3 = fA(t_ctr);
        A4 = fA(t_ctr+u4_dt);
        A5 = fA(t_ctr+u5_dt);
        w3_2 = 128/225; % w3*2
        w4_2 = 13*sqrt(70);
        w5_2 = (322-w4_2)/900; % w5*2
        w4_2 = (322+w4_2)/900; % w4*2
        g0 = dt.*(w3_2*A3+w4_2*(A2+A4)+w5_2*(A1+A5));
        g1 = dt.*(w4_2*u4*(A4-A2)+w5_2*u5*(A5-A1));
        P2_3 = -1/2;
        P2_4 = (3*u4^2-1)/2;
        P2_5 = (3*u5^2-1)/2;
        g2 = dt.*(w3_2*P2_3*A3+w4_2*P2_4*(A2+A4)+w5_2*P2_5*(A1+A5));
        P3_4 = 3\(5*u4*P2_4-2*u4);
        P3_5 = 3\(5*u5*P2_5-2*u5);
        g3 = dt.*(w4_2*P3_4*(A4-A2)+w5_2*P3_5*(A5-A1));
        P4_3 = -3/4*P2_3;
        P4_4 = 4\(7*u4*P3_4-3*P2_4);
        P4_5 = 4\(7*u5*P3_5-3*P2_5);
        g4 = dt.*(w3_2*P4_3*A3+w4_2*P4_4*(A2+A4)+w5_2*P4_5*(A1+A5));
        % Note: x([7 19 22 33 38 42 54]) is all ones; these factors are
        % omitted. Also, x(48)==-1 and x(48) is omitted.
        C1 = cmt(x(1)*g0+x(2)*g2, ...
            x(3)*g1+x(4)*g3);
        C2 = cmt(x(5)*g0+x(6)*g2+C1, ...
            x(8)*g1+x(9)*g3);
        C3 = cmt(x(10)*g0+x(11)*g2, ...
            x(12)*g0+x(13)*g2+x(14)*g4+x(15)*C1);
        C4 = cmt(x(16)*g0+x(17)*g2+x(18)*g4+C2, ...
            x(20)*g1+x(21)*g3+C3);
        C5 = cmt(x(23)*g0+x(24)*g2+x(25)*g4, ...
            x(26)*g0+x(27)*g2+x(28)*g4+x(29)*C4);
        C6 = cmt(x(30)*g0+x(31)*g2+x(32)*g4+C2+x(34)*C4, ...
            x(35)*g1+x(36)*g3+x(37)*C3+C5);
        C7 = cmt(x(39)*g0+x(40)*g2+x(41)*g4+C2+x(43)*C4+x(44)*C6, ...
            x(45)*g0+x(46)*g2+x(47)*g4-C1+x(49)*C4+x(50)*C6);
        C8 = cmt(x(51)*g0+x(52)*g2+x(53)*g4+C1+x(55)*C2+x(56)*C4 ...
            +x(57)*C6, ...
            x(58)*g1+x(59)*g3+x(60)*C3+x(61)*C5+x(62)*C7);
        C9 = cmt(x(63)*g0+x(64)*g2+x(65)*g4+x(66)*C1+x(67)*C2+x(68)*C4 ...
            +x(69)*C6+x(70)*C8, ...
            x(71)*g1+x(72)*g3+x(73)*C3+x(74)*C5+x(75)*C7);
        Omega = g0+x(76)*C1+x(77)*C2+x(78)*C4+x(79)*C6+x(80)*C8+x(81)*C9;
    otherwise
        error('explde:Validation', ...
            ['In magnus_exp(fA,t0,t1,ord): ' ...
            'Invalid ord (must be 2, 4, 6, 8, or 10).']);
end
end % magnus_exp

function C = cmt(A,B) % commutator
C = pagemtimes(A,B)-pagemtimes(B,A);
end


