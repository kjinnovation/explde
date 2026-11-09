function [F,F_] = explde(A,t,ord,F0)
%EXPLDE Exponential-matrix solution of a system of linear differential
%equations.
%   This function numerically solves a first-order, linear, and homogeneous
%   matrix differential equation with either constant or non-constant
%   coefficients. The solution employs the Magnus exponent method. (See
%   help magnus_exp.)
%   Code and documentation: https://github.com/kjinnovation/explde
%
%   The differential equation is of the form
%     F'(t) = A(t)*F(t)
%   with the initial value of F(t) specified at some point t = t0. A and F
%   are matrix functions of scalar argument t. The solution has the form
%     F(t) = explde(A,[t0,t])*F(t0)
%   Alternatively,
%     F(t) = explde(A,[t0,t],ord,F(t0))
%   (This is more efficient if F has fewer columns than rows. The first
%   form would be more efficient if F has more columns than rows.)
%
%   For a constant matrix A and scalar t, explde(A,[0,t]) is functionally
%   equivalent to expm(A*t), and explde(A) is equivalent to expm(A). For a
%   function handle A, explde(A,[0,t]) is expm(Omega(t)), where Omega(t) is
%   a quadrature approximation to the Magnus exponent (Gauss-Legendre
%   quadrature) with approximation order ord = 2, 4, 6, 8, or 10. (The
%   default is 6.)
%
%   The above code example performs a single integration step from 0 to t.
%   A multi-step integration can be performed by making the second explde
%   argument a vector defining the integration intervals, e.g.,
%     F(t(end)) = explde(A,t)*F(t(1))
%
%   A, t, and F can generally be multidimensional arrays that are
%   size-compatible in dimensions 3, 4, ... for paged operation.
%   "Size-compatible" means size-matched with implicit size expansion in
%   singleton dimensions:
%   https://www.mathworks.com/help/matlab/matlab_prog/compatible-array-sizes-for-basic-operations.html
%   For constant array A, explde(A) is equivalent to pageexpm(A).
%
%   pageexpm.m and magnus_exp.m must be on the MATLAB path; otherwise
%   explde will fail.
%   https://www.mathworks.com/matlabcentral/fileexchange/107959-pageexpm
%
% Syntax:
%    F = explde(A);
%    F = explde(A,t,ord,F0);
%    [F,F_] = explde(A,...);
%
% Inputs:
%
%   A: one of the following:
%     * a size-[M,M,...], non-empty numeric array (isfloat(A), M>0). This
%     is for the constant-coefficient case. A general multidimensional
%     array A (ndims(A)>2) is used for parallelizing multiple integrations.
%     * a function handle, which accepts a size-[1,1,...] argument t and
%     returns a size-[M,M,...] numeric array A(t).
%
%   t (optional, default = [0,1]): a numeric, real size-[1,Nt] vector or
%   size-[1,Nt,...] array with Nt>=2. t must be size-compatible with A in
%   dimensions 3, 4, ... . For vector t, F is integrated from t(1) to
%   t(end). The t vector defines the algorithm integration steps, and F is
%   calculated at all t points. A general multidimensional array t
%   (ndims(t)>2) is used for parallelizing multiple integrations.
%
%   ord (optional, default = 6): either 2, 4, 6, 8, or 10; approximation
%   order for the Magnus exponent. The approximation error for t step dt is
%   of order dt^(ord+1). (The number of steps is generally inversely
%   proportional to dt, so the cumulative integration error over all steps
%   is of order dt^ord.)
%
%   F0 (optional, default = eye(M,M)): a size-[M,N,...], non-empty numeric
%   array (isfloat(F0), N>0). F0 is the initial value of F at t(1) (or at
%   t(1,1,:,:,...) for non-vector t).
%
% Outputs:
%
%   F: size-[M,N,...] array, the equation solution at t(end) (or
%   t(1,end,:,:,...) for non-vector t).
%
%   F_ (optional): cell vector, size-[1,Nt]. F_{j} is the equation solution
%   at t(j) (or t(1,j,:,:,...) for non-vector t).

% version 09-Nov-2026, Ken Johnson, ken@kjinnovation.com

% Argument defaults and validation:
narginchk(1,4)
if nargin<2
    t = [0,1];
elseif ~(isfloat(t) && isreal(t) && ~isempty(t) ...
        && size(t,1)==1 && size(t,2)>=2)
    error('explde:Validation', ...
        'In explde(A,t,...): Invalid t (type or size).');
end
size_ = size(t);
Nt = size_(2);
size_(2) = 1; % t page sizes
if nargin<3
    ord = 6;
elseif ~(isscalar(ord) && any(ord==2:2:10))
    error('explde:Validation', ...
        'In explde(A,t,ord,...): Invalid ord.');
end
fA = [];
if isa(A,'function_handle')
    fA = A;
    t0 = reshape(t(1,1,:),size_);
    A = fA(t0);
end
if ~(isfloat(A) && ~isempty(A) && size(A,1)==size(A,2))
    error('explde:Validation', ...
        'In explde(A,...): Invalid A (type or size).');
end
if nargin<4
    F0 = []; % alias for eye(M,M)
elseif ~(isfloat(F0) && ~isempty(F0) && size(F0,1)==size(A,2))
    error('explde:Validation', ...
        'In explde(A,t,ord,F0): Invalid F0 (type or size).');
end
if ismatrix(A)+ismatrix(t)+ismatrix(F0)<=1
    % At least 2 of A, t, and F0 are paged. Check size compatibility in
    % dimensions 3, 4, ... .
    n = max([ndims(A),ndims(t),ndims(F0)]); % n>=3
    pagesize_A = size(A,3:n);
    pagesize_t = size(t,3:n);
    pagesize_F0 = size(F0,3:n);
    if ~ismatrix(A) && ~ismatrix(t)
        j = find(pagesize_A~=1 & pagesize_t~=1 ...
            & pagesize_A~=pagesize_t,1);
        if ~isempty(j)
            error('explde:Validation', ...
                ['In explde(A,t,...): A and t are not ' ...
                'size-compatible in dimension ' num2str(j+2) '.']);
        end
    end
    if ~ismatrix(A) && ~ismatrix(F0)
        j = find(pagesize_A~=1 & pagesize_F0~=1 ...
            & pagesize_A~=pagesize_F0,1);
        if ~isempty(j)
            error('explde:Validation', ...
                ['In explde(A,t,ord,F0): A and F0 are not ' ...
                'size-compatible in dimension ' num2str(j+2) '.']);
        end
    end
    if ~ismatrix(t) && ~ismatrix(F0)
        j = find(pagesize_t~=1 & pagesize_F0~=1 ...
            & pagesize_t~=pagesize_F0,1);
        if ~isempty(j)
            error('explde:Validation', ...
                ['In explde(A,t,ord,F0): t and F0 are not ' ...
                'size-compatible in dimension ' num2str(j+2) '.']);
        end
    end
end
if nargout<2
    F_ = [];
else
    F_ = cell(1,Nt);
    if ~isempty(F0)
        F_{1} = F0;
    else
        F_{1} = eye(size(A,1));
    end
end
F = F0;
t1 = reshape(t(1,1,:),size_);
if isempty(fA) && isempty(F_)
    t0 = t1;
    t1 = reshape(t(1,end,:),size_);
    dt = t1-t0; % full integration range
    U = pageexpm(A.*dt);
    if isempty(F) % alias for eye(M,M)
        F = U;
    else
        F = pagemtimes(U,F);
    end
else
    for j = 1:Nt-1
        % Define the integration interval ([t0,t1]).
        t0 = t1;
        t1 = reshape(t(1,j+1,:),size_);
        dt = t1-t0; % integration step
        if isempty(fA)
            % constant-A case
            U = pageexpm(A.*dt);
        else
            % variable-A case
            Omega = magnus_exp(fA,t0,t1,ord);
            U = pageexpm(Omega);
        end
        if isempty(F) % alias for eye(M,M)
            F = U;
        else
            F = pagemtimes(U,F);
        end
        if ~isempty(F_)
            F_{j+1} = F; %#ok<AGROW>
        end
    end
end
end % explde
