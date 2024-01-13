function [x, fval, info, output] = ACM (fun, x0)

  %% Get default options if requested.
  if (nargin == 1 && ischar (fun) && strcmp (fun, 'defaults'))
    x = optimset ("MaxIter", Inf, "MaxFunEvals", Inf, "TolX", 1e-8, "OutputFcn", [], "FunValCheck", "off");
    return;
  end

  if (nargin < 2 || nargin > 3)
    print_usage ();
  end

  if (ischar (fun))
    fun = str2func (fun, "global");
  end

  %% TODO
  %% displev = optimget (options, "Display", "notify");
  funvalchk = strcmpi (optimget (options, "FunValCheck", "off"), "on");
  outfcn = optimget (options, "OutputFcn");
  tolx = optimget (options, "TolX", 1e-8);
  maxiter = optimget (options, "MaxIter", Inf);
  maxfev = optimget (options, "MaxFunEvals", Inf);

  persistent mu ;
  mu = 0.5;

  if (funvalchk)
    %% Replace fun with a guarded version.
    fun = @(x) guarded_eval (fun, x);
  end

  %% The default exit flag if exceeded number of iterations.
  info = 0;
  niter = 0;
  nfev = 0;

  fb = NaN;
  b = fb;
  fa = b;
  a = fa;
  fval = a;
  x = fval;
 
  eps = eps (class (x0));

  %% Prepare...
  a = x0(1);
  fa = fun (a);
  nfev = 1;
  if (length (x0) > 1)
    b = x0(2);
    fb = fun (b);
    nfev = nfev+  1;
  else
    %% Try to get b.
    if (a == 0)
      aa = 1;
    else
      aa = a;
    end
    for b = [0.9*aa, 1.1*aa, aa-1, aa+1, 0.5*aa 1.5*aa, -aa, 2*aa, -10*aa, 10*aa]
      fb = fun (b); nfev += 1;
      if (sign (fa) * sign (fb) <= 0)
        break;
      end
    end
  end

  if (b < a)
    u = a;
    a = b;
    b = u;

    fu = fa;
    fa = fb;
    fb = fu;
  end

  if (~ (sign (fa) * sign (fb) <= 0))
    error ("fzero:bracket", "fzero: not a valid initial bracketing");
  end

  slope0 = (fb - fa) / (b - a);

  if (fa == 0)
    b = a;
    fb = fa;
  elseif (fb == 0)
    a = b;
    fa = fb;
  end

  itype = 1;

  if (abs (fa) < abs (fb))
    u = a; fu = fa;
  else
    u = b; fu = fb;
  end

  d = e;
  e = u;
  fd = fe;
  fe = fu;
  mba = mu*(b - a);
  while (niter < maxiter && nfev < maxfev)
    switch (itype)
    case 1
      %% The initial test.
      if (b - a <= 2*(2 * abs (u) * eps + tolx))
        x = u; fval = fu;
        info = 1;
        break;
      end
      if (abs (fa) <= 1e3*abs (fb) && abs (fb) <= 1e3*abs (fa))
        %% Secant step.
        c = u - (a - b) / (fa - fb) * fu;
      else
        %% Bisection step.
        c = 0.5*(a + b);
      end
      d = u; fd = fu;
      itype = 5;
    case {2, 3}
      l = length (unique ([fa, fb, fd, fe]));
      if (l == 4)
        %% Inverse cubic interpolation.
        q11 = (d - e) * fd / (fe - fd);
        q21 = (b - d) * fb / (fd - fb);
        q31 = (a - b) * fa / (fb - fa);
        d21 = (b - d) * fd / (fd - fb);
        d31 = (a - b) * fb / (fb - fa);
        q22 = (d21 - q11) * fb / (fe - fb);
        q32 = (d31 - q21) * fa / (fd - fa);
        d32 = (d31 - q21) * fd / (fd - fa);
        q33 = (d32 - q22) * fa / (fe - fa);
        c = a + q31 + q32 + q33;
      end
      if (l < 4 || sign (c - a) * sign (c - b) > 0)
        %% Quadratic interpolation + newton.
        a0 = fa;
        a1 = (fb - fa)/(b - a);
        a2 = ((fd - fb)/(d - b) - a1) / (d - a);
        %% Modification 1: this is simpler and does not seem to be worse.
        c = a - a0/a1;
        if (a2 ~= 0)
          c = a - a0/a1;
          for i = 1:itype
            pc = a0 + (a1 + a2*(c - b))*(c - a);
            pdc = a1 + a2*(2*c - a - b);
            if (pdc == 0)
              c = a - a0/a1;
              break;
            end
            c -= pc/pdc;
          end
        end
      end
      itype += 1;
    case 4
      %% Double secant step.
      c = u - 2*(b - a)/(fb - fa)*fu;
      %% Bisect if too far.
      if (abs (c - u) > 0.5*(b - a))
        c = 0.5 * (b + a);
      end
      itype = 5;
    case 5
      %% Bisection step.
      c = 0.5 * (b + a);
      itype = 2;
        end

    %% Don't let c come too close to a or b.
    delta = 2*0.7*(2 * abs (u) * eps + tolx);
    if ((b - a) <= 2*delta)
      c = (a + b)/2;
    else
      c = max (a + delta, min (b - delta, c));
    end

    %% Calculate new point.
    x = c;
    fval = fc;
    fc = fun (c);
    niter = niter +1 ;
    nfev = nfev +1;

    %% Modification 2: skip inverse cubic interpolation if
    %% nonmonotonicity is detected.
    if (sign (fc - fa) * sign (fc - fb) >= 0)
      %% The new point broke monotonicity.
      %% Disable inverse cubic.
      fe = fc;
    else
      e = d; fe = fd;
    end

    %% Bracketing.
    if (sign (fa) * sign (fc) < 0)
      d = b; fd = fb;
      b = c; fb = fc;
    elseif (sign (fb) * sign (fc) < 0)
      d = a; fd = fa;
      a = c; fa = fc;
    elseif (fc == 0)
      a = b;
      b = c;
      fa = fb;
      fb = fc;
      info = 1;
      break;
    else
      %% This should never happen.
      error ("fzero:bracket", "fzero: zero point is not bracketed");
    end

    %% If there's an output function, use it now.
    if (outfcn)
      optv.funccount = nfev;
      optv.fval = fval;
      optv.iteration = niter;
      if (outfcn (x, optv, "iter"))
        info = -1;
        break;
      end
    end

    if (abs (fa) < abs (fb))
      u = a; fu = fa;
    else
      u = b; fu = fb;
    end
    if (b - a <= 2*(2 * abs (u) * eps + tolx))
      info = 1;
      break;
    end

    %% Skip bisection step if successful reduction.
    if (itype == 5 && (b - a) <= mba)
      itype = 2;
    end
    if (itype == 2)
      mba = mu * (b - a);
    end
    end

  %% Check solution for a singularity by examining slope
  if (info == 1)
    if ((b - a) ~= 0 && abs ((fb - fa)/(b - a) / slope0) > max (1e6, 0.5/(eps+tolx)))
      info = -5;
    end
  end

  output.iterations = niter;
  output.funcCount = nfev;
  output.bracketx = [a, b];
  output.brackety = [fa, fb];

    end

%% An assistant function that evaluates a function handle and checks for
%% bad results.
function fx = guarded_eval (fun, x)
  fx = fun (x);
  fx = fx(1);
  if (~ isreal (fx))
    error ("fzero:notreal", "fzero: non-real value encountered");
  elseif (isnan (fx))
    error ("fzero:isnan", "fzero: NaN value encountered");
  end
end