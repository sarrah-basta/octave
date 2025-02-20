########################################################################
##
## Copyright (C) 2013-2022 The Octave Project Developers
##
## See the file COPYRIGHT.md in the top-level directory of this
## distribution or <https://octave.org/copyright/>.
##
## This file is part of Octave.
##
## Octave is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <https://www.gnu.org/licenses/>.
##
########################################################################

## -*- texinfo -*-
## @deftypefn {} {@var{solution} =} integrate_adaptive (@var{@@stepper}, @var{order}, @var{@@fcn}, @var{tspan}, @var{x0}, @var{options})
##
## This function file can be called by an ODE solver function in order to
## integrate the set of ODEs on the interval @var{[t0, t1]} with an adaptive
## timestep.
##
## The function returns a structure @var{solution} with two fields: @var{t}
## and @var{y}.  @var{t} is a column vector and contains the time stamps.
## @var{y} is a matrix in which each column refers to a different unknown
## of the problem and the row number is the same as the @var{t} row number.
## Thus, each row of the matrix @var{y} contains the values of all unknowns at
## the time value contained in the corresponding row in @var{t}.
##
## The first input argument must be a function handle or inline function
## representing the stepper, i.e., the function responsible for step-by-step
## integration.  This function discriminates one method from the others.
##
## The second input argument is the order of the stepper.  It is needed
## to compute the adaptive timesteps.
##
## The third input argument is a function handle or inline function that
## defines the ODE:
##
## @ifhtml
##
## @example
## @math{y' = f(t,y)}
## @end example
##
## @end ifhtml
## @ifnothtml
## @math{y' = f(t,y)}.
## @end ifnothtml
##
## The fourth input argument is the time vector which defines the integration
## interval, i.e., @var{[tspan(1), tspan(end)]} and all intermediate elements
## are taken as times at which the solution is required.
##
## The fifth argument represents the initial conditions for the ODEs and the
## last input argument contains some options that may be needed for the
## stepper.
##
## @end deftypefn

function solution = integrate_adaptive (stepper, order, fcn, tspan, x0,
                                        options)

  fixed_times = numel (tspan) > 2;

  t_new = t_old = t = tspan(1);
  x_new = x_old = x = x0(:);

  ## Get first initial timestep
  dt = options.InitialStep;
  if (isempty (dt))
    dt = starting_stepsize (order, fcn, t, x,
                            options.AbsTol, options.RelTol,
                            strcmp (options.NormControl, "on"),
                            options.funarguments);
  endif

  dir = options.direction;
  dt = dir * min (abs (dt), options.MaxStep);

  options.comp = 0.0;

  ## Factor multiplying the stepsize guess
  facmin = 0.8;
  facmax = 1.5;
  fac = 0.38^(1/(order+1));  # formula taken from Hairer

  ## Initialize the OutputFcn
  if (options.haveoutputfunction)
    if (! isempty (options.OutputSel))
      solution.retout = x(options.OutputSel,end);
    else
      solution.retout = x;
    endif
    feval (options.OutputFcn, tspan, solution.retout, "init",
           options.funarguments{:});
  endif

  ## Initialize the EventFcn
  have_EventFcn = false;
  if (! isempty (options.Events))
    have_EventFcn  = true;
    ode_event_handler (options.Events, tspan(1), x,
                       "init", options.funarguments{:});
  endif

  if (options.havenonnegative)
    nn = options.NonNegative;
  endif

  solution.cntloop = 0;
  solution.cntcycles = 0;
  solution.cntsave = 2;
  solution.unhandledtermination = true;
  ireject = 0;

  NormControl = strcmp (options.NormControl, "on");
  k_vals = [];
  iout = istep = 1;

  while (dir * t_old < dir * tspan(end))

    ## Compute integration step from t_old to t_new = t_old + dt
    [t_new, options.comp] = kahan (t_old, options.comp, dt);
    [t_new, x_new, x_est, new_k_vals] = ...
      stepper (fcn, t_old, x_old, dt, options, k_vals, t_new);

    solution.cntcycles += 1;

    if (options.havenonnegative)
      x_new(nn, end) = abs (x_new(nn, end));
      x_est(nn, end) = abs (x_est(nn, end));
    endif

    err = AbsRel_norm (x_new, x_old, options.AbsTol, options.RelTol,
                       NormControl, x_est);

    ## Accept solution only if err <= 1.0
    if (err <= 1)

      solution.cntloop += 1;
      ireject = 0;              # Clear reject counter

      ## if output time steps are fixed
      if (fixed_times)

        t_caught = find ((dir * tspan(iout:end) > dir * t_old)
                         & (dir * tspan(iout:end) <= dir * t_new));
        t_caught = t_caught + iout - 1;

        if (! isempty (t_caught))
          t(t_caught) = tspan(t_caught);
          iout = max (t_caught);
          x(:, t_caught) = ...
            runge_kutta_interpolate (order, [t_old t_new], [x_old x_new],
                                     t(t_caught), new_k_vals, dt, fcn,
                                     options.funarguments);

          istep += 1;

          ## Call Events function only if a valid result has been found.
          ## Stop integration if eventbreak is true.
          if (have_EventFcn)
            break_loop = false;
            for idenseout = 1:numel (t_caught)
              id = t_caught(idenseout);
              td = t(id);
              solution.event = ...
                ode_event_handler (options.Events, t(id), x(:, id), [],
                                   options.funarguments{:});
              if (! isempty (solution.event{1}) && solution.event{1} == 1)
                t(id) = solution.event{3}(end);
                t = t(1:id);
                x(:, id) = solution.event{4}(end, :).';
                x = x(:,1:id);
                solution.unhandledtermination = false;
                break_loop = true;
                break;
              endif
            endfor
            if (break_loop)
              break;
            endif
          endif

          ## Call OutputFcn only if a valid result has been found.
          ## Stop integration if function returns true.
          if (options.haveoutputfunction)
            cnt = options.Refine + 1;
            approxtime = linspace (t_old, t_new, cnt);
            approxvals = interp1 ([t_old, t(t_caught), t_new],
                                  [x_old, x(:, t_caught), x_new] .',
                                  approxtime, "linear") .';
            if (isvector (approxvals) && ! isscalar (approxtime))
              approxvals = approxvals.';
            endif
            if (! isempty (options.OutputSel))
              approxvals = approxvals(options.OutputSel, :);
            endif
            stop_solve = false;
            for ii = 1:numel (approxtime)
              stop_solve = feval (options.OutputFcn,
                                  approxtime(ii), approxvals(:, ii), [],
                                  options.funarguments{:});
              if (stop_solve)
                break;  # break from inner loop
              endif
            endfor
            if (stop_solve)  # Leave main loop
              solution.unhandledtermination = false;
              break;
            endif
          endif

        endif

      else   # not fixed times

        t(++istep)  = t_new;
        x(:, istep) = x_new;
        iout = istep;

        ## Call Events function only if a valid result has been found.
        ## Stop integration if eventbreak is true.
        if (have_EventFcn)
          solution.event = ...
            ode_event_handler (options.Events, t(istep), x(:, istep), [],
                               options.funarguments{:});
          if (! isempty (solution.event{1}) && solution.event{1} == 1)
            t(istep) = solution.event{3}(end);
            x(:, istep) = solution.event{4}(end, :).';
            solution.unhandledtermination = false;
            break;
          endif
        endif

        ## Call OutputFcn only if a valid result has been found.
        ## Stop integration if function returns true.
        if (options.haveoutputfunction)
          cnt = options.Refine + 1;
          approxtime = linspace (t_old, t_new, cnt);
          approxvals = interp1 ([t_old, t_new],
                                [x_old, x_new] .',
                                approxtime, "linear") .';
          if (isvector (approxvals) && ! isscalar (approxtime))
            approxvals = approxvals.';
          endif
          if (! isempty (options.OutputSel))
            approxvals = approxvals(options.OutputSel, :);
          endif
          stop_solve = false;
          for ii = 1:numel (approxtime)
            stop_solve = feval (options.OutputFcn,
                                approxtime(ii), approxvals(:, ii), [],
                                options.funarguments{:});
            if (stop_solve)
              break;  # break from inner loop
            endif
          endfor
          if (stop_solve)  # Leave main loop
            solution.unhandledtermination = false;
            break;
          endif
        endif

      endif

      ## move to next time-step
      t_old = t_new;
      x_old = x_new;
      k_vals = new_k_vals;

    else  # error condition

      ireject += 1;

      ## Stop solving if, in the last 5,000 steps, no successful valid
      ## value has been found.
      if (ireject >= 5_000)
        error (["integrate_adaptive: Solving was not successful. ", ...
                " The iterative integration loop exited at time", ...
                " t = %f before the endpoint at tend = %f was reached. ", ...
                " This happened because the iterative integration loop", ...
                " did not find a valid solution at this time stamp. ", ...
                " Try to reduce the value of 'InitialStep' and/or", ...
                " 'MaxStep' with the command 'odeset'.\n"],
               t_old, tspan(end));
      endif

    endif

    ## Compute next timestep, formula taken from Hairer
    err += eps;  # avoid divisions by zero
    dt *= min (facmax, max (facmin, fac * (1 / err)^(1 / (order + 1))));
    dt = dir * min (abs (dt), options.MaxStep);
    if (! (abs (dt) > eps (t(end))))
      break;
    endif

    ## Make sure we don't go past tpan(end)
    dt = dir * min (abs (dt), abs (tspan(end) - t_old));

  endwhile

  ## Check if integration of the ode has been successful
  if (dir * t(end) < dir * tspan(end))
    if (solution.unhandledtermination == true)
      warning ("integrate_adaptive:unexpected_termination",
               [" Solving was not successful. ", ...
                " The iterative integration loop exited at time", ...
                " t = %f before the endpoint at tend = %f was reached. ", ...
                " This may happen if the stepsize becomes too small. ", ...
                " Try to reduce the value of 'InitialStep'", ...
                " and/or 'MaxStep' with the command 'odeset'."],
               t(end), tspan(end));
    else
      warning ("integrate_adaptive:unexpected_termination",
               ["Solver was stopped by a call of 'break'", ...
                " in the main iteration loop at time", ...
                " t = %f before the endpoint at tend = %f was reached. ", ...
                " This may happen because the @odeplot function", ...
                " returned 'true' or the @event function returned 'true'."],
               t(end), tspan(end));
    endif
  endif

  ## Set up return structure
  solution.t = t(:);
  solution.x = x.';

endfunction
