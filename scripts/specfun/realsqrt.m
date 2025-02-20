########################################################################
##
## Copyright (C) 2008-2022 The Octave Project Developers
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
## @deftypefn {} {@var{y} =} realsqrt (@var{x})
## Return the real-valued square root of each element of @var{x}.
##
## If any element results in a complex return value @code{realsqrt} aborts and
## issues an error.
## @seealso{sqrt, realpow, reallog}
## @end deftypefn

function y = realsqrt (x)

  if (nargin < 1)
    print_usage ();
  elseif (iscomplex (x) || any (x(:) < 0))
    error ("realsqrt: produced complex result");
  else
    y = sqrt (x);
  endif

endfunction


%!assert (realsqrt (1:5), sqrt (1:5))
%!test
%! x = rand (10, 10);
%! assert (realsqrt (x), sqrt (x));

%!error <Invalid call> realsqrt ()
%!error <produced complex result> realsqrt (-1)
