########################################################################
##
## Copyright (C) 1994-2022 The Octave Project Developers
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
## @deftypefn {} {@var{y} =} sinc (@var{x})
## Compute the sinc function.
##
## Return
## @tex
## $ \sin (\pi x)/(\pi x)$.
## @end tex
## @ifnottex
## sin (pi*x) / (pi*x).
## @end ifnottex
## @end deftypefn

function y = sinc (x)

  if (nargin != 1)
    print_usage ();
  endif

  y = ones (size (x));

  i = (x != 0);

  if (any (i(:)))
    t = pi * x(i);
    y(i) = sin (t) ./ t;
  endif

endfunction


%!assert (sinc (0), 1)
%!assert (sinc (1), 0,1e-6)
%!assert (sinc (1/2), 2/pi, 1e-6)

## Test input validation
%!error <Invalid call> sinc ()
