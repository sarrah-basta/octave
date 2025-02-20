########################################################################
##
## Copyright (C) 1996-2022 The Octave Project Developers
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
## @deftypefn {} {@var{q} =} empirical_inv (@var{x}, @var{data})
## For each element of @var{x}, compute the quantile (the inverse of the CDF)
## at @var{x} of the empirical distribution obtained from the
## univariate sample @var{data}.
## @end deftypefn

function q = empirical_inv (x, data)

  if (nargin != 2)
    print_usage ();
  endif

  if (! isvector (data))
    error ("empirical_inv: DATA must be a vector");
  endif

  q = discrete_inv (x, data, ones (size (data)));

endfunction


%!shared x,v,y
%! x = [-1 0 0.1 0.5 1 2];
%! v = 0.1:0.2:1.9;
%! y = [NaN v(1) v(1) v(end/2) v(end) NaN];
%!assert (empirical_inv (x, v), y, eps)

## Test class of input preserved
%!assert (empirical_inv ([x, NaN], v), [y, NaN], eps)
%!assert (empirical_inv (single ([x, NaN]), v), single ([y, NaN]), eps)
%!assert (empirical_inv ([x, NaN], single (v)), single ([y, NaN]), eps)

## Test input validation
%!error <Invalid call> empirical_inv ()
%!error <Invalid call> empirical_inv (1)
%!error empirical_inv (1, ones (2))
