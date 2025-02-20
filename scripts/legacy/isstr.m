########################################################################
##
## Copyright (C) 2003-2022 The Octave Project Developers
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
## @deftypefn {} {@var{tf} =} isstr (@var{x})
## This function is obsolete.  Use @code{ischar} instead.
## @seealso{ischar}
## @end deftypefn

## At one time, Matlab docs stated that this function is obsolete and would be
## removed in some future version.  Now users are told that it should be
## avoided, but there is no mention of possible future removal.

function tf = isstr (varargin)

  persistent warned = false;
  if (! warned)
    warned = true;
    warning ("Octave:legacy-function",
             "isstr is obsolete; please use ischar instead");
  endif

  tf = ischar (varargin{:});

endfunction
