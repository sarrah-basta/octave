########################################################################
##
## Copyright (C) 1995-2022 The Octave Project Developers
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
## @deftypefn {} {@var{str} =} date ()
## Return the current date as a character string in the form DD-MMM-YYYY@.
##
## For example:
##
## @example
## @group
## date ()
##   @result{} 20-Aug-1993
## @end group
## @end example
## @seealso{now, clock, datestr, localtime}
## @end deftypefn

function str = date ()

  str = strftime ("%d-%b-%Y", localtime (time ()));

endfunction


%!assert (strcmp (date (), strftime ("%d-%b-%Y", localtime (time ()))))
