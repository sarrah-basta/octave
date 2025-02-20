########################################################################
##
## Copyright (C) 2009-2022 The Octave Project Developers
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
## @deftypefn {} {@var{names} =} __all_opts__ (@dots{})
## Internal function.
##
## Query all options from all known optimization functions and return a
## list of possible values.
## @end deftypefn

function names = __all_opts__ (varargin)

  ## This variable is filled by the auto-generated PKG_ADD script at
  ## Octave startup.
  persistent saved_names = {};

  ## do not clear this function
  mlock ();

  ## guard against recursive calls.
  persistent recursive = false;

  if (recursive)
    names = {};
  elseif (nargin == 0)
    names = saved_names;
  else
    ## Query all options from all known functions.  These may call optimset,
    ## which will in turn call us, but we won't answer.
    recursive = true;
    names = saved_names;
    for i = 1:nargin
      try
        opts = optimset (varargin{i});
        fcn = fieldnames (opts).';
        names = [names, fcn];
      catch
        ## throw the error as a warning.
        warning (lasterr ());
      end_try_catch
    endfor
    names = unique (names);
    [lnames, idx] = unique (tolower (names));
    if (length (lnames) < length (names))
      ## This is bad.
      error ("__all_opts__: duplicate options with inconsistent case");
    else
      names = names(idx);
    endif
    saved_names = names;
    recursive = false;
  endif

endfunction


## No test needed for internal helper function.
%!assert (1)
