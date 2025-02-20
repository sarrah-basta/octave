////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2022 The Octave Project Developers
//
// See the file COPYRIGHT.md in the top-level directory of this
// distribution or <https://octave.org/copyright/>.
//
// This file is part of Octave.
//
// Octave is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, see
// <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////

#if ! defined (octave_mex_private_h)
#define octave_mex_private_h 1

#include "octave-config.h"

OCTAVE_NAMESPACE_BEGIN

extern OCTINTERP_API octave_value_list
call_mex (octave_mex_function& mex_fcn, const octave_value_list& args,
          int nargout_arg);

OCTAVE_NAMESPACE_END

#endif
