# Copyright (C) 1996,1999 Auburn University.  All Rights Reserved
#
# This file is part of Octave. 
#
# Octave is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the 
# Free Software Foundation; either version 2, or (at your option) any 
# later version. 
# 
# Octave is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
# for more details.
# 
# You should have received a copy of the GNU General Public License 
# along with Octave; see the file COPYING.  If not, write to the Free 
# Software Foundation, 59 Temple Place, Suite 330, Boston, MA 02111 USA. 
 
function sys = syssub(...)
# 
# sys = syssub(Gsys{,Hsys, ...})
#
# returns transfer function sys = Gsys - Hsys - ...
#
# Method: sysgroup used to connect systems in parallel
# The input vector is connected to all systems; the outputs of all
# systems are connected to a summing junction with the first system's
# outputs added, all other outputs subtracted.
# Returned system input/output signal names are those of Gsys.
# Example: syssub(Gsys, Hsys) results in
#
#                  ________
#             ----|  Gsys  |---
#        u   |    ----------  +|         
#        -----                (_)----> y
#            |     ________   -|
#             ----|  Hsys  |---
#                  --------

# Written by John Ingram July 1996
# updated for variable numbers of input arguments by July 1999 A. S. Hodel

  save_val = implicit_str_to_num_ok;	# save for later
  implicit_str_to_num_ok = 1;

  if(nargin < 1)
    usage("syssub: sys = syssub(Gsys{,Hsys,...})");
  endif

  # collect all arguments
  arglist = list();
  va_start();
  for kk=1:nargin
    arglist(kk) = va_arg();
    if(!is_struct(nth(arglist,kk)))
      error("syssub: argument %d is not a data structure",kk);
    endif
  endfor           

  # check system dimensions
  [n,nz,mg,pg,Gyd] = sysdimensions(nth(arglist,1));
  for kk=2:nargin
    [n,nz,mh,ph,Hyd] = sysdimensions(nth(arglist,kk));
    if(mg != mh)
      error("arg 1 has %d inputs; arg %d has vs %d inputs",mg,kk,mh);
    elseif(pg != ph)
      error("arg 1 has %d outputs; arg %d has vs %d outputs",pg,kk,ph);
    elseif(norm(Gyd - Hyd))
      warning("cannot add a discrete output to a continuous output");
      error("Output type mismatch: arguments 1 and %d\n",kk);
    endif
  endfor

  # perform the subtract
  if(nargin == 2)
    Gsys = nth(arglist,1);   Hsys = nth(arglist,2);
    if( strcmp(sysgettype(Gsys),"tf") | strcmp(sysgettype(Hsys),"tf") )
      # see if subtracting  transfer functions with identical denominators
      [Gnum,Gden,GT,Gin,Gout] = sys2tf(Gsys);
      [Hnum,Hden,HT,Hin,Hout] = sys2tf(Hsys);
      if(length(Hden) == length(Gden) )
        if( (Hden == Gden) & (HT == GT) )
          sys = tf2sys(Gnum+Hnum,Gden,GT,Gin,Gout);
          return
        endif
        # if not, we go on and do the usual thing...
      endif
    endif
  
    # make sure in ss form
    Gsys = sysupdate(Gsys,"ss");
    Hsys = sysupdate(Hsys,"ss");
  
    # change signal names to avoid warning messages from sysgroup
    Gsys = syssetsignals(Gsys,"in",sysdefioname(length(Gin),"Gin_u"));
    Gsys = syssetsignals(Gsys,"out",sysdefioname(length(Gout),"Gout_u"));
    Hsys = syssetsignals(Hsys,"in",sysdefioname(length(Hin),"Hin_u"));
    Hsys = syssetsignals(Hsys,"out",sysdefioname(length(Hout),"Hout_u"));
    
    sys = sysgroup(Gsys,Hsys);
  
    eyin = eye(mg);
    eyout = eye(pg);
  
    sys = sysscale(sys,[eyout -eyout],[eyin;eyin],Gout,Gin);
  
  else
    # multiple systems (or a single system); combine together one by one
    sys = nth(arglist,1);
    for kk=2:length(arglist)
      sys = syssub(sys,nth(arglist,kk));
    endfor
  endif
  
endfunction
