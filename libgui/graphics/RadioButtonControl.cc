/*

Copyright (C) 2011-2020 The Octave Project Developers

See the file COPYRIGHT.md in the top-level directory of this distribution
or <https://octave.org/COPYRIGHT.html/>.


This file is part of Octave.

Octave is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Octave is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Octave; see the file COPYING.  If not, see
<https://www.gnu.org/licenses/>.

*/

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <QRadioButton>

#include "ButtonGroup.h"
#include "RadioButtonControl.h"
#include "Container.h"
#include "QtHandlesUtils.h"

#include "octave-qobject.h"

namespace QtHandles
{

  RadioButtonControl*
  RadioButtonControl::create (octave::base_qobject& oct_qobj,
                              octave::interpreter& interp,
                              const graphics_object& go)
  {
    Object *parent = parentObject (interp, go);

    if (parent)
      {
        Container *container = parent->innerContainer ();

        if (container)
          return new RadioButtonControl (oct_qobj, interp, go,
                                         new QRadioButton (container));
      }

    return nullptr;
  }

  RadioButtonControl::RadioButtonControl (octave::base_qobject& oct_qobj,
                                          octave::interpreter& interp,
                                          const graphics_object& go,
                                          QRadioButton *radio)
    : ButtonControl (oct_qobj, interp, go, radio)
  {
    Object *parent = parentObject (interp, go);
    ButtonGroup *btnGroup = dynamic_cast<ButtonGroup *>(parent);
    if (btnGroup)
      btnGroup->addButton (radio);

    radio->setAutoFillBackground (true);
    radio->setAutoExclusive (false);
  }

  RadioButtonControl::~RadioButtonControl (void)
  { }

};
