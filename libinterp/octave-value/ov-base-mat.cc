////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 1996-2022 The Octave Project Developers
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

// This file should not include config.h.  It is only included in other
// C++ source files that should have included config.h before including
// this file.

#include <ostream>
#include <sstream>

#include "Array-util.h"

#include "Cell.h"
#include "errwarn.h"
#include "ovl.h"
#include "oct-map.h"
#include "ov-base.h"
#include "ov-base-mat.h"
#include "ov-base-scalar.h"
#include "pr-output.h"

template <typename MT>
octave_value
octave_base_matrix<MT>::subsref (const std::string& type,
                                 const std::list<octave_value_list>& idx)
{
  octave_value retval;

  switch (type[0])
    {
    case '(':
      retval = do_index_op (idx.front ());
      break;

    case '{':
    case '.':
      {
        std::string nm = type_name ();
        error ("%s cannot be indexed with %c", nm.c_str (), type[0]);
      }
      break;

    default:
      panic_impossible ();
    }

  return retval.next_subsref (type, idx);
}

template <typename MT>
octave_value
octave_base_matrix<MT>::subsasgn (const std::string& type,
                                  const std::list<octave_value_list>& idx,
                                  const octave_value& rhs)
{
  octave_value retval;

  switch (type[0])
    {
    case '(':
      {
        if (type.length () == 1)
          retval = numeric_assign (type, idx, rhs);
        else if (isempty ())
          {
            // Allow conversion of empty matrix to some other type in
            // cases like
            //
            //  x = []; x(i).f = rhs

            if (type[1] != '.')
              error ("invalid assignment expression");

            octave_value tmp = octave_value::empty_conv (type, rhs);

            retval = tmp.subsasgn (type, idx, rhs);
          }
        else
          {
            std::string nm = type_name ();
            error ("in indexed assignment of %s, last lhs index must be ()",
                   nm.c_str ());
          }
      }
      break;

    case '{':
    case '.':
      {
        if (! isempty ())
          {
            std::string nm = type_name ();
            error ("%s cannot be indexed with %c", nm.c_str (), type[0]);
          }

        octave_value tmp = octave_value::empty_conv (type, rhs);

        retval = tmp.subsasgn (type, idx, rhs);
      }
      break;

    default:
      panic_impossible ();
    }

  return retval;
}

template <typename MT>
octave_value
octave_base_matrix<MT>::do_index_op (const octave_value_list& idx,
                                     bool resize_ok)
{
  octave_value retval;

  octave_idx_type n_idx = idx.length ();

  int nd = m_matrix.ndims ();
  const MT& cmatrix = m_matrix;

  // If we catch an indexing error in index_vector, we flag an error in
  // index k.  Ensure it is the right value before each idx_vector call.
  // Same variable as used in the for loop in the default case.

  octave_idx_type k = 0;

  try
    {
      switch (n_idx)
        {
        case 0:
          warn_empty_index (type_name ());
          retval = m_matrix;
          break;

        case 1:
          {
            octave::idx_vector i = idx (0).index_vector ();

            // optimize single scalar index.
            if (! resize_ok && i.is_scalar ())
              retval = cmatrix.checkelem (i(0));
            else
              retval = MT (m_matrix.index (i, resize_ok));
          }
          break;

        case 2:
          {
            octave::idx_vector i = idx (0).index_vector ();

            k=1;
            octave::idx_vector j = idx (1).index_vector ();

            // optimize two scalar indices.
            if (! resize_ok && i.is_scalar () && j.is_scalar ())
              retval = cmatrix.checkelem (i(0), j(0));
            else
              retval = MT (m_matrix.index (i, j, resize_ok));
          }
          break;

        default:
          {
            Array<octave::idx_vector> idx_vec (dim_vector (n_idx, 1));
            bool scalar_opt = n_idx == nd && ! resize_ok;
            const dim_vector dv = m_matrix.dims ();

            for (k = 0; k < n_idx; k++)
              {
                idx_vec(k) = idx(k).index_vector ();

                scalar_opt = (scalar_opt && idx_vec(k).is_scalar ());
              }

            if (scalar_opt)
              retval = cmatrix.checkelem (conv_to_int_array (idx_vec));
            else
              retval = MT (m_matrix.index (idx_vec, resize_ok));
          }
          break;
        }
    }
  catch (octave::index_exception& ie)
    {
      // Rethrow to allow more info to be reported later.
      ie.set_pos_if_unset (n_idx, k+1);
      throw;
    }

  return retval;
}

/*
%% This behavior is required for Matlab compatibility.
%!test
%! a = [1, 2; 3, 4];
%! assert (a(), a);
*/

template <typename MT>
void
octave_base_matrix<MT>::assign (const octave_value_list& idx, const MT& rhs)
{
  octave_idx_type n_idx = idx.length ();

  // If we catch an indexing error in index_vector, we flag an error in
  // index k.  Ensure it is the right value before each idx_vector call.
  // Same variable as used in the for loop in the default case.

  octave_idx_type k = 0;

  try
    {
      switch (n_idx)
        {
        case 0:
          panic_impossible ();
          break;

        case 1:
          {
            octave::idx_vector i = idx (0).index_vector ();

            m_matrix.assign (i, rhs);
          }
          break;

        case 2:
          {
            octave::idx_vector i = idx (0).index_vector ();

            k = 1;
            octave::idx_vector j = idx (1).index_vector ();

            m_matrix.assign (i, j, rhs);
          }
          break;

        default:
          {
            Array<octave::idx_vector> idx_vec (dim_vector (n_idx, 1));

            for (k = 0; k < n_idx; k++)
              idx_vec(k) = idx(k).index_vector ();

            m_matrix.assign (idx_vec, rhs);
          }
          break;
        }
    }
  catch (octave::index_exception& ie)
    {
      // Rethrow to allow more info to be reported later.
      ie.set_pos_if_unset (n_idx, k+1);
      throw;
    }

  // Clear cache.
  clear_cached_info ();
}

template <typename MT>
MatrixType
octave_base_matrix<MT>::matrix_type (const MatrixType& typ) const
{
  delete m_typ;
  m_typ = new MatrixType (typ);
  return *m_typ;
}

template <typename MT>
void
octave_base_matrix<MT>::assign (const octave_value_list& idx,
                                typename MT::element_type rhs)
{
  octave_idx_type n_idx = idx.length ();

  int nd = m_matrix.ndims ();

  MT mrhs (dim_vector (1, 1), rhs);

  // If we catch an indexing error in index_vector, we flag an error in
  // index k.  Ensure it is the right value before each idx_vector call.
  // Same variable as used in the for loop in the default case.

  octave_idx_type k = 0;

  try
    {
      switch (n_idx)
        {
        case 0:
          panic_impossible ();
          break;

        case 1:
          {
            octave::idx_vector i = idx (0).index_vector ();

            // optimize single scalar index.
            if (i.is_scalar () && i(0) < m_matrix.numel ())
              m_matrix(i(0)) = rhs;
            else
              m_matrix.assign (i, mrhs);
          }
          break;

        case 2:
          {
            octave::idx_vector i = idx (0).index_vector ();

            k = 1;
            octave::idx_vector j = idx (1).index_vector ();

            // optimize two scalar indices.
            if (i.is_scalar () && j.is_scalar () && nd == 2
                && i(0) < m_matrix.rows () && j(0) < m_matrix.columns ())
              m_matrix(i(0), j(0)) = rhs;
            else
              m_matrix.assign (i, j, mrhs);
          }
          break;

        default:
          {
            Array<octave::idx_vector> idx_vec (dim_vector (n_idx, 1));
            bool scalar_opt = n_idx == nd;
            const dim_vector dv = m_matrix.dims ().redim (n_idx);

            for (k = 0; k < n_idx; k++)
              {
                idx_vec(k) = idx(k).index_vector ();

                scalar_opt = (scalar_opt && idx_vec(k).is_scalar ()
                              && idx_vec(k)(0) < dv(k));
              }

            if (scalar_opt)
              {
                // optimize all scalar indices.  Don't construct
                // an index array, but rather calc a scalar index directly.
                octave_idx_type n = 1;
                octave_idx_type j = 0;
                for (octave_idx_type i = 0; i < n_idx; i++)
                  {
                    j += idx_vec(i)(0) * n;
                    n *= dv (i);
                  }
                m_matrix(j) = rhs;
              }
            else
              m_matrix.assign (idx_vec, mrhs);
          }
          break;
        }
    }
  catch (octave::index_exception& ie)
    {
      // Rethrow to allow more info to be reported later.
      ie.set_pos_if_unset (n_idx, k+1);
      throw;
    }

  // Clear cache.
  clear_cached_info ();
}

template <typename MT>
void
octave_base_matrix<MT>::delete_elements (const octave_value_list& idx)
{
  octave_idx_type len = idx.length ();

  Array<octave::idx_vector> ra_idx (dim_vector (len, 1));

  for (octave_idx_type i = 0; i < len; i++)
    ra_idx(i) = idx(i).index_vector ();

  m_matrix.delete_elements (ra_idx);

  // Clear cache.
  clear_cached_info ();
}

template <typename MT>
octave_value
octave_base_matrix<MT>::resize (const dim_vector& dv, bool fill) const
{
  MT retval (m_matrix);
  if (fill)
    retval.resize (dv, 0);
  else
    retval.resize (dv);
  return retval;
}

// Return true if this matrix has all true elements (non-zero, not NA/NaN).
template <typename MT>
bool
octave_base_matrix<MT>::is_true (void) const
{
  bool retval = false;
  dim_vector dv = m_matrix.dims ();
  int nel = dv.numel ();

  if (nel > 0)
    {
      MT t1 (m_matrix.reshape (dim_vector (nel, 1)));

      if (t1.any_element_is_nan ())
        octave::err_nan_to_logical_conversion ();

      if (nel > 1)
        warn_array_as_logical (dv);

      boolNDArray t2 = t1.all ();

      retval = t2(0);
    }

  return retval;
}

template <typename MT>
bool
octave_base_matrix<MT>::print_as_scalar (void) const
{
  dim_vector dv = dims ();

  return (dv.all_ones () || dv.any_zero ());
}

template <typename MT>
void
octave_base_matrix<MT>::print (std::ostream& os, bool pr_as_read_syntax)
{
  print_raw (os, pr_as_read_syntax);
  newline (os);
}

template <typename MT>
void
octave_base_matrix<MT>::print_info (std::ostream& os,
                                    const std::string& prefix) const
{
  m_matrix.print_info (os, prefix);
}

template <typename MT>
void
octave_base_matrix<MT>::short_disp (std::ostream& os) const
{
  if (m_matrix.isempty ())
    os << "[]";
  else if (m_matrix.ndims () == 2)
    {
      // FIXME: should this be configurable?
      octave_idx_type max_elts = 10;
      octave_idx_type elts = 0;

      octave_idx_type nel = m_matrix.numel ();

      octave_idx_type nr = m_matrix.rows ();
      octave_idx_type nc = m_matrix.columns ();

      os << '[';

      for (octave_idx_type i = 0; i < nr; i++)
        {
          for (octave_idx_type j = 0; j < nc; j++)
            {
              std::ostringstream buf;
              octave_print_internal (buf, m_matrix(j*nr+i));
              std::string tmp = buf.str ();
              std::size_t pos = tmp.find_first_not_of (' ');
              if (pos != std::string::npos)
                os << tmp.substr (pos);
              else if (! tmp.empty ())
                os << tmp[0];

              if (++elts >= max_elts)
                goto done;

              if (j < nc - 1)
                os << ", ";
            }

          if (i < nr - 1 && elts < max_elts)
            os << "; ";
        }

    done:

      if (nel <= max_elts)
        os << ']';
    }
  else
    os << "...";
}

template <typename MT>
float_display_format
octave_base_matrix<MT>::get_edit_display_format (void) const
{
  return make_format (m_matrix);
}

template <typename MT>
std::string
octave_base_matrix<MT>::edit_display (const float_display_format& fmt,
                                      octave_idx_type i,
                                      octave_idx_type j) const
{
  std::ostringstream buf;
  octave_print_internal (buf, fmt, m_matrix(i, j));
  return buf.str ();
}

template <typename MT>
octave_value
octave_base_matrix<MT>::fast_elem_extract (octave_idx_type n) const
{
  if (n < m_matrix.numel ())
    return m_matrix(n);
  else
    return octave_value ();
}

template <typename MT>
bool
octave_base_matrix<MT>::fast_elem_insert (octave_idx_type n,
                                          const octave_value& x)
{
  if (n < m_matrix.numel ())
    {
      // Don't use builtin_type () here to avoid an extra VM call.
      typedef typename MT::element_type ET;
      const builtin_type_t btyp = class_to_btyp<ET>::btyp;
      if (btyp == btyp_unknown) // Dead branch?
        return false;

      // Set up the pointer to the proper place.
      void *here = reinterpret_cast<void *> (&m_matrix(n));
      // Ask x to store there if it can.
      return x.get_rep ().fast_elem_insert_self (here, btyp);
    }
  else
    return false;
}
