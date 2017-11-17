/*

Copyright (C) 1993-2017 John W. Eaton
Copyright (C) 2009 VZLU Prague

This file is part of Octave.

Octave is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Octave is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Octave; see the file COPYING.  If not, see
<http://www.gnu.org/licenses/>.

*/

#if ! defined (octave_symtab_h)
#define octave_symtab_h 1

#include "octave-config.h"

#include <deque>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <string>

#include "glob-match.h"
#include "lo-regexp.h"
#include "oct-refcount.h"

class tree_argument_list;
class octave_user_function;

#include "fcn-info.h"
#include "ov.h"
#include "ovl.h"
#include "symrec.h"
#include "symscope.h"
#include "workspace-element.h"

namespace octave
{
  class OCTINTERP_API symbol_table
  {
  public:

    typedef octave::symbol_record symbol_record;
    typedef octave::symbol_scope scope;
    typedef octave::fcn_info fcn_info;

    typedef scope::context_id context_id;

    symbol_table (void)
      : m_global_symbols (), m_fcn_table (), m_class_precedence_table (),
        m_parent_map (), m_global_scope (new scope ("global scope")),
        m_top_scope (new symbol_scope ("top scope")),
        m_current_scope (m_top_scope)
      { }

    // No copying!

    symbol_table (const symbol_table&) = delete;

    symbol_table& operator = (const symbol_table&) = delete;

    ~symbol_table (void)
      {
        delete m_top_scope;
        delete m_global_scope;
      }

    symbol_scope *global_scope (void) { return m_global_scope; }
    symbol_scope *top_scope (void) { return m_top_scope; }

    symbol_scope *current_scope (void) { return m_current_scope; }

    symbol_scope *require_current_scope (const std::string& who)
    {
      if (! m_current_scope)
        error ("%s: missing scope", who.c_str ());

      return m_current_scope;
    }

    context_id current_context (void) const
    {
      return m_current_scope ? m_current_scope->current_context () : 0;
    }

    void set_scope (symbol_scope *sid)
    {
      set_scope_and_context (sid, 0);
    }

    void set_scope_and_context (symbol_scope *sid, context_id context)
    {
      if (sid == m_global_scope)
        error ("can't set scope to global");

      m_current_scope = sid;

      if (m_current_scope)
        m_current_scope->set_context (context);
    }

    symbol_record find_symbol (const std::string& name, symbol_scope *sid)
    {
      return sid ? sid->find_symbol (name) : symbol_record ();
    }

    symbol_record find_symbol (const std::string& name)
    {
      return find_symbol (name, m_current_scope);
    }

    void inherit (symbol_scope *recipient_scope, symbol_scope *donor_scope)
    {
      if (recipient_scope)
        recipient_scope->inherit (donor_scope);
    }

    void inherit (symbol_scope *recipient_scope)
    {
      inherit (recipient_scope, m_current_scope);
    }

    bool at_top_level (void) { return m_current_scope == m_top_scope; }

    // Find a value corresponding to the given name in the table.
    octave_value
      find (const std::string& name,
            const octave_value_list& args = octave_value_list (),
            bool skip_variables = false,
            bool local_funcs = true);

    void assign (const std::string& name, const octave_value& value, bool force_add)
    {
      if (m_current_scope)
        m_current_scope->assign (name, value, force_add);
    }

    void assign (const std::string& name,
                 const octave_value& value = octave_value ())
    {
      if (m_current_scope)
        m_current_scope->assign (name, value);
    }

    octave_value varval (const std::string& name) const
    {
      return (m_current_scope
              ? m_current_scope->varval (name) : octave_value ());
    }

    void
      global_assign (const std::string& name,
                     const octave_value& value = octave_value ())
    {
      global_symbols_iterator p = m_global_symbols.find (name);

      if (p == m_global_symbols.end ())
        m_global_symbols[name] = value;
      else
        p->second = value;
    }

    octave_value& global_varref (const std::string& name)
    {
      global_symbols_iterator p = m_global_symbols.find (name);

      return (p == m_global_symbols.end ()
              ? m_global_symbols[name] : p->second);
    }

    octave_value global_varval (const std::string& name) const
    {
      global_symbols_const_iterator p = m_global_symbols.find (name);

      return (p != m_global_symbols.end ()) ? p->second : octave_value ();
    }

    void
      top_level_assign (const std::string& name,
                        const octave_value& value = octave_value ())
    {
      m_top_scope->assign (name, value);
    }

    octave_value top_level_varval (const std::string& name) const
    {
      return m_top_scope->varval (name);
    }

    bool
      is_built_in_function_name (const std::string& name)
    {
      octave_value val = find_built_in_function (name);

      return val.is_defined ();
    }

    octave_value
      find_method (const std::string& name, const std::string& dispatch_type)
    {
      fcn_table_const_iterator p = m_fcn_table.find (name);

      if (p != m_fcn_table.end ())
        {
          octave_value fcn = p->second.find_method (dispatch_type);

          if (! fcn.is_defined ())
            fcn = find_submethod (name, dispatch_type);

          return fcn;
        }
      else
        {
          fcn_info finfo (name);

          octave_value fcn = finfo.find_method (dispatch_type);

          if (! fcn.is_defined ())
            fcn = find_submethod (name, dispatch_type);

          if (fcn.is_defined ())
            m_fcn_table[name] = finfo;

          return fcn;
        }
    }

    octave_value
      find_submethod (const std::string& name, const std::string& dispatch_type);

    octave_value
      find_built_in_function (const std::string& name)
    {
      fcn_table_const_iterator p = m_fcn_table.find (name);

      return (p != m_fcn_table.end ()
              ? p->second.find_built_in_function () : octave_value ());
    }

    octave_value
      find_autoload (const std::string& name)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      return (p != m_fcn_table.end ()
              ? p->second.find_autoload () : octave_value ());
    }

    octave_value builtin_find (const std::string& name);

    octave_value
    fcn_table_find (const std::string& name,
                    const octave_value_list& args = octave_value_list (),
                    bool local_funcs = true);

    octave_value
      find_function (const std::string& name,
                     const octave_value_list& args = octave_value_list (),
                     bool local_funcs = true);

    octave_value find_user_function (const std::string& name)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      return (p != m_fcn_table.end ()
              ? p->second.find_user_function () : octave_value ());
    }

    octave_value find_cmdline_function (const std::string& name)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      return (p != m_fcn_table.end ()
              ? p->second.find_cmdline_function () : octave_value ());
    }

    void install_cmdline_function (const std::string& name,
                                   const octave_value& fcn)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      if (p != m_fcn_table.end ())
        {
          fcn_info& finfo = p->second;

          finfo.install_cmdline_function (fcn);
        }
      else
        {
          fcn_info finfo (name);

          finfo.install_cmdline_function (fcn);

          m_fcn_table[name] = finfo;
        }
    }

    // Install local function FCN named NAME.  FILE_NAME is the name of
    // the file containing the local function.

    void install_local_function (const std::string& name,
                                 const octave_value& fcn,
                                 const std::string& file_name)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      if (p != m_fcn_table.end ())
        {
          fcn_info& finfo = p->second;

          finfo.install_local_function (fcn, file_name);
        }
      else
        {
          fcn_info finfo (name);

          finfo.install_local_function (fcn, file_name);

          m_fcn_table[name] = finfo;
        }
    }

    void install_user_function (const std::string& name,
                                const octave_value& fcn)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      if (p != m_fcn_table.end ())
        {
          fcn_info& finfo = p->second;

          finfo.install_user_function (fcn);
        }
      else
        {
          fcn_info finfo (name);

          finfo.install_user_function (fcn);

          m_fcn_table[name] = finfo;
        }
    }

    void install_built_in_function (const std::string& name,
                                    const octave_value& fcn)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      if (p != m_fcn_table.end ())
        {
          fcn_info& finfo = p->second;

          finfo.install_built_in_function (fcn);
        }
      else
        {
          fcn_info finfo (name);

          finfo.install_built_in_function (fcn);

          m_fcn_table[name] = finfo;
        }
    }

    void clear_all (bool force = false)
    {
      if (m_current_scope)
        {
          m_current_scope->clear_variables ();
          m_current_scope->clear_global_pattern ("*");
        }

      clear_functions (force);
    }

    // This is written as two separate functions instead of a single
    // function with default values so that it will work properly with
    // unwind_protect.

    void clear_functions (bool force = false)
    {
      fcn_table_iterator p = m_fcn_table.begin ();

      while (p != m_fcn_table.end ())
        (p++)->second.clear (force);
    }

    void clear_function (const std::string& name)
    {
      clear_user_function (name);
    }

    void clear_symbol (const std::string& name)
    {
      // FIXME: are we supposed to do both here?

      if (m_current_scope)
        m_current_scope->clear_variable (name);

      clear_function (name);
    }

    void clear_function_pattern (const std::string& pat)
    {
      glob_match pattern (pat);

      fcn_table_iterator p = m_fcn_table.begin ();

      while (p != m_fcn_table.end ())
        {
          if (pattern.match (p->first))
            (p++)->second.clear_user_function ();
          else
            p++;
        }
    }

    void clear_symbol_pattern (const std::string& pat)
    {
      // FIXME: are we supposed to do both here?

      if (m_current_scope)
        m_current_scope->clear_variable_pattern (pat);

      clear_function_pattern (pat);
    }

    void clear_user_function (const std::string& name)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      if (p != m_fcn_table.end ())
        {
          fcn_info& finfo = p->second;

          finfo.clear_user_function ();
        }
      // FIXME: is this necessary, or even useful?
      // else
      //   error ("clear: no such function '%s'", name.c_str ());
    }

    // This clears oct and mex files, including autoloads.
    void clear_dld_function (const std::string& name)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      if (p != m_fcn_table.end ())
        {
          fcn_info& finfo = p->second;

          finfo.clear_autoload_function ();
          finfo.clear_user_function ();
        }
    }

    void clear_mex_functions (void)
    {
      fcn_table_iterator p = m_fcn_table.begin ();

      while (p != m_fcn_table.end ())
        (p++)->second.clear_mex_function ();
    }

    void erase_global (const std::string& name);

    void erase_global_pattern (const glob_match& pattern);

    bool set_class_relationship (const std::string& sup_class,
                                 const std::string& inf_class);

    bool is_superiorto (const std::string& a, const std::string& b);

    void alias_built_in_function (const std::string& alias,
                                  const std::string& name)
    {
      octave_value fcn = find_built_in_function (name);

      if (fcn.is_defined ())
        {
          fcn_info finfo (alias);

          finfo.install_built_in_function (fcn);

          m_fcn_table[alias] = finfo;
        }
      else
        panic ("alias: '%s' is undefined", name.c_str ());
    }

    void install_built_in_dispatch (const std::string& name,
                                    const std::string& klass)
    {
      fcn_table_iterator p = m_fcn_table.find (name);

      if (p != m_fcn_table.end ())
        {
          fcn_info& finfo = p->second;

          finfo.install_built_in_dispatch (klass);
        }
      else
        error ("install_built_in_dispatch: '%s' is undefined", name.c_str ());
    }

    void push_context (void)
    {
      if (m_current_scope == m_top_scope)
        error ("invalid call to symtab::push_context");

      if (m_current_scope)
        m_current_scope->push_context ();
    }

    // This is written as two separate functions instead of a single
    // function with default values so that it will work properly with
    // unwind_protect.

    void pop_context (void)
    {
      if (m_current_scope == m_top_scope)
        error ("invalid call to symtab::pop_context");

      if (m_current_scope)
        m_current_scope->pop_context ();
    }

    // For unwind_protect where a pointer argument is needed.

    void pop_context (void *) { pop_context (); }

    std::list<symbol_record> glob (const std::string& pattern)
    {
      return (m_current_scope
              ? m_current_scope->glob (pattern) : std::list<symbol_record> ());
    }

    std::list<symbol_record>
      glob_global_variables (const std::string& pattern)
    {
      std::list<symbol_record> retval;

      glob_match pat (pattern);

      for (const auto& nm_val : m_global_symbols)
        {
          // We generate a list of symbol_record objects so that the results from
          // glob_variables and glob_global_variables may be handled the same
          // way.
          if (pat.match (nm_val.first))
            retval.push_back (symbol_record (m_global_scope,
                                             nm_val.first, nm_val.second,
                                             symbol_record::global));
        }

      return retval;
    }

    std::list<symbol_record>
      regexp_global_variables (const std::string& pattern)
    {
      std::list<symbol_record> retval;

      octave::regexp pat (pattern);

      for (const auto& nm_val : m_global_symbols)
        {
          // We generate a list of symbol_record objects so that the results from
          // regexp_variables and regexp_global_variables may be handled the same
          // way.
          if (pat.is_match (nm_val.first))
            retval.push_back (symbol_record (m_global_scope,
                                             nm_val.first, nm_val.second,
                                             symbol_record::global));
        }

      return retval;
    }

    std::list<symbol_record> glob_variables (const string_vector& patterns)
    {
      std::list<symbol_record> retval;

      if (! m_current_scope)
        return retval;

      size_t len = patterns.numel ();

      for (size_t i = 0; i < len; i++)
        {
          std::list<symbol_record> tmp = m_current_scope->glob (patterns[i]);

          retval.insert (retval.begin (), tmp.begin (), tmp.end ());
        }

      return retval;
    }

    std::list<symbol_record> regexp_variables (const string_vector& patterns)
    {
      std::list<symbol_record> retval;

      if (! m_current_scope)
        return retval;

      size_t len = patterns.numel ();

      for (size_t i = 0; i < len; i++)
        {
          std::list<symbol_record> tmp = m_current_scope->regexp (patterns[i]);

          retval.insert (retval.begin (), tmp.begin (), tmp.end ());
        }

      return retval;
    }

    std::list<std::string> user_function_names (void)
    {
      std::list<std::string> retval;

      for (const auto& nm_finfo : m_fcn_table)
        {
          if (nm_finfo.second.is_user_function_defined ())
            retval.push_back (nm_finfo.first);
        }

      if (! retval.empty ())
        retval.sort ();

      return retval;
    }

    std::list<std::string> global_variable_names (void)
    {
      std::list<std::string> retval;

      for (const auto& nm_val : m_global_symbols)
        retval.push_back (nm_val.first);

      retval.sort ();

      return retval;
    }

    std::list<std::string> top_level_variable_names (void)
    {
      return (m_top_scope
              ? m_top_scope->variable_names () : std::list<std::string> ());
    }

    std::list<std::string> variable_names (void)
    {
      return (m_current_scope
              ? m_current_scope->variable_names () : std::list<std::string> ());
    }

    std::list<std::string> built_in_function_names (void)
    {
      std::list<std::string> retval;

      for (const auto& nm_finfo : m_fcn_table)
        {
          octave_value fcn = nm_finfo.second.find_built_in_function ();

          if (fcn.is_defined ())
            retval.push_back (nm_finfo.first);
        }

      if (! retval.empty ())
        retval.sort ();

      return retval;
    }

    std::list<std::string> cmdline_function_names (void)
    {
      std::list<std::string> retval;

      for (const auto& nm_finfo : m_fcn_table)
        {
          octave_value fcn = nm_finfo.second.find_cmdline_function ();

          if (fcn.is_defined ())
            retval.push_back (nm_finfo.first);
        }

      if (! retval.empty ())
        retval.sort ();

      return retval;
    }

    octave_value dump (void) const;

    void add_to_parent_map (const std::string& classname,
                            const std::list<std::string>& parent_list)
    {
      m_parent_map[classname] = parent_list;
    }

    std::list<std::string>
      parent_classes (const std::string& dispatch_type)
    {
      std::list<std::string> retval;

      const_parent_map_iterator it = m_parent_map.find (dispatch_type);

      if (it != m_parent_map.end ())
        retval = it->second;

      for (const auto& nm : retval)
        {
          // Search for parents of parents and append them to the list.

          // FIXME: should we worry about a circular inheritance graph?

          std::list<std::string> parents = parent_classes (nm);

          if (! parents.empty ())
            retval.insert (retval.end (), parents.begin (), parents.end ());
        }

      return retval;
    }

    octave_user_function * get_curr_fcn (void)
    {
      return m_current_scope ? m_current_scope->function () : nullptr;
    }

    void cleanup (void);

    fcn_info * get_fcn_info (const std::string& name)
    {
      fcn_table_iterator p = m_fcn_table.find (name);
      return p != m_fcn_table.end () ? &p->second : nullptr;
    }

  private:

    typedef std::map<std::string, octave_value>::const_iterator
      global_symbols_const_iterator;
    typedef std::map<std::string, octave_value>::iterator
      global_symbols_iterator;

    typedef std::map<std::string, fcn_info>::const_iterator
      fcn_table_const_iterator;
    typedef std::map<std::string, fcn_info>::iterator
      fcn_table_iterator;

    // Map from names of global variables to values.
    std::map<std::string, octave_value> m_global_symbols;

    // Map from function names to function info (private
    // functions, class constructors, class methods, etc.)
    // Note that subfunctions are defined in the scope that contains
    // them.
    std::map<std::string, fcn_info> m_fcn_table;

    // Map from class names to set of classes that have lower
    // precedence.
    std::map<std::string, std::set<std::string>> m_class_precedence_table;

    typedef std::map<std::string, std::set<std::string>>::const_iterator
      class_precedence_table_const_iterator;
    typedef std::map<std::string, std::set<std::string>>::iterator
      class_precedence_table_iterator;

    // Map from class names to parent class names.
    std::map<std::string, std::list<std::string>> m_parent_map;

    typedef std::map<std::string, std::list<std::string>>::const_iterator
      const_parent_map_iterator;
    typedef std::map<std::string, std::list<std::string>>::iterator
      parent_map_iterator;

    symbol_scope *m_global_scope;
    symbol_scope *m_top_scope;

    symbol_scope *m_current_scope;

    octave_value dump_fcn_table_map (void) const;
  };

  extern bool out_of_date_check (octave_value& function,
                                 const std::string& dispatch_type = "",
                                 bool check_relative = true);

  extern OCTINTERP_API std::string
  get_dispatch_type (const octave_value_list& args);

  extern OCTINTERP_API std::string
  get_dispatch_type (const octave_value_list& args,
                     builtin_type_t& builtin_type);
}

#endif
