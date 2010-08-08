// Copyright (C) 2010, INRIA
// Author(s): Vivien Mallet
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.


#ifndef VERDANDI_FILE_SHARE_OPS_CXX

#include "ops/Ops.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    /*! Nothing is performed. The Lua state is set to NULL.
     */
    Ops::Ops(): ::Ops::Ops()
    {
    }


    //! Main constructor.
    /*! The Lua configuration file is loaded and run. An exception may be
      raised during this evaluation.
      \param[in] file_path path to the configuration file.
    */
    Ops::Ops(string file_path): ::Ops::Ops(file_path)
    {
    }


    //! Destructor.
    /*! Destroys the Lua state object.
     */
    Ops::~Ops()
    {
    }


    //////////////////
    // MAIN METHODS //
    //////////////////


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[in] constraint constraint that the entry value must satisfy.
      \param[in] default_value default value for the entry in case it is not
      found in the configuration file.
      \param[out] value value of the entry.
    */
    template<class TD, class T>
    void
    Ops::Set(string name, string constraint,
             const TD& default_value, T& value)
    {
        SetValue(name, constraint, default_value, true, value);
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[in] constraint constraint that the entry value must satisfy.
      \param[out] value value of the entry.
    */
    template<class T>
    void Ops::Set(string name, string constraint, T& value)
    {
        SetValue(name, constraint, value, false, value);
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[out] value value of the entry.
    */
    template <class T>
    void Ops::Set(string name, T& value)
    {
        SetValue(name, "", value, false, value);
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[in] constraint constraint that the entry value must satisfy.
      \param[in] default_value default value for the entry in case it is not
      found in the configuration file.
      \return The value of the entry.
    */
    template<class T>
    T Ops::Get(string name, string constraint, const T& default_value)
    {
        T value;
        SetValue(name, constraint, default_value, true, value);
        return value;
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \param[in] constraint constraint that the entry value must satisfy.
      \return The value of the entry.
    */
    template<class T>
    T Ops::Get(string name, string constraint)
    {
        T value;
        SetValue(name, constraint, value, false, value);
        return value;
    }


    //! Retrieves a value from the configuration file.
    /*!
      \param[in] name name of the entry.
      \return The value of the entry.
    */
    template <class T>
    T Ops::Get(string name)
    {
        T value;
        SetValue(name, "", value, false, value);
        return value;
    }


    //! Checks whether \a name is of type 'T'.
    /*! On exit, the value of the entry (if it exists) is on the stack.
      \param[in] name the name of the entry whose type is checked.
      \return True if the entry is of type 'T', false otherwise.
      \note The prefix is prepended to \a name. If \a name does not exist, an
      exception is raised.
    */
    template<class T>
    bool Ops::Is(string name)
    {
        T value;
        return IsParam(name, value);
    }


    //! Retrieves a value and checks if it satisfies given constraints.
    /*! If the entry is not found, the default value is returned (if any). If
      the vector is given in a file, the Lua entry should provide the file
      name. The file is read with method 'Read' of \a value.
      \param[in] name name of the entry.
      \param[in] constraint constraint to be satisfied.
      \param[in] default_value default value.
      \param[in] with_default is there a default value? If not, \a
      default_value is ignored.
      \param[out] value the value of the entry named \a name.
      \note The default value may not satisfy the constraint.
    */
    template<class T, class Allocator>
    void Ops::SetValue(string name, string constraint,
                       const Seldon::Vector<T, VectFull, Allocator>&
                       default_value,
                       bool with_default,
                       Seldon::Vector<T, VectFull, Allocator>& value)
    {
        if (!this->Exists(name))
        {
            if (with_default)
                value = default_value;
            else
                throw Error("SetValue",
                            "The " + Entry(name) + " was not found.");
        }
        else if (this->Is<string>(name))
        {
            string filename;
            SetValue(name, "", "", false, filename);
            value.Read(filename);
            if (constraint.empty())
                return;
            for (int i = 0; i < value.GetLength(); i++)
                if (!CheckConstraintOnValue(to_str(value(i)), constraint))
                    throw Error("SetValue",
                                "The entry "
                                + Entry(name + "[" + to_str(i + 1) + "]")
                                + " does not satisfy the constraint:\n"
                                + Constraint(constraint));
        }
        else
        {
            std::vector<T> default_data(default_value.GetLength()), data;
            for (int i = 0; i < default_value.GetLength(); i++)
                default_data[i] = default_value(i);
            SetValue(name, constraint, default_data, with_default, data);
            value.Reallocate(int(data.size()));
            for (size_t i = 0; i < data.size(); i++)
                value(i) = data[i];
        }
    }


    //! Checks whether \a name is a table of 'T'.
    /*!
      \param[in] name the name of the entry whose type is checked.
      \param[in] value anything: it is used to determine the type.
      \return True if the entry is a 'Vector<T>', false otherwise.
      \note The prefix is prepended to \a name. If \a name does not exist, an
      exception is raised.
    */
    template<class T, class Allocator>
    bool Ops::IsParam(string name,
                      Seldon::Vector<T, VectFull, Allocator>& value)
    {
        string str_value;
        std::vector<T> vect_value;
        return this->IsParam(name, str_value)
            || this->IsParam(name, vect_value);
    }


} // namespace Verdandi.


#define VERDANDI_FILE_SHARE_OPS_CXX
#endif

