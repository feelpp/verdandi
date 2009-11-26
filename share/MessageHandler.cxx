// Copyright (C) 2009, INRIA
// Author(s): Claire Mouton
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


#ifndef VERDANDI_FILE_MESSAGEHANDLER_CXX


#include "MessageHandler.hxx"


namespace Verdandi
{


    /////////////////////
    // MESSAGE HANDLER //
    /////////////////////


    //! Adds a new object in the recipient list.
    /*!
      \tparam R type of the recipient object.
      \param[in] recipient the string describing the object to add.
      \param[in] object reference to the recipient object.
      \param[in] pointer the pointer to the method to add.
    */
    template <class R>
    void MessageHandler
    ::AddRecipient(string recipient, R& object,
                   MessageHandler::function_pointer pointer)
    {
        pair<void*, function_pointer>
            pointer_pair(reinterpret_cast<void*>(&object), pointer);
        recipient_map_[recipient].push_back(pointer_pair);
        recipient_map_["all"].push_back(pointer_pair);
    }


    //! Adds a new object in the recipient list.
    /*!
      \param[in] recipient the string describing the object to add.
      \param[in] object pointer to the recipient object.
      \param[in] pointer the pointer to the method to add.
    */
    void MessageHandler
    ::AddRecipient(string recipient, void* object,
                   MessageHandler::function_pointer pointer)
    {
        pair<void*, function_pointer> pointer_pair(object, pointer);
        recipient_map_[recipient].push_back(pointer_pair);
        recipient_map_["all"].push_back(pointer_pair);
    }


    //! Sends a message to a recipient.
    /*!
      \param[in] recipient the recipient of the message.
      \param[in] message the string containing the message.
    */
    void MessageHandler
    ::Send(string recipient, string message)
    {
#ifndef VERDANDI_IGNORE_MESSAGE
        if (recipient_map_.count(recipient) == 0)
            throw ErrorArgument("MessageHandler::Send",
                                string("The object \"") + recipient +
                                "\" is not part of the recipient list.");

        recipient_list my_list = recipient_map_[recipient];
        SendToList(my_list, message);
#endif
    }


    //! Sends a message to a recipient and specifies the sender.
    /*!
      \param[in] sender the object sending the message.
      \param[in] recipient the recipient of the message.
      \param[in] message the string containing the message.
      \tparam Sender the type of the sending object.
    */
    template <class Sender>
    void MessageHandler
    ::Send(const Sender& sender, string recipient, string message)
    {
#ifndef VERDANDI_IGNORE_MESSAGE
        if (recipient_map_.count(recipient) == 0)
            throw ErrorArgument("MessageHandler::Send",
                                string("The object \"") + recipient +
                                "\" is not part of the recipient list.");

        recipient_list my_list = recipient_map_[recipient];
        SendToList(my_list, "Message from " + sender.GetName() + ": "
                   + message);
#endif
    }


    //! Sends a message to a list of recipients.
    /*!
      \param[in] my_list the list of recipients.
      \param[in] message the string containing the message.
    */
    void MessageHandler::SendToList(recipient_list& my_list, string message)
    {
#ifndef VERDANDI_IGNORE_MESSAGE
        recipient_list::iterator my_iterator;
        for(my_iterator = my_list.begin(); my_iterator != my_list.end();
            ++my_iterator)
            (*my_iterator->second)(my_iterator->first, message);
#endif
    }


    MessageHandler::recipient_map MessageHandler::recipient_map_;


} // namespace Verdandi.


#define VERDANDI_FILE_MESSAGEHANDLER_CXX
#endif