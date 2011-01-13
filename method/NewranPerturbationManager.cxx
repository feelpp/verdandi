// Copyright (C) 2010 INRIA
// Author(s): Vivien Mallet, Anne Tilloy
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
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_METHOD_NEWRANPERTURBATIONMANAGER_CXX


#include "NewranPerturbationManager.hxx"
#include "BasePerturbationManager.cxx"
#include "share/LockFile.cxx"

namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    NewranPerturbationManager
    ::NewranPerturbationManager(): BasePerturbationManager(), urng_(NULL)
    {
    }


    //! Main constructor.
    /*! Builds the manager and reads option keys in the configuration file.
      \param[in] configuration_file configuration file.
    */
    NewranPerturbationManager
    ::NewranPerturbationManager(string configuration_file):
        BasePerturbationManager(), urng_(NULL)
    {
        Initialize(configuration_file);
    }


    //! Destructor.
    NewranPerturbationManager::~NewranPerturbationManager()
    {
        if (urng_ != NULL)
            delete urng_;
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the manager.
    /*!
      \param[in] configuration_file configuration file.
    */
    void NewranPerturbationManager::Initialize(string configuration_file)
    {
        Ops configuration_stream(configuration_file);

        configuration_stream.SetPrefix("perturbation_manager.newran.");

        configuration_stream.Set("seed_path", seed_path_);
        NEWRAN::Random::SetDirectory(seed_path_.c_str());
        urng_ = new NEWRAN::LGM_mixed;
        NEWRAN::Random::Set(*urng_);
        if (!Lock(seed_path_ + "lock"))
            throw ErrorIO("NewranPerturbationManager::Initialize(string)",
                          "Unable to create the Newran lock file \""
                          + seed_path_ + "lock\".");
        NEWRAN::Random::CopySeedFromDisk();
    }


    //! Reinitializes the manager.
    /*! Locks and reloads the seed. */
    void NewranPerturbationManager::Reinitialize()
    {
        Finalize();

        if (!Lock(seed_path_ + "lock"))
            throw ErrorIO("NewranPerturbationManager::Reinitialize()",
                          "Unable to create the Newran lock file \""
                          + seed_path_ + "lock\".");

        NEWRAN::Random::CopySeedFromDisk();
    }


    //! Copies the seed to disk.
    /*! Saves and unlocks the seed. */
    void NewranPerturbationManager::Finalize()
    {
        NEWRAN::Random::CopySeedToDisk();

        if (!Unlock(seed_path_ + "lock"))
            throw ErrorIO("NewranPerturbationManager::Finalize()",
                          "Unable to remove the Newran lock file \""
                          + seed_path_ + "lock\".");
    }


    //! Generates a random number with normal distribution.
    /*!
      \param[in] mean mean of the normal distribution.
      \param[in] variance variance of the random variable.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal distribution, any random value lies in
      \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$ is the mean
      of the random variable and \f$ \sigma \f$ is its standard deviation.
      \return A random number following the previously described normal
      distribution.
    */
    double NewranPerturbationManager
    ::Normal(double mean, double variance,
             Vector<double, VectFull>& parameter)
    {
        NEWRAN::Normal N;
        double value = N.Next();
        if (parameter.GetLength() == 2)
            while (value < parameter(0) || value > parameter(1))
                value = N.Next();
        else if (parameter.GetLength() != 0)
            throw ErrorArgument("NewranPerturbationManager"
                                "::Normal(double, double, Vector)",
                                "The vector of parameters should be either "
                                "empty or of length 2, but it contains "
                                + to_str(parameter.GetLength())
                                + " element(s).");
        return mean + sqrt(variance) * value;
    }


    //! Generates a vector of random numbers with normal distribution.
    /*! Each component of the random vector is generated independently.
      \param[in] mean mean of the normal distribution.
      \param[in] variance variance of the random variable.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal distribution, any random value lies in
      \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$ is the mean
      of the random variable and \f$ \sigma \f$ is its standard deviation.
      \param[out] output the generated random vector.
    */
    void NewranPerturbationManager
    ::Normal(double mean, double variance,
             Vector<double, VectFull>& parameter, Vector<double>& sample)
    {
        if (parameter.GetLength() != 0 && parameter.GetLength() != 2)
            throw ErrorArgument("NewranPerturbationManager"
                                "::Normal(double, double, Vector, Vector)",
                                "The vector of parameters should be either "
                                "empty or of length 2, but it contains "
                                + to_str(parameter.GetLength())
                                + " element(s).");

        NEWRAN::Normal N;
        double value;
        int size = sample.GetSize();
        for (int i = 0; i < size; i++)
        {
            value = N.Next();
            if (parameter.GetLength() == 2)
                while (value < parameter(0) || value > parameter(1))
                    value = N.Next();
            sample(i) = mean + sqrt(variance) * value;
        }
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_NEWRANPERTURBATIONMANAGER_CXX
#endif
