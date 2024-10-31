        // Check if all values in convergence list are true, then stop and write
        if (std::all_of(
                std::begin(convergedList), 
                std::end(convergedList), 
                [](bool i)
                    {
                        return i;
                    }
        ))
        {
            Info << "Simulation has reached convergence criterion!" << endl;
            runTime.writeAndEnd();
        }