const bool debugging(false);

// Read probesControl
dictionary probesControlDict = mesh.solutionDict().subDict("probesControl");

bool doProbesControl(false);
probesControlDict.readIfPresent("doProbesControl", doProbesControl);

if (doProbesControl)
{
    if (debugging)
    {
        Info<< "Probes control on! " << endl;
    }

    // Reading locationList
    const List<vector> locationList
    (
        probesControlDict.lookup("probeLocations")
    );

    if (debugging)
    {
        Info<< "Probes location: " << locationList << endl;
    }

    // Iterate over locationList
    forAll(locationList, locationi)
    {
        // Get location and cellID
        const vector location = locationList[locationi];
        const label celli = mesh.findCell(location);

        // do it only on the processor that HAS this point
        if (celli >= 0) 
        {
            Info<< "Probe " << locationi << ": "
                << "p = "   << p[celli] << " Pa" 
                << endl;
        }

        if (debugging)
        {
            Info<< "Probes cellID: " << celli << endl;
        }
    }
}
else
{
    if (debugging)
    {
        Info<< "Probes control off!.." << endl;
    };
}