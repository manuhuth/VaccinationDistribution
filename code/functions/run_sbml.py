import amici


def get_model_and_solver_from_sbml(path_sbml, model_name, model_directory):
    """Get model and solver objects from a SBML file.

     Parameters
     ----------
     path_sbml : str
         Path to sbml file that contains the model.

     model_name : str
         Name of the model.

     model_directory : str
         Path and name of where model directory should be stored.

     Returns
     -------
    model_solver : dict
        Dictionary containing model and solver as keys and as values
        the respective objects from the sbml file.

    """

    filename = path_sbml + ".xml"
    sbml_importer = amici.SbmlImporter(filename)
    sbml_importer.sbml2amici(model_name, model_directory)

    model_module = amici.import_model_module(model_name, model_directory)
    model = model_module.getModel()
    solver = model.getSolver()

    model_solver = {"model": model, "solver": solver}

    return model_solver


def model_run(model, solver, timepoints, set_parameter=None):
    """Run model and return outputs for specification.

    Parameters
    ----------
    model : amici.Model
        Model to run.

    solver : amici.Solver
        Solver for amici model.

    timepoints : array
        Defining the grid for which model should be computed.

    set_parameter : dict
        Allows to change parameter values by passing a dictionary. Keys must
        be the names of the parameters and the values the magnitudes of the
        parameters.

    Returns
    -------
    rdata : numpy.ReturnDataView
        Output data.

    """

    if set_parameter is not None:
        for keys in set_parameter.keys():
            model.setParameterByName(keys, set_parameter[keys])

    model.setTimepoints(timepoints)
    rdata = amici.runAmiciSimulation(model, solver)

    return rdata
