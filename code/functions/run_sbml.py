import amici

def get_model_and_solver_from_sbml(path_sbml, model_name, model_directory):
#TODO incorporate optimizer options
    #path = 'test.xml'
    #name = 'test'
    #directory = 'model_dir'
    filename = path_sbml + '.xml'

    sbml_importer = amici.SbmlImporter(filename)

    sbml_importer.sbml2amici(model_name, model_directory)
    
    # load the model module
    model_module = amici.import_model_module(model_name, model_directory)
    # instantiate model
    model = model_module.getModel()
    # instantiate solver
    solver = model.getSolver()
    
    return {'model' : model, 'solver' : solver}


def model_run(model, solver, timepoints, set_parameter=None):
   
    
    if set_parameter is not None:
        for keys in set_parameter.keys():
            model.setParameterByName(keys , set_parameter[keys])
    
    model.setTimepoints(timepoints)
    rdata = amici.runAmiciSimulation(model, solver)
    
    return rdata