import amici
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sbml_importer = amici.SbmlImporter('test.xml')

model_name = 'test'
model_dir = 'model_dir'
sbml_importer.sbml2amici(model_name, model_dir)

# load the model module
model_module = amici.import_model_module(model_name, model_dir)
# instantiate model
model = model_module.getModel()
# instantiate solver
solver = model.getSolver()


model.setParameterByName('beta' , 1)
model.setParameterByName('gamma', 0.13)


# set timepoints
model.setTimepoints(np.linspace(0, 100, 100))
rdata = amici.runAmiciSimulation(model, solver)

rdata.x

cols = ['Susceptible', 'Infectious', 'Recovered']
df = pd.DataFrame(rdata.x, columns = cols)

for compart in cols:
    plt.plot(df[compart], label=compart)
    
plt.legend()
    






