import pandas as pd
import matplotlib.pyplot as plt

def get_state_trajectory_data_frame(results, model):
    states = model.getStateIds()
    trajectories = results['x']
    df_trajectories = pd.DataFrame(trajectories, columns=(states)) 
    return df_trajectories

def plot_states(results, model, xlabel='$t$ (d)', ylabel='$x_i(t)$', title='State trajectories', state_ids=None):
    fig, ax = plt.subplots()
    df_trajectories = get_state_trajectory_data_frame(results, model)
    
    if state_ids is None:
        states = df_trajectories.columns
    else:
        states = state_ids
    
    for index in states:
            label = index
            ax.plot(results['t'], df_trajectories[index], label = label)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.set_title(title)
            
    return fig, ax

def get_substates(model, substrings):
    states = model.getStateIds()
    substates = get_states_by_substrings(states=states, substrings=substrings)
    return substates

def get_states_by_substrings(states, substrings):
    states_removed = [
        x for x in substrings if x in states
    ]
    return states_removed

    