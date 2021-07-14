import pandas as pd
import numpy as np


def create_inflow_from_data(url="https://opendata.ecdc.europa.eu/covid19/vaccine_tracker/xlsx/data.xlsx",
                            column_interest = "NumberDosesReceived",
                            group_id = "ALL",
                            column_group = 'TargetGroup',
                            sheet_name = "data",
                            date_variable = "YearWeekISO",
                            country_variable = "ReportingCountry",
                            population_variable = "Population",
                            vaccine_variable = "Vaccine",
                            vaccine_allocation_dict = {"COM" : "vac1", "MOD" : "vac1", "AZ" : "vac2", "JANSS" : "vac2"}, #do not use unknown
                            number_decision_periods = 10,
                            weeks_decision_period = 2,
                            model_population = 160000000,
                            format_df = "long"):
    
    df_raw = pd.read_excel(url,
                           sheet_name=sheet_name,
                           skipfooter=0)  
    vaccines_to_use = vaccine_allocation_dict.keys()
    df_all = df_raw.loc[(df_raw[column_group] == group_id)*(df_raw[vaccine_variable].isin(vaccines_to_use))]
    
    column_to_keep = [date_variable, country_variable, column_interest, vaccine_variable]
    df_column_reduced = df_all[column_to_keep]
    df_column_reduced["vaccine_model"] = df_column_reduced[vaccine_variable].map(vaccine_allocation_dict)
    
    df_grouped = df_column_reduced.groupby([date_variable, "vaccine_model"]).sum()
    
    df_wide = df_grouped.pivot_table(index=[date_variable], 
                        columns='vaccine_model', 
                        values=column_interest)
    
    df_non_zero = df_wide[(df_wide != 0).all(1)].reset_index()
    
    length_df = df_non_zero.shape[0]
    
    if length_df >= number_decision_periods*2:
        df_periods = df_non_zero.iloc[range((number_decision_periods*2))]
    else:
        diff = number_decision_periods * 2 - length_df
        df_periods = df_non_zero.append( df_non_zero.iloc[[-1]*diff] ).reset_index()
    
    
    sequence = np.floor(np.array(range(0,number_decision_periods * 2)) / 2)
    df_periods["decision_period"] = sequence
    population_size = df_raw[population_variable].unique().sum()
    linear_pop_scaling = model_population / population_size
    
    df_supply_model = df_periods.groupby("decision_period").sum() / weeks_decision_period / 7 * linear_pop_scaling
    #df_supply_model = df_supply_model.drop(labels=["index"], axis = 1)
    if format_df == "long":
        df_supply_model = df_supply_model.unstack()

    return df_supply_model




