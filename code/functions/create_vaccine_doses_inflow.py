import pandas as pd
import numpy as np


def create_inflow_from_data(
    url="https://opendata.ecdc.europa.eu/covid19/vaccine_tracker/xlsx/data.xlsx",
    url_uk="https://coronavirus.data.gov.uk/api/v1/data?filters=areaType=overview&structure=%7B%22areaType%22:%22areaType%22,%22areaName%22:%22areaName%22,%22areaCode%22:%22areaCode%22,%22date%22:%22date%22,%22newPeopleVaccinatedFirstDoseByPublishDate%22:%22newPeopleVaccinatedFirstDoseByPublishDate%22,%22cumPeopleVaccinatedFirstDoseByPublishDate%22:%22cumPeopleVaccinatedFirstDoseByPublishDate%22%7D&format=csv",
    column_interest="SecondDose",
    group_id="ALL",
    column_group="TargetGroup",
    sheet_name="data",
    date_variable="YearWeekISO",
    country_variable="ReportingCountry",
    population_variable="Population",
    vaccine_variable="Vaccine",
    vaccine_allocation_dict={
        "COM": "vac1",
        "MOD": "vac1",
        "AZ": "vac2",
        "JANSS": "vac2",
    },  # do not use unknown
    number_decision_periods=10,
    weeks_decision_period=4,
    model_population=160000000,
    format_df="long",
    percentage_mRNA=57 / 187,
    percentage_vector=130 / 187,
):
    col = "newPeopleVaccinatedFirstDoseByPublishDate"
    df_uk = (pd.read_csv(url_uk, skipfooter=0))[["date", col]]
    df_uk = df_uk.iloc[::-1].reset_index(drop=True).dropna()
    weeks = np.repeat(range(2, int(np.floor(df_uk.shape[0] / 7)) + 2), 7)
    df_uk["week"] = np.concatenate(
        [weeks, np.repeat(np.max(weeks) + 1, df_uk.shape[0] - len(weeks))]
    )
    df_uk_grouped = df_uk.groupby(["week"]).sum()
    insert_df = pd.DataFrame([0, 0, 0], index=[0, 1, 2], columns=[col])
    df_uk_full = insert_df.append(df_uk_grouped).reset_index(drop=True)
    df_uk_full["vac1"] = df_uk_full[col] * percentage_mRNA
    df_uk_full["vac2"] = df_uk_full[col] * percentage_vector
    df_uk_wide = df_uk_full[["vac1", "vac2"]]

    df_raw = pd.read_excel(url, sheet_name=sheet_name, skipfooter=0)
    vaccines_to_use = vaccine_allocation_dict.keys()
    df_all = df_raw.loc[
        (df_raw[column_group] == group_id)
        * (df_raw[vaccine_variable].isin(vaccines_to_use))
    ]

    column_to_keep = [
        date_variable,
        country_variable,
        column_interest,
        vaccine_variable,
    ]
    df_column_reduced = df_all[column_to_keep]
    df_column_reduced["vaccine_model"] = df_column_reduced[vaccine_variable].map(
        vaccine_allocation_dict
    )

    df_grouped = df_column_reduced.groupby([date_variable, "vaccine_model"]).sum()

    df_wide = df_grouped.pivot_table(
        index=[date_variable], columns="vaccine_model", values=column_interest
    )
    df_wide_sum = df_wide.reset_index(drop=True) + df_uk_wide
    
    

    df_non_zero = df_wide_sum[(df_wide_sum > 200).all(1)].reset_index(drop=True)

    length_df = df_non_zero.shape[0]

    if length_df >= number_decision_periods * weeks_decision_period:
        df_periods = df_non_zero.iloc[range((number_decision_periods * weeks_decision_period))]
    else:
        diff = int(number_decision_periods * weeks_decision_period - length_df)
        df_periods = df_non_zero.append(df_non_zero.iloc[[-1] * diff]).reset_index(
            drop=True
        )

    sequence = np.repeat(
        list(range(int(number_decision_periods))), int(weeks_decision_period)
    )
    df_periods["decision_period"] = sequence
    population_size = df_raw[population_variable].unique().sum() + 66647.1 * 1000
    linear_pop_scaling = model_population / population_size

    df_supply_model = (
        df_periods.groupby("decision_period").sum()
        / weeks_decision_period
        / 7
        * linear_pop_scaling
    )
    # df_supply_model = df_supply_model.drop(labels=["index"], axis = 1)
    if format_df == "long":
        df_supply_model = df_supply_model.unstack()

    return df_supply_model

def create_inflow_from_data_sum(
    url="https://opendata.ecdc.europa.eu/covid19/vaccine_tracker/xlsx/data.xlsx",
    url_uk="https://coronavirus.data.gov.uk/api/v1/data?filters=areaType=overview&structure=%7B%22areaType%22:%22areaType%22,%22areaName%22:%22areaName%22,%22areaCode%22:%22areaCode%22,%22date%22:%22date%22,%22newPeopleVaccinatedFirstDoseByPublishDate%22:%22newPeopleVaccinatedFirstDoseByPublishDate%22,%22cumPeopleVaccinatedFirstDoseByPublishDate%22:%22cumPeopleVaccinatedFirstDoseByPublishDate%22%7D&format=csv",
    column_interest="SecondDose",
    group_id="ALL",
    column_group="TargetGroup",
    sheet_name="data",
    date_variable="YearWeekISO",
    country_variable="ReportingCountry",
    population_variable="Population",
    vaccine_variable="Vaccine",
    vaccine_allocation_dict={
        "COM": "vac1",
        "MOD": "vac1",
        "AZ": "vac2",
        "JANSS": "vac2",
    },  # do not use unknown
    number_decision_periods=10,
    weeks_decision_period=4,
    model_population=160000000,
    format_df="long",
    percentage_mRNA=57 / 187,
    percentage_vector=130 / 187,
):
    col = "newPeopleVaccinatedFirstDoseByPublishDate"
    df_uk = (pd.read_csv(url_uk, skipfooter=0))[["date", col]]
    df_uk = df_uk.iloc[::-1].reset_index(drop=True).dropna()
    weeks = np.repeat(range(2, int(np.floor(df_uk.shape[0] / 7)) + 2), 7)
    df_uk["week"] = np.concatenate(
        [weeks, np.repeat(np.max(weeks) + 1, df_uk.shape[0] - len(weeks))]
    )
    df_uk_grouped = df_uk.groupby(["week"]).sum()
    insert_df = pd.DataFrame([0, 0, 0], index=[0, 1, 2], columns=[col])
    df_uk_full = insert_df.append(df_uk_grouped).reset_index(drop=True)
    df_uk_full["vac1"] = df_uk_full[col] * percentage_mRNA
    df_uk_full["vac2"] = df_uk_full[col] * percentage_vector
    df_uk_wide = df_uk_full[["vac1", "vac2"]]

    df_raw = pd.read_excel(url, sheet_name=sheet_name, skipfooter=0)
    vaccines_to_use = vaccine_allocation_dict.keys()
    df_all = df_raw.loc[
        (df_raw[column_group] == group_id)
        * (df_raw[vaccine_variable].isin(vaccines_to_use))
    ]

    column_to_keep = [
        date_variable,
        country_variable,
        column_interest,
        vaccine_variable,
    ]
    df_column_reduced = df_all[column_to_keep]
    df_column_reduced["vaccine_model"] = df_column_reduced[vaccine_variable].map(
        vaccine_allocation_dict
    )

    df_grouped = df_column_reduced.groupby([date_variable, "vaccine_model"]).sum()

    df_wide = df_grouped.pivot_table(
        index=[date_variable], columns="vaccine_model", values=column_interest
    )
    df_wide_sum = df_wide.reset_index(drop=True) + df_uk_wide
    
    return {"europe":df_wide.reset_index(drop=True),
            "uk":df_uk_wide}