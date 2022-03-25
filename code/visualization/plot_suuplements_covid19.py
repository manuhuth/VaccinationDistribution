def plot_vacc_available_initial_pop(color_vac1 = "steelblue",
                                    color_vac2 = "C1",
                                    label_vac1 = "mRNA",
                                    label_vac2 = "Vector",
                                    linspace = np.linspace(0, 140, 5000 ),scale=1000000,
                                    figsize=(16,10)):
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    gs = GridSpec(1, 2, figure=fig)

    y_position = 1.08
    y_size = 24
    count_plot = 65

    time =  pd.date_range("2021-01-01", "2022-02-25", freq="7259s")

    vaccine_available = pd.DataFrame(
            {
                "vac1": np.repeat(
                    np.array(list(vaccine_inflow.values())[0 : (number_yy - 1)]), 1000
                ),
                "vac2": np.repeat(
                    np.array(
                        list(vaccine_inflow.values())[(number_yy - 1) : (2 * number_yy)]
                    ),
                    1000,
                ),
                "t": np.linspace(0, length, len(linspace)),
            }
        )  
    ax = fig.add_subplot(gs[0, 0])
    
    ax.plot(
            time,
            vaccine_available["vac1"] / scale,
            color=color_vac1,
            label=label_vac1,
    )
    ax.set_ylabel("Doses per day \nin millions")
    ax.fill_between(
            time,
            vaccine_available["vac1"] / scale,
            vaccine_available["vac2"] / scale,
            color=color_vac1,
            step="pre",
            alpha=0.6,
    )
    ax.plot(
            time,
            vaccine_available["vac2"] / scale,
            color=color_vac2,
            label=label_vac2,
    )
    ax.fill_between(
            time,
            vaccine_available["vac2"] / scale,
            color=color_vac2,
            step="pre",
            alpha=0.6,
    )
    #ax.set_xlabel("Weeks")
    ax.xaxis.set_tick_params(rotation=45)
    ax.yaxis.grid(alpha=0.6)
    
    ax.set_title("Available vaccines")
    ax.legend()
    ax.text(
        -0.05,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
    )
    count_plot += 1


    ax = fig.add_subplot(gs[0, 1])

    susceptible = (
            start_population["susceptible_countryA_vac0_t0"]
            + start_population["susceptible_countryB_vac0_t0"]
            + start_population["susceptible_countryC_vac0_t0"]
            + start_population["susceptible_countryD_vac0_t0"]
    )
    infected = (
            infectious_t0["infectious_countryA_vac0_virus1_t0"]
            + infectious_t0["infectious_countryB_vac0_virus1_t0"]
            + infectious_t0["infectious_countryC_vac0_virus1_t0"]
            + infectious_t0["infectious_countryD_vac0_virus1_t0"]
    )
    recovered = (
            recovered_t0["recovered_countryA_vac0_virus1_t0"]
            + recovered_t0["recovered_countryB_vac0_virus1_t0"]
            + recovered_t0["recovered_countryC_vac0_virus1_t0"]
            + recovered_t0["recovered_countryD_vac0_virus1_t0"]
    )

    sums = {"susceptible": susceptible, "infectious": infected, "recovered": recovered}
    dicts = [start_population, infectious_t0, recovered_t0]
    category_names = ["Belgium", "France", "Germany", "Uk"]
    states = ["susceptible", "infectious", "recovered"]
    results = {}
    for i in range(len(states)):
        d = dicts[i]
        ph_str = ""
        if states[i] != "susceptible":
            ph_str = "_virus1"
        results[states[i].capitalize()] = np.round(
            np.array(
                [
                    d[f"{states[i]}_countryA_vac0{ph_str}_t0"],
                    d[f"{states[i]}_countryB_vac0{ph_str}_t0"],
                    d[f"{states[i]}_countryC_vac0{ph_str}_t0"],
                    d[f"{states[i]}_countryD_vac0{ph_str}_t0"],
                ]
            )
            / sums[states[i]],
            2,
        )

    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = ["steelblue", "C1", "seagreen", "firebrick"] #plt.get_cmap("RdYlGn")(np.linspace(0.15, 0.85, data.shape[1]))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())
    ax.set_title("Relative initial populations")
    index = 0
    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        rects = ax.barh(
            labels, widths, left=starts, height=0.5, label=colname, color=category_colors[index], alpha = 0.4
        )

        #r, g, b, _ = color
        #text_color = "black" if r * g * b < 0.5 else "darkgrey"
        if colname != "Belgium":
            ax.bar_label(
                rects,
                label_type="center",
                fmt="%.2f%%",
                color=text_color,
                fontsize="small",
                padding=0,
            )
        # ax.legend(ncol=2, fontsize="small")
        ax.legend(ncol=len(category_names), bbox_to_anchor=(0, -0.2), loc="lower left")
        index += 1

    ax.text(
        -0.05,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
    )
    
    
    fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/supplements_vacc_inflow",
    bbox_inches="tight",
    )


    
    
plot_vacc_available_initial_pop(figsize=(20, 7))



#################################################################################################
def plot_used_NPIs():
        fig = plt.figure(constrained_layout=True, figsize=(20,14))
        gs = GridSpec(2, 2, figure=fig)
        
        
        grid_data=np.linspace(0, end_data, int(total_grid * end_data / length))
        grid_sim=np.linspace(
            end_data, length, total_grid - int(total_grid * end_data / length) - 1
        )
        text_x=end_data / 7 / 3 - 2
        text_y=0.07
        ylim=[0, 100]
        text_str="Vaccination \nperiod"
        text_lockdown_x=((length - end_data) / 2 + end_data) / 7 - 7
        text_lockdown_y=0.03
        text_lockdown_str="Constant \nNPIs"
        y_position = 1.08
        y_size = 24
        count_plot = 64

        time =  pd.date_range("2021-01-01", "2022-02-25", freq="6049s")
        colors = ["steelblue", "C1", "seagreen", "firebrick"]
        for j in range(len(countries)):
            if j <=1:
                ax = fig.add_subplot(gs[0, j])
            else:
                ax = fig.add_subplot(gs[1, j-2])
            count_plot += 1
            ax.text(
                -0.05,
                y_position,
                chr(count_plot),
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                weight="bold",
                size=y_size,
                )

            if j in [0, 2]:
                ax.set_ylabel("Degree of NPIs")
            country = areas[j]
            array = np.array([par_R[x] for x in par_R.keys() if country in x])
            spline_R = get_spline(
                array,
                periods=number_xx_R - 1,
                length=length / number_xx_R,
                total_length=length,
                grid_points=total_grid,
                transform=True,
            )
    
            #ax.set_xlabel("Weeks")
            ax.plot(time, 100 - 100*spline_R, color=colors[j])
 
            #ax.text(text_x, text_y, text_str)
            #ax.text(text_lockdown_x, text_lockdown_y, text_lockdown_str)
            ax.set_ylim(ylim)
            if j <3:
                ax.set_title(countries[j].capitalize())
            if j == 3:
                ax.set_title("United Kingdom")
            
            ax.fill_between(
                time,
                100 - 100*spline_R,
                step="pre",
                alpha=0.4,
                color=colors[j],
            )
            ax.xaxis.set_tick_params(rotation=45)
            ax.yaxis.grid(alpha=0.6)
 
            if j == 0:
                ax.legend(loc="upper right")
            fig.tight_layout(pad=1.5)
            fig.savefig(
            "/home/manuel/Documents/VaccinationDistribution/paper/images/supplements_npis_used",
            bbox_inches="tight",
            )

plot_used_NPIs()
