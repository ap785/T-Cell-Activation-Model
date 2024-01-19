import matplotlib.pyplot as plt
import numpy as np

# Define time vector and convert to hours
t = np.linspace(0, 20 * 3600, 1000)
t_hours = t / 3600

# Define cytokines and their data
cytokines = ['IFN-g', 'IL-2', 'TNF-a']
cytokine_data = {
    'IFN-g': [x_no_costim[:, 3], x_with_cd2[:, 3], x_with_icam[:, 3], x_with_cd28[:, 3], x_with_cd27[:, 3], x_with_41bb[:, 3]],
    'IL-2': [x_no_costim_IL2[:, 3], x_with_cd2_IL2[:, 3], x_with_icam_IL2[:, 3], x_with_cd28_IL2[:, 3], x_with_cd27_IL2[:, 3], x_with_41bb_IL2[:, 3]],
    'TNF-a': [x_no_costim_TNF[:, 3], x_with_cd2_TNF[:, 3], x_with_icam_TNF[:, 3], x_with_cd28_TNF[:, 3], x_with_cd27_TNF[:, 3], x_with_41bb_TNF[:, 3]],
}

# Loop through cytokines
for cytokine in cytokines:
    # Create a new subplot for each cytokine
    fig, ax = plt.subplots(figsize=(10, 6))

    # Loop through costimulation scenarios
    scenarios = ['No Costim', 'CD2', 'ICAM', 'CD28', 'CD27', '4-1BB']
    for label, data in zip(scenarios, cytokine_data[cytokine]):
        cumulative_output = np.cumsum(data)  # Cumulative sum for the specific cytokine and costimulation
        ax.bar(label, cumulative_output[-1], label=label)  # Use the last value for cumulative output

    ax.set_title(f'Total Output - {cytokine}')
    ax.set_ylabel('Cumulative Output O')
    ax.set_xlabel('Costimulation Scenarios')
    ax.legend()

    plt.show()

