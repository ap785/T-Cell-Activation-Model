import matplotlib.pyplot as plt
import numpy as np

# Define cytokines and their data
cytokines = ['IFN-g', 'IL-2', 'TNF-a']
cytokine_data = {
    'IFN-g': [x_no_costim[:, 3], x_with_cd2[:, 3], x_with_icam[:, 3], x_with_cd28[:, 3], x_with_cd27[:, 3], x_with_41bb[:, 3]],
    'IL-2': [x_no_costim_IL2[:, 3], x_with_cd2_IL2[:, 3], x_with_icam_IL2[:, 3], x_with_cd28_IL2[:, 3], x_with_cd27_IL2[:, 3], x_with_41bb_IL2[:, 3]],
    'TNF-a': [x_no_costim_TNF[:, 3], x_with_cd2_TNF[:, 3], x_with_icam_TNF[:, 3], x_with_cd28_TNF[:, 3], x_with_cd27_TNF[:, 3], x_with_41bb_TNF[:, 3]],
}

# Check the number of cytokines and time intervals
num_cytokines = len(cytokines)
num_intervals = 2  # First 12 Hours and Last 8 Hours

# Create subplots for each cytokine
fig, axs = plt.subplots(num_cytokines, num_intervals, figsize=(12, 14), sharex=True, sharey='row')

# Define the color for the first 12 hours (lilac)
first_12_hours_color = ['thistle']

# Define y-axis limits for each cytokine
ylim_dict = {
    'IFN-g': (0, 0.000020),
    'IL-2': (0, 0.000025),
    'TNF-a': (0, 0.000025)
}

# Loop through cytokines and time intervals
for i, cytokine in enumerate(cytokines):
    for j, time_interval in enumerate(['First 12 Hours', 'Last 8 Hours']):
        ax = axs[i, j]

        # Calculate the average rate of production for the selected time interval for each scenario
        is_first_12_hours = (t_hours <= 12) if time_interval == 'First 12 Hours' else (t_hours > 12)
        avg_dO_dt = [np.mean(np.gradient(data, t)[is_first_12_hours]) for data in cytokine_data[cytokine]]

        # Create bar chart for the average rate of cytokine production with reduced gap
        scenarios = ['pMHC', '+ CD2', '+ ICAM-1', '+ CD28', '+ CD27', '+ 4-1BB']
        bar_width = 0.9  # Adjust the width to reduce the gap
        index = np.arange(len(scenarios))

        # Set color for the first 12 hours
        if time_interval == 'First 12 Hours':
            bars = ax.bar(index, avg_dO_dt, bar_width, color=first_12_hours_color)
        else:
            bars = ax.bar(index, avg_dO_dt, bar_width, color='lightskyblue')

        ax.set_title(f'{cytokine} Rate - {time_interval}')
        ax.set_ylabel('Rate of Cytokine Production')
        
        # Set y-axis limits for each cytokine
        ax.set_ylim(ylim_dict[cytokine])
        
        ax.set_xticks(index)
        ax.set_xticklabels(scenarios)

        # Add values to the bars in scientific notation
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f'{height:.1e}', xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords='offset points', ha='center', va='bottom')

# Adjust spacing between subplots
plt.tight_layout()

# Save or display the figure
plt.savefig('cytokine_rate_plots_2x3.pdf', dpi=300)  # Save the figure as an image file
plt.show()  # Display the figure

