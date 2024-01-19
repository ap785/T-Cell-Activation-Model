import matplotlib.pyplot as plt
import numpy as np

# Define cytokines and their data
cytokines = ['IFN-gamma', 'IL-2', 'TNF-alpha']  # Updated cytokine names
cytokine_data = {
    'IFN-gamma': [x_no_costim[:, 3], x_with_cd2[:, 3], x_with_icam[:, 3], x_with_cd28[:, 3], x_with_cd27[:, 3], x_with_41bb[:, 3]],
    'IL-2': [x_no_costim_IL2[:, 3], x_with_cd2_IL2[:, 3], x_with_icam_IL2[:, 3], x_with_cd28_IL2[:, 3], x_with_cd27_IL2[:, 3], x_with_41bb_IL2[:, 3]],
    'TNF-alpha': [x_no_costim_TNF[:, 3], x_with_cd2_TNF[:, 3], x_with_icam_TNF[:, 3], x_with_cd28_TNF[:, 3], x_with_cd27_TNF[:, 3], x_with_41bb_TNF[:, 3]],
}

# Check the number of cytokines and time intervals
num_cytokines = len(cytokines)
num_intervals = 2  # First 12 Hours and Last 8 Hours

# Create a grid of subplots for each cytokine
fig, axes = plt.subplots(1, num_cytokines, figsize=(15, 5))  # Adjust the figsize as needed

# Loop through cytokines
for i, cytokine in enumerate(cytokines):
    ax = axes[i]
    ax.set_title(cytokine)
    
    # Initialize lists to store data
    x_data = []  # Rate in the first 12 hours
    y_data = []  # Rate in the last 8 hours

    # Loop through scenarios
    scenarios = ['pMHC', '+ CD2', '+ LFA-1', '+ CD28', '+ CD27', '+ 4-1BB']
    
    # Define specific colors for each scenario
    scenario_colors = {
        'pMHC': '#808080',
        '+ CD2': '#217DB6',
        '+ LFA-1': '#664C99',
        '+ CD28': '#34A755',
        '+ CD27': '#E73735',
        '+ 4-1BB': '#F08234',
}
    

    
    for scenario in scenarios:
        # Extract data for the current cytokine and scenario
        data = cytokine_data[cytokine][scenarios.index(scenario)]
        
        # Calculate the average rate of production for the first 12 hours and last 8 hours
        avg_dO_dt_first_12 = np.mean(np.gradient(data[:int(len(data)/2)], t[:int(len(t)/2)]))
        avg_dO_dt_last_8 = np.mean(np.gradient(data[int(len(data)/2):], t[int(len(t)/2):]))
        
        x_data.append(avg_dO_dt_first_12)
        y_data.append(avg_dO_dt_last_8)

        # Get the specific color for the scenario
        color = scenario_colors.get(scenario, 'black')
        
        # Increase the size of each point (adjust the value as needed)
        point_size = 70
        
        # Create a scatter plot for the point with the specific color and larger size
        ax.scatter(avg_dO_dt_first_12, avg_dO_dt_last_8, c=[color], s=point_size, marker='o', label=scenario)
        
        # Add text annotation with the scenario name
#         ax.annotate(scenario, (avg_dO_dt_first_12, avg_dO_dt_last_8), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=8, color=color)

        # Add text annotation with the scenario name
        label_offset = 10  # Offset for the label text
        ax.annotate(scenario, (avg_dO_dt_first_12, avg_dO_dt_last_8),
                    textcoords="offset points", xytext=(0, label_offset), ha='center', fontsize=8, color=color)
    
    
    
    # Adjust label positions to avoid overlaps
    adjust_label_positions(ax, x_data, y_data)
        
    ax.set_xlabel('Production Rate (Early Collection)')
    ax.set_ylabel('Production Rate (Late Collection)')

    # Set axis limits to go from 0 to 2.5
    ax.set_xlim(0, 0.00003)
    ax.set_ylim(0, 0.00003)

    # Set equal aspect ratio
    ax.set_aspect('equal')
    
    # Plot a dashed line through x = y
    ax.plot([0, 2.5], [0, 2.5], linestyle='--', color='gray')

    ax.grid(True)
#     ax.legend()

# Adjust spacing between subplots
plt.tight_layout()

# Save the scatter plots as an EPS image file
plt.savefig('cytokine_scatter_plots.eps', format='eps')
plt.show()

