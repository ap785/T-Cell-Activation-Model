import matplotlib.pyplot as plt

costim_data = {
    'No Costim': {'data': x_no_costim[:, 3], 'color': 'orange'},
    'CD2 Costim': {'data': x_with_cd2[:, 3], 'color': 'green'},
    'CD28 Costim': {'data': x_with_icam[:, 3], 'color': 'purple'},
    'CD28 Costim': {'data': x_with_cd28[:, 3], 'color': 'red'},
    'CD28 Costim': {'data': x_with_cd27[:, 3], 'color': 'yellow'},
    '4-1BB Costim': {'data': x_with_41bb[:, 3], 'color': 'blue'}
}

plt.figure(figsize=(10, 6))

for label, params in costim_data.items():
    plt.plot(t_hours, params['data'], label=label, color=params['color'])

plt.xlabel('Time (s)')
plt.ylabel('Output O')
plt.title('Comparison of Outputs from Different Costimulatory Receptors')
plt.legend()
plt.grid(False)

# Remove top and right spines from the plot
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.show()

