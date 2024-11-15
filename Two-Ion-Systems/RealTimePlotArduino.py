import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
xs = []
ys = []



def animate(i, xs, ys):

    # Read temperature (Celsius) from TMP102
    rand = np.random.rand()

    # Add x and y to lists
    xs.append(i)
    ys.append(rand)

    # Limit x and y lists to 20 items
    xs = xs[-20:]
    ys = ys[-20:]

    # Draw x and y lists
    ax.clear()
    ax.plot(xs, ys)

    # Format plot
    plt.xticks(rotation=45, ha='right')
    plt.subplots_adjust(bottom=0.30)
    plt.title('TMP102 Temperature over Time')
    plt.ylabel('Temperature (deg C)')

ani = animation.FuncAnimation(fig, animate, fargs=(xs, ys), interval=1000)
plt.show()