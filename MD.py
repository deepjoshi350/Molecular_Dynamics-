import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as ani

# Initialising Particles

N = 100
L = 20
positions = []
grid = np.arange(0, L, 2)
while len(positions) < N:
    x = np.random.choice(grid)
    y = np.random.choice(grid)

    if (x, y) not in positions:
        positions.append((x, y))

positions = np.array(positions)
#print(positions[0:20])

def Distance(positions, L):

    N = positions.shape[0]
    dx_matrix = np.zeros((N, N))
    dy_matrix = np.zeros((N, N))
    dr_matrix = np.zeros((N, N))

# for periodic boundary

    for i in range(N):
        for j in range(i + 1, N):
            dx = positions[j, 0] - positions[i, 0]
            dy = positions[j, 1] - positions[i, 1]

            dx -= L * np.round(dx / L)
            dy -= L * np.round(dy / L)

            dr = np.sqrt(dx**2 + dy**2)

          
            dx_matrix[i, j] = dx
            dy_matrix[i, j] = dy
            dr_matrix[i, j] = dr
             
            
            dr_matrix[j, i] = dr
            dx_matrix[j, i] = -dx
            dy_matrix[j, i] = -dy

# stacking both matrices for easier loops

    r_vectors = np.stack((dx_matrix, dy_matrix), axis=-1)

    return dr_matrix, r_vectors

def net_forces(dr_matrix, r_vectors, positions, sigma=1.0, epsilon=1.0):
    N = dr_matrix.shape[0]
    net_forces = np.zeros((N, 2))

    for i in range(N):
        total_force = np.array([0.0, 0.0])

        for j in range(N):

            if i == j:
               continue

            dr = dr_matrix[i, j]
            r_vector = r_vectors[i, j]
            # if dr < 0.5:
             # magnitude = 24*epsilon * (
             # 2*(sigma/dr)**12 - (sigma/dr)**6
             # ) / dr**2

             # print("dr =", dr)
             # print("r_vector =", r_vector)
             # print("magnitude =", magnitude)
             # print("force =", magnitude * r_vector)
             # raise ValueError("force sign")
            
            

            magnitude = 24*epsilon * (2 * (sigma / dr)**12 - (sigma / dr)**6)/ dr**2
            total_force -= magnitude * r_vector
        net_forces[i] = total_force

    return net_forces

# Total energy

def total_energy(positions, velocities, L, Mass, sigma=1.0, epsilon=1.0):



    dr_matrix, _ = Distance(positions, L)

    # Kinetic Energy
    KE = 0.5 * Mass * np.sum(velocities**2)

    # Potential Energy
    PE = 0.0
    N = len(positions)

    for i in range(N):
        for j in range(i+1, N):

            r = dr_matrix[i,j]

            PE = PE + 4 * epsilon * (
                (sigma/r)**12 -
                (sigma/r)**6
            )

    return KE + PE, KE, PE

# For velocity verlet

def velocity_verlet_step(positions, V, A, dt, L, mass=1.0):
    positions = positions + V * dt + 0.5 * A * dt**2
    positions = positions % L
    #positions += V * dt + 0.5 * A * dt**2
    #positions %= L

    dr_matrix, r_vectors = Distance(positions, L)
    forces = net_forces(dr_matrix, r_vectors, positions)
    A_new = forces / mass

    V = V + 0.5 * (A + A_new) * dt
    A = A_new

    return positions, V, A



# defining velocity matrix, dt and steps

#V = np.zeros((N,2))
kT = 1.0
Mass = 1.0

# choosing initial velocities from maxwell distribution
V = np.random.normal(loc=0,
    scale=np.sqrt(kT/Mass),
    size=(N,2))

# Remove COM drift
V -= np.mean(V, axis=0)
dt = 0.0005

Steps = 20000


# Initial values of relative distance and force matrices
dr_matrix, r_vectors = Distance(positions, L)
forces = net_forces(dr_matrix, r_vectors, positions)
A = forces / Mass

energy_history = []
KE_history = []
PE_history = []

position_history = []
#velocity_history = []

E, KE, PE = total_energy(positions,V,L,Mass)




# Main Loop

save_every = 10

for step in range(Steps):

       positions, V, A = velocity_verlet_step(positions, V, A, dt, L, Mass)
    # dr_matrix, _ = Distance(positions, L)

       if step % save_every == 0:
         position_history.append(positions.copy())
     #velocity_history.append(V.copy())

         energy_history.append(E)
         KE_history.append(KE)
         PE_history.append(PE)
         E, KE, PE = total_energy( positions, V, L, Mass)

    



position_history = np.array(position_history)
#velocity_history = np.array(velocity_history)
energy_history = np.array(energy_history)
KE_history = np.array(KE_history)
PE_history = np.array(PE_history)

#print(np.min(PE_history))
#print(np.max(PE_history))

# plotting energies

plt.figure(figsize=(8,5))

plt.plot(energy_history, label='Total Energy')
plt.plot(KE_history, label='Kinetic Energy')
plt.plot(PE_history, label='Potential Energy')

plt.xlabel('Step')
plt.ylabel('Energy')
plt.legend()
plt.grid()

plt.show()

#plt.plot(min_r)
#plt.xlabel("Step")
#plt.ylabel("Minimum pair distance")
#plt.grid()
#plt.show()

# Relative energy error

E0 = energy_history[0]

relative_error = np.abs((energy_history - E0) / E0)

plt.figure(figsize=(8,5))
plt.plot(relative_error)

plt.xlabel('Step')
plt.ylabel('Relative Energy Error')
plt.yscale('log')
plt.grid()

plt.show()

# for animation

from matplotlib.animation import FuncAnimation

def animation(position_history, L):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_xlabel('X position')
    ax.set_ylabel('Y position')
    ax.set_title('Molecular Dynamics Simulation')
    ax.grid(True, linestyle='--', alpha=0.5)

    # scatter plot for particles and lines for trajectories
    scat = ax.scatter([], [], color='blue', s=30)
    lines = [ax.plot([], [], color='gray', alpha=0.5, linewidth=1)[0]
         for _ in range(N)]

    def init():
        scat.set_offsets(np.empty((0, 2)))
        for line in lines:
            line.set_data([], [])
        return [scat] + lines

    def update(frame):

        pos = position_history[frame]
        scat.set_offsets(pos)
        #print(frame)

        for i in range(N):
            x = position_history[:frame+1, i, 0]
            y = position_history[:frame+1, i, 1]


            dx = np.diff(x)
            dy = np.diff(y)
            jumps = (np.abs(dx) > L/2) | (np.abs(dy) > L/2)

            if np.any(jumps):

                last_jump = np.where(jumps)[0][-1]
                lines[i].set_data(x[last_jump+1:], y[last_jump+1:])
            else:
                lines[i].set_data(x, y)

        return [scat] + lines

    #anim = FuncAnimation(fig, update, frames=len(position_history), init_func=init, blit=True)
    #return anim
    #plt.close()

#from IPython.display import HTML

#ani = animation(position_history, L)

# Save GIF

#ani.save(
    #"md_simulation.gif",
    #writer="pillow",
   # fps=30
#)

# Preview inside Colab

#from IPython.display import Image

#Image(filename="md_simulation.gif")
# print("Animation saved")

#from google.colab import files

#files.download("md_simulation.gif")

# print(position_history)
