"""
Code initially taken from : https://stackoverflow.com/questions/49710067/generating-multiple-moving-objects-in-the-tkinter-canvas
Then modified for this simulation of molecular dynamics.
The data is taken from the function simulation_data(t, step, steps), whose output is:
positions: nd.array(n_atoms, coordinates, time_step)
velocities: nd.array(n_atoms, coordinates, time_step)
Here we only use the positions. They are received by the class MyCanvas(window, positions)
This executes an animation that iterates over all time_step coordinate
"""

import tkinter as tk
import time

class Ball:

    def __init__(self, atom_index, positions, time,coord_1,coord_2):
        self.xpos = positions[atom_index, coord_1, time]
        self.ypos = positions[atom_index, coord_2, time]
        print(self.xpos, self.ypos)

class MyCanvas(tk.Canvas):
    """
    Class that creates an animation.

    Parameters:
    -----------
        master: tkinter.Tk() object
            Window where animation is executed
        positions: nd.array
            Positions array with shape (n_atoms, coordinates, time_steps)
    """

    def __init__(self, master, positions, coord_1, coord_2, sleep_time):
        """

        """
        super().__init__(master, width=1200, height=400, bg="snow2", bd=0, highlightthickness=0, relief="ridge")
        self.pack()

        self.number_atoms = len(positions)
        self.balls = []
        self.bs = []
        self.sleep_time = sleep_time
        self.time = 0
        self.atom_index = 0
        self.coord_1 = coord_1
        self.coord_2 = coord_2
        self.positions = positions
        self.steps = len(self.positions[0,0])

        self.run()

    def check_steps(self):
        """
        Class that checks if time step has reached the end of the simulation.
        If so, it starts again.
        """
        if (self.time+1)%self.steps == 0:
            self.time = 0

    def run(self):
        """
        Execute time step.
        This function is called recursively using self.after(time(ms), self.run) .
        """
        self.check_steps()

        colors = ["red", "blue", "green", "yellow"]
        #Create atoms at position corresponding at time_step = self.time
        for _ in range(self.number_atoms):
            ball = Ball(self.atom_index, self.positions, self.time, self.coord_1, self.coord_2)
            self.atom_index += 1
            self.balls.append(ball)
            self.bs.append(self.create_oval(ball.xpos - 10,
                                            ball.ypos - 10,
                                            ball.xpos + 10,
                                            ball.ypos + 10,
                                            tags="atoms",
                                            fill=colors[(self.atom_index-1)%4]))

        self.atom_index = 0
        self.time += 1

        print("time step: "+str(self.time))

        super().update()
        super().delete("atoms")
        time.sleep(self.sleep_time)
        super().update()

        self.after(1, self.run)
