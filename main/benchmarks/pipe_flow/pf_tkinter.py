from Tkinter import *
import main.schemes.d2q9 as lb
import main.schemes.grid as gd
import numpy as np
import main.schemes.plot_and_save as ps
import os, inspect


class EntryWindow:
    def __init__(self):
        self.master = Tk()
        self.errormsg = Label(text='', fg='red')
        self.errormsg.pack()

        row1 = Frame(self.master)
        row1.pack(side=TOP, fill=X, padx=5, pady=5)
        row2 = Frame(self.master)
        row2.pack(side=TOP, fill=X, padx=5, pady=5)
        row3 = Frame(self.master)
        row3.pack(side=TOP, fill=X, padx=5, pady=5)
        row4 = Frame(self.master)
        row4.pack(side=TOP, fill=X, padx=5, pady=5)
        row5 = Frame(self.master)
        row5.pack(side=TOP, fill=X, padx=5, pady=5)
        row6 = Frame(self.master)
        row6.pack(side=TOP, fill=X, padx=5, pady=5)

        self.velocity_entry = Entry(row1,
                                    textvariable=StringVar(self.master, value='0.1'))
        self.pressure_entry = Entry(row2,
                                     textvariable=StringVar(self.master, value='1'))
        self.viscosity_entry = Entry(row3,
                                     textvariable=StringVar(self.master, value='0.1'))
        self.height_entry = Entry(row4,
                                  textvariable=StringVar(self.master, value='50'))
        self.width_entry = Entry(row5,
                                 textvariable=StringVar(self.master, value='100'))
        self.time_entry = Entry(row6,
                                textvariable=StringVar(self.master, value='1000'))

        Label(row1, text='Velocity', width=15).pack(side=LEFT)
        Label(row2, text='Pressure', width=15).pack(side=LEFT)
        Label(row3, text='Viscosity', width=15).pack(side=LEFT)
        Label(row4,text='Height', width=15).pack(side=LEFT)
        Label(row5,text='Length', width=15).pack(side=LEFT)
        Label(row6,text='Time', width=15).pack(side=LEFT)

        self.reynolds_label = Label(text='')
        self.velocity_entry = Entry(row1,
                                    textvariable=StringVar(self.master, value='0.1'),
                                    validate="focusout",
                                    validatecommand=(self.master.register(self.validate_number), '%P', 'velocity'),
                                    invalidcommand=(self.master.register(self.invalid_number), '%P', 'velocity'))
        self.velocity_entry.pack(side=RIGHT)
        self.pressure_entry = Entry(row2,
                                     textvariable=StringVar(self.master, value='1'),
                                     validate="focusout",
                                     validatecommand=(self.master.register(self.validate_number), '%P', 'pressure'),
                                     invalidcommand=(self.master.register(self.invalid_number), '%P', 'pressure'))
        self.pressure_entry.pack(side=RIGHT)
        self.viscosity_entry = Entry(row3,
                                     textvariable=StringVar(self.master, value='0.1'),
                                     validate="focusout",
                                     validatecommand=(self.master.register(self.validate_number), '%P', 'viscosity'),
                                     invalidcommand=(self.master.register(self.invalid_number), '%P', 'viscosity'))
        self.viscosity_entry.pack(side=RIGHT)
        self.height_entry = Entry(row4,
                                  textvariable=StringVar(self.master, value='50'),
                                  validate="focusout",
                                  validatecommand=(self.master.register(self.validate_number), '%P', 'height'),
                                  invalidcommand=(self.master.register(self.invalid_number), '%P', 'height'))
        self.height_entry.pack(side=RIGHT)
        self.width_entry = Entry(row5,
                                 textvariable=StringVar(self.master, value='100'),
                                 validate="focusout",
                                 validatecommand=(self.master.register(self.validate_number), '%P', 'length'),
                                 invalidcommand=(self.master.register(self.invalid_number), '%P', 'length'))
        self.width_entry.pack(side=RIGHT)
        self.time_entry = Entry(row6,
                                  textvariable=StringVar(self.master, value='1000'),
                                  validate="focusout",
                                  validatecommand=(self.master.register(self.validate_number), '%P', 'time'),
                                  invalidcommand=(self.master.register(self.invalid_number), '%P', 'time'))
        self.time_entry.pack(side=RIGHT)
        self.reynolds_label.pack()
        Button(self.master, text="OK", command=self.create_window).pack()
        mainloop()

    def validate_number(self, p, parameter):
        self.errormsg.config(text='')
        if parameter == 'velocity':
            x = False
            try:
                x = (float(p) <= 0.2)
                ren = float(self.height_entry.get()) * \
                      float(self.velocity_entry.get()) / \
                      float(self.viscosity_entry.get())
                self.reynolds_label.config(text='Re: ' + str(ren))
                # print str(ren) + ' vel'
            except ValueError:
                pass
            return x  # always return True or False
        if parameter == 'viscosity':
            x = False
            try:
                x = (float(p) >= 0.03)
                ren = float(self.height_entry.get()) * \
                      float(self.velocity_entry.get()) / \
                     float(self.viscosity_entry.get())
                self.reynolds_label.config(text='Re: ' + str(ren))
                # print str(ren) + ' vis'
            except ValueError:
                pass
            return x  # always return True or False
        if parameter == 'height':
            x = False
            try:
                x = (float(p) <= 500)
                ren = float(self.height_entry.get()) * \
                      float(self.velocity_entry.get()) / \
                     float(self.viscosity_entry.get())
                # print str(ren) + ' len'
                self.reynolds_label.config(text='Re: ' + str(ren))
            except ValueError:
                pass
            return x  # always return True or False

    def invalid_number(self, p, parameter):
        if parameter == 'velocity':
            try:
                if float(p) > 0.2:
                    self.errormsg.config(text='Velocity is > 0.2')
            except ValueError:
                self.errormsg.config(text='Not a number')
        if parameter == 'pressure':
            try:
                if float(p) < 0:
                    self.errormsg.config(text='Pressure is < 0')
            except ValueError:
                self.errormsg.config(text='Not a number')
        if parameter == 'viscosity':
            try:
                if float(p) < 0.03:
                    self.errormsg.config(text='Viscosity is < 0.03')
            except ValueError:
                self.errormsg.config(text='Not a number')
        if parameter == 'height':
            try:
                if float(p) > 500:
                    self.errormsg.config(text='Height is > 500')
            except ValueError:
                self.errormsg.config(text='Not a number')
        if parameter == 'length':
            try:
                if float(p) > 500:
                    self.errormsg.config(text='Length is > 500')
            except ValueError:
                self.errormsg.config(text='Not a number')
        if parameter == 'time':
            try:
                if float(p) <=0:
                    self.errormsg.config(text='Time is <= 0')
            except ValueError:
                self.errormsg.config(text='Not a number')

    def create_window(self):
        entries = (('Velocity', float(self.velocity_entry.get())),
                   ('Pressure', float(self.pressure_entry.get())),
                   ('Viscosity', float(self.viscosity_entry.get())),
                   ('Width', int(self.height_entry.get())),
                   ('Length', int(self.width_entry.get())),
                   ('Time', int(self.time_entry.get())))
        self.master.destroy()
        new_window = MainWindow(entries)


class MainWindow(object):
    def __init__(self, entries):

        self.start_pipe_flow(entries)

    def start_pipe_flow(self, entries):
        print entries
        n = entries[4][1]  # width  - 0x
        m = entries[3][1]  # length - 0y
        v = entries[2][1]  # viscosity
        t_f = entries[5][1]  # final time
        tau = v * 3. + 0.5
        u_x = entries[0][1]
        p = entries[1][1]
        reynolds = n * u_x / v
        user_grid = gd.Grid(n, m)
        # print isinstance(velocity,lv.LatticeVelocity)
        # boundary conditions
        bc = self.pipe_vp_bc(grid=user_grid, v_left=u_x, p_right=p)
        lat_bol = lb.D2Q9(grid=user_grid,
                          iterations=t_f,
                          boundary=bc,
                          relaxation_time=tau,
                          reynolds=reynolds)

        path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        ps.save_all(lat_bol, path=path)
        ps.plot_streamlines(lat_bol, path=path)

    def pipe_pp_bc(self, grid, p_right, p_left):
        boundaries = []
        u = gd.LatticeVelocity(0, 0)  # zero velocity
        assert isinstance(p_right, float), "Incorrect type"
        assert isinstance(p_left, float), "Incorrect type"
        for i in range(0, grid.width):
            boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=u))
            boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=u))
        for i in range(0, grid.height):
            boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_pressure', value=p_left))
            boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_pressure', value=p_right))
        return boundaries

    def pipe_vp_bc(self, grid, p_right, v_left):
        boundaries = []
        u = gd.LatticeVelocity(0, 0)  # zero velocity
        assert isinstance(p_right, float), "Incorrect type"
        assert isinstance(v_left, float), "Incorrect type"

        for i in range(0, grid.width):
            boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=u))
            boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=u))
        for i in range(0, grid.height):
            boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_velocity', value=gd.LatticeVelocity(v_left, 0)))
            boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_pressure', value=p_right))
        return boundaries

    def pipe_vv_bc(self, grid, v_right, v_left):
        boundaries = []
        u = gd.LatticeVelocity(0, 0)  # zero velocity
        assert isinstance(v_right, float), "Incorrect type"
        assert isinstance(v_left, float), "Incorrect type"

        for i in range(0, grid.width):
            boundaries.append(gd.Boundary(grid.top_faces[i], 'zoe_he_velocity', value=u))
            boundaries.append(gd.Boundary(grid.bottom_faces[i], 'zoe_he_velocity', value=u))
        for i in range(0, grid.height):
            boundaries.append(gd.Boundary(grid.left_faces[i], 'zoe_he_velocity', value=gd.LatticeVelocity(v_left, 0)))
            boundaries.append(gd.Boundary(grid.right_faces[i], 'zoe_he_velocity', value=gd.LatticeVelocity(v_right, 0)))
        return boundaries

if __name__ == "__main__":
    kv = EntryWindow()


