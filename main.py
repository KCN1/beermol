# version 0.4 alpha

from sys import argv
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
from vtkmodules.vtkCommonCore import (
    vtkPoints,
    vtkUnsignedCharArray
)
from vtkmodules.vtkCommonDataModel import (
    vtkCellArray,
    vtkLine,
    vtkPolyData
)
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkGlyph3DMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)


class RefreshWindow:

    def __init__(self, filename):
        self.filename = filename
        self.reset_flag = False
        self.renderer = vtkRenderer()
        self.renderer.SetBackground(0.0, 0.0, 0.0)
        self.window = vtkRenderWindow()
        self.window.SetSize(800, 800)
        self.window.AddRenderer(self.renderer)

    def from_gauss(self):

        n, description, p, el = 0, '', [], []

        with open(self.filename, 'r') as file:

            flag = 'ready'  # start/stop flag for reading data
            descr_flag = 'ready'  # start/stop flag for reading description
            geometries = 0  # to store all geometries (for next version)
            n0, p0, el0 = 0, [], []  # current number of atoms, coordinates, elements

            atom_numbers = (
                'n', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F',
                'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
                'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
                'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
                'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr',
                'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
                'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',
                'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
                'Fm', 'Md', 'No', 'Lr'
            )

            for line in file:

                if 'Coordinates' in line:
                    flag = 'set'
                    geometries += 1

                elif flag == 'set' and '-------' in line:
                    flag = 'go'

                elif flag == 'go' and '-------' in line:
                    flag = 'ready'
                    n = n0  # last number of atoms
                    p = p0.copy()  # last coordinates
                    el = el0.copy()  # last list of els
                    n0, p0, el0 = 0, [], []

                elif flag == 'go' and '-------' not in line:
                    s = line.split()
                    p0.append((float(s[-3]), float(s[-2]), float(s[-1])))
                    el0.append(atom_numbers[int(s[1])])
                    n0 += 1
                else:
                    pass

                if 'l101.exe' in line:
                    descr_flag = 'set'
                elif descr_flag == 'set':
                    descr_flag = 'go'
                elif descr_flag == 'go':
                    descr_flag = 'ready'
                    description = line.strip()
                else:
                    pass

        return n, description, el, p

    # def from_xyz(self):
    #
    #     n, description, p, el = 0, '', [], []
    #
    #     with open(self.filename, 'r') as file:
    #         n = int(file.readline().strip())
    #         description = file.readline().strip()
    #
    #         # create a coordinate array
    #         # create an element array
    #
    #         for i in range(n):
    #             s = file.readline().split()
    #             p_i = (float(s[-3]), float(s[-2]), float(s[-1]))
    #             p.append(p_i)
    #             el.append(s[0])
    #
    #     return n, description, el, p

    def from_trj_xyz(self):

        n, description, p, el = 0, '', [], []
        geometries = 0
        n0, p0, el0 = 0, [], []  # current number of atoms, coordinates, elements

        with open(self.filename, 'r') as file:
            line = file.readline()

            while line:
                n0 = int(line.strip())
                description = file.readline().strip()

                # create a coordinate array
                # create an element array

                for i in range(n0):
                    s = file.readline().split()
                    p_i = (float(s[-3]), float(s[-2]), float(s[-1]))
                    p0.append(p_i)
                    el0.append(s[0])

                n = n0
                p = p0.copy()
                el = el0.copy()
                n0, p0, el0 = 0, [], []
                geometries += 1
                line = file.readline()

        return n, description, el, p

    # noinspection PyUnusedLocal
    def read_xyz(self, caller=None, event=None):
        """Creates separate mappers for lines and spheres from a coordinate list"""

        if self.filename[-4:] == '.xyz':
            n, description, el, p = self.from_trj_xyz()
        elif self.filename[-4:] in ('.log', '.out'):
            n, description, el, p = self.from_gauss()
        else:
            n, description, el, p = 0, '', [], []
        points = vtkPoints()

        if not n:
            self.reset_flag = True

        for atom in p:
            points.InsertNextPoint(atom)

        hyper_colors = {
            'C': 'Cyan',
            'N': 'Blue',
            'O': 'Red',
            'F': 'Yellow',
            'Na': 'Purple',
            'P': 'Yellow',
            'S': 'Yellow',
            'Cl': 'Yellow',
            'K': 'Purple',
            'Fe': 'Red',
            'Co': 'Blue',
            'Cu': 'Green',
            'Br': 'Yellow',
            'I': 'Red',
            'Au': 'Yellow'
        }

        # define element colors
        el_color = ['White'] * n
        for i in range(n):
            if el[i] in hyper_colors:
                el_color[i] = hyper_colors[el[i]]
            else:
                pass

        radii = {
            'H': 0.35, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
            'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
            'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07,
            'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76,
            'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.61,
            'Fe': 1.52, 'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22,
            'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20,
            'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75,
            'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42,
            'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39,
            'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40, 'Cs': 2.44,
            'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01,
            'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94,
            'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87,
            'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.70, 'W': 1.62, 'Re': 1.51,
            'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
            'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 'At': 1.50,
            'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06,
            'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87, 'Am': 1.80,
            'Cm': 1.69, 'Bk': 1.68, 'Cf': 1.68, 'Es': 1.65, 'Fm': 1.67,
            'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'n': 0.0
        }

        # create a connectivity array
        bond_length = 2.1
        connectivity = []
        connected = []

        for i in range(n):
            for j in range(i + 1, n):
                if el[i] in radii and el[j] in radii:
                    bond_length = (radii[el[i]] + radii[el[j]]) * 1.1
                if (p[i][0] - p[j][0]) ** 2 + (p[i][1] - p[j][1]) ** 2 + (p[i][2] - p[j][2]) ** 2 < bond_length ** 2:
                    # add virtual points in the middle of each bond
                    pv = ((p[i][0] + p[j][0]) / 2, (p[i][1] + p[j][1]) / 2, (p[i][2] + p[j][2]) / 2)
                    p.append(pv)
                    points.InsertNextPoint(pv)
                    # connect last pair of bonded atoms to the last virtual point
                    connectivity.append([i, len(p) - 1])
                    connectivity.append([j, len(p) - 1])
                    # update the list of connected atoms
                    connected.append(i)
                    connected.append(j)

        # Create a list of unbound atoms
        unconnected = []
        unconnected_points = vtkPoints()

        for i in range(n):
            if i not in connected:
                unconnected.append(i)
                unconnected_points.InsertNextPoint(p[i])

        # Create a cell array to store the lines in and add the lines to it
        lines = vtkCellArray()

        for atom in connectivity:
            line = vtkLine()
            line.GetPointIds().SetId(0, atom[0])
            line.GetPointIds().SetId(1, atom[1])
            lines.InsertNextCell(line)

        # Create a polydata to store everything in
        lines_poly_data = vtkPolyData()

        # Add the points to the dataset
        lines_poly_data.SetPoints(points)

        # Add the lines to the dataset
        lines_poly_data.SetLines(lines)

        # Create spheres for unbound atoms
        sphere = vtkSphereSource()
        sphere.SetPhiResolution(21)
        sphere.SetThetaResolution(21)
        sphere.SetRadius(.05)

        unbound_point = vtkPolyData()

        # Set the points we created as polydata
        unbound_point.SetPoints(unconnected_points)
        # Create a vtkUnsignedCharArray container and store the colors in it
        named_colors = vtkNamedColors()
        colors = vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        for atom in connectivity:
            curr_color_index = atom[0]
            curr_color = el_color[curr_color_index]
            try:
                colors.InsertNextTupleValue(named_colors.GetColor3ub(curr_color))
            except AttributeError:
                # For compatibility with new VTK generic data arrays.
                colors.InsertNextTypedTuple(named_colors.GetColor3ub(curr_color))

        # Same for unbound atoms:
        unbound_colors = vtkUnsignedCharArray()
        unbound_colors.SetNumberOfComponents(3)
        for atom in unconnected:
            curr_color_index = atom
            curr_color = el_color[curr_color_index]
            try:
                unbound_colors.InsertNextTupleValue(named_colors.GetColor3ub(curr_color))
            except AttributeError:
                # For compatibility with new VTK generic data arrays.
                unbound_colors.InsertNextTypedTuple(named_colors.GetColor3ub(curr_color))

        # Color the lines. SetScalars() automatically associates the values in the data array passed as parameter to
        # the elements in the same indices of the cell data array on which it is called.
        lines_poly_data.GetCellData().SetScalars(colors)

        # Color the unbound atoms.
        unbound_point.GetCellData().SetScalars(unbound_colors)

        # Create a mapper for lines
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(lines_poly_data)
        # for spheres
        mapper1 = vtkGlyph3DMapper()
        mapper1.SetSourceConnection(sphere.GetOutputPort())
        mapper1.SetInputData(unbound_point)
        # Create an actor for lines
        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetLineWidth(2)
        # for spheres
        actor1 = vtkActor()
        actor1.SetMapper(mapper1)
        actor1.GetProperty().LightingOff()

        # Create a renderer and a window
        for xactor in self.renderer.GetActors():
            self.renderer.RemoveActor(xactor)

        self.renderer.AddActor(actor)
        self.renderer.AddActor(actor1)

        # if initial file is empty, reset the renderer when the molecule appears
        if self.reset_flag and n:

            self.renderer.ResetCamera()
            self.reset_flag = False

        self.window.SetWindowName(description)
        self.window.Render()


def main():
    # load coordinates from xyz file
    if len(argv) > 1:
        filename = argv[1]
    else:
        print('Enter filename or drag and drop your .xyz or Gaussian .log file here:')
        filename = input().strip("""'" """)

    refresher = RefreshWindow(filename)

    refresher.read_xyz()

    # Create an interactor and an observer
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(refresher.window)
    style = vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(style)
    interactor.Initialize()

    interactor.CreateRepeatingTimer(1000)
    # noinspection PyTypeChecker
    interactor.AddObserver('TimerEvent', refresher.read_xyz)
    interactor.Start()


if __name__ == "__main__":
    main()
