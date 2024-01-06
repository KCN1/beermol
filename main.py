# main process does not stop if window is closed while running read_xyz()
# adjust timer duration depending on the size, add window description
# add a monitor for when the molecule appears

from gc import collect as gc_collect
from sys import argv as sys_argv, exit as sys_exit
from math import dist
from scipy.spatial import KDTree
import vtk
from PyQt5 import QtCore, QtWidgets
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class MoleculeRenderer:

    _atom_numbers = (
        'Bq', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F',
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
    _hyper_colors = {
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
            'Au': 'Yellow',
            'Bq': 'Pink'
        }
    _radii = {
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
            'Md': 1.73, 'No': 1.76, 'Lr': 1.61, 'Bq': 0.0
        }

    def __init__(self, filename):
        self.color_dict = self._hyper_colors
        self.filename = filename
        self.ren = vtk.vtkRenderer()

    def from_orca_gauss(self):

        n, description, p, el = 0, '', [], []
        flag = 'ready'  # start/stop flag for reading the data
        log_format = ''
        geometries = 0  # to store all geometries (for next version)
        n0, p0, el0 = 0, [], []  # current number of atoms, coordinates, elements

        with open(self.filename, 'r') as file:

            for line in file:   # guess the log file format
                if 'Gaussian' in line:
                    log_format = 'gaussian'
                    break
                elif 'O   R   C   A' in line:
                    log_format = 'orca'
                    break
                elif 'Priroda' in line:
                    log_format = 'priroda'
                    break
                # elif 'GAMESS' in line:
                #     log_format = 'gamess'

            line = file.readline()
            if log_format == 'orca':    # try to find the description
                while line and 'END OF INPUT' not in line:
                    if '|  1> #' in line:
                        description = line.lstrip('|  1> #').strip()
                    line = file.readline()
            elif log_format == 'gaussian':
                curr_line, prev_line, target_line = '', '', ''
                while line and 'Charge' not in line and 'orientation' not in line:
                    target_line = prev_line
                    prev_line = curr_line
                    curr_line = line
                    line = file.readline()
                if '-' in prev_line:
                    description = target_line
            elif log_format == 'priroda':
                while line and 'atoms' not in line:
                    if 'molecule input:' in line:
                        description = line.split("'")[1]
                    line = file.readline()
            else:
                description = 'Unknown format'

            for line in file:
                if flag == 'ready':
                    if log_format == 'orca' and 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                        flag = 'set'
                    elif log_format == 'gaussian' and 'Coordinates (Angstroms)' in line:
                        flag = 'set'
                    elif log_format == 'priroda' and 'Atomic Coordinates:' in line:
                        flag = 'go'
                elif flag == 'set' and '-------' in line:
                    flag = 'go'
                elif flag == 'go':
                    if line.strip() and '-------' not in line and '#' not in line:
                        s = line.split()
                        p0.append((float(s[-3]), float(s[-2]), float(s[-1])))
                        if log_format == 'gaussian':
                            el0.append(self._atom_numbers[int(s[1])])
                        elif log_format == 'priroda':
                            el0.append(self._atom_numbers[int(s[0])])
                        elif log_format == 'orca':
                            el0.append(s[0])
                        n0 += 1
                    else:
                        flag = 'ready'
                        n = n0  # last number of atoms
                        p = p0.copy()  # last coordinates
                        el = el0.copy()  # last list of els
                        n0, p0, el0 = 0, [], []
                        geometries += 1

        return n, description, el, p

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

    def read_xyz(self):
        """Creates separate mappers for lines and spheres from a list of coordinates"""

        # read and translate with dictionary
        if self.filename[-4:] == '.xyz':
            n, description, el, p = self.from_trj_xyz()
        elif self.filename[-4:] in ('.log', '.out'):
            n, description, el, p = self.from_orca_gauss()
        else:
            n, description, el, p = 0, '', [], []

        # create a points array for atoms
        points = vtk.vtkPoints()
        for atom in p:
            points.InsertNextPoint(atom)

        # define element colors
        el_color = ['White'] * n
        for i in range(n):
            if el[i] in self.color_dict:
                el_color[i] = self.color_dict[el[i]]

        # build a KD-Tree to search for pairs of atoms within max_bond distance (sum of max(radii) + 10%)
        radii_incr = [1.1 * self._radii[atom] for atom in el]

        if n:
            max_bond = 2 * max(radii_incr)
            kd_tree = KDTree(p, leafsize=200)
            prelim_pairs = kd_tree.query_pairs(r=max_bond)
        else:
            prelim_pairs = []

        # create a connectivity array
        connectivity = []
        unconnected = set(range(n))

        for (i, j) in prelim_pairs:
            if dist(p[i], p[j]) < radii_incr[i] + radii_incr[j]:
                # add virtual points in the middle of each bond
                pv = ((p[i][0] + p[j][0]) / 2, (p[i][1] + p[j][1]) / 2, (p[i][2] + p[j][2]) / 2)
                p.append(pv)
                points.InsertNextPoint(pv)
                # connect last pair of bonded atoms to the last virtual point
                connectivity.append((i, len(p) - 1))
                connectivity.append((j, len(p) - 1))
                # update the list of unbound atoms
                unconnected.discard(i)
                unconnected.discard(j)

        # Create an array of unbound atoms
        unconnected_points = vtk.vtkPoints()

        for i in unconnected:
            unconnected_points.InsertNextPoint(p[i])

        # Create a cell array to store the lines
        lines = vtk.vtkCellArray()

        for (i, j) in connectivity:
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, i)
            line.GetPointIds().SetId(1, j)
            lines.InsertNextCell(line)

        # Create a polydata to store points and lines in
        lines_poly_data = vtk.vtkPolyData()
        # Add the points to the dataset
        lines_poly_data.SetPoints(points)
        # Add the lines to the dataset
        lines_poly_data.SetLines(lines)
        # Create spheres for unbound atoms
        sphere = vtk.vtkSphereSource()
        sphere.SetPhiResolution(21)
        sphere.SetThetaResolution(21)
        sphere.SetRadius(.05)
        # Create a polydata to store unbound points in
        unbound_point = vtk.vtkPolyData()
        # Set the points we created as polydata
        unbound_point.SetPoints(unconnected_points)

        # Create a vtkUnsignedCharArray container and store the colors in it
        named_colors = vtk.vtkNamedColors()
        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)

        for (i, j) in connectivity:
            curr_color = el_color[i]
            try:
                colors.InsertNextTupleValue(named_colors.GetColor3ub(curr_color))
            except AttributeError:
                # For compatibility with new VTK generic data arrays.
                colors.InsertNextTypedTuple(named_colors.GetColor3ub(curr_color))

        # Same for unbound atoms:
        unbound_colors = vtk.vtkUnsignedCharArray()
        unbound_colors.SetNumberOfComponents(3)

        for i in unconnected:
            curr_color = el_color[i]
            try:
                unbound_colors.InsertNextTupleValue(named_colors.GetColor3ub(curr_color))
            except AttributeError:
                # For compatibility with new VTK generic data arrays.
                unbound_colors.InsertNextTypedTuple(named_colors.GetColor3ub(curr_color))

        # Color the lines. SetScalars() automatically associates the values in the data array passed as parameter
        # to the elements in the same indices of the cell data array on which it is called.
        lines_poly_data.GetCellData().SetScalars(colors)
        # Color the unbound atoms.
        unbound_point.GetCellData().SetScalars(unbound_colors)

        # Create mappers and actors for lines and spheres
        mapper_lines = vtk.vtkPolyDataMapper()
        mapper_lines.SetInputData(lines_poly_data)

        mapper_spheres = vtk.vtkGlyph3DMapper()
        mapper_spheres.SetSourceConnection(sphere.GetOutputPort())
        mapper_spheres.SetInputData(unbound_point)

        actor_lines = vtk.vtkActor()
        actor_spheres = vtk.vtkActor()

        actor_lines.SetMapper(mapper_lines)
        actor_lines.GetProperty().SetLineWidth(2)

        actor_spheres.SetMapper(mapper_spheres)
        actor_spheres.GetProperty().LightingOff()

        for old_actor in self.ren.GetActors():
            self.ren.RemoveActor(old_actor)

        self.ren.AddActor(actor_lines)
        self.ren.AddActor(actor_spheres)
        gc_collect()


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, parent=None, molecule=None):
        QtWidgets.QMainWindow.__init__(self, parent)

        self.molecule = molecule
        self.frame = QtWidgets.QFrame()

        self.vl = QtWidgets.QVBoxLayout()
        self.vtkWidget = MyQVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget)

        self.window = self.vtkWidget.GetRenderWindow()
        self.window.AddRenderer(self.molecule.ren)
        self.iren = self.window.GetInteractor()

        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)
        self.show()
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        self.iren.Initialize()

        self.repeat_time = 10000
        self.reset_flag = True
        self.refresh()  # here to add a monitor for when the molecule appears
        self.vtkWidget.CreateRepeatingTimer(self.repeat_time)
        self.vtkWidget.AddObserver('TimerEvent', self.refresh)

    def refresh(self, caller=None, event=None):
        self.molecule.read_xyz()
        if self.reset_flag:
            self.molecule.ren.ResetCamera()
            self.reset_flag = False
        self.window.Render()


class MyQVTKRenderWindowInteractor(QVTKRenderWindowInteractor):
    _TimerDuration = 100

    def CreateTimer(self, obj, event):
        self._Timer.start(self._TimerDuration)  # self._Timer.start(10) in original

    def CreateRepeatingTimer(self, duration):
        self._TimerDuration = duration
        super().GetRenderWindow().GetInteractor().CreateRepeatingTimer(duration)


def main():

    app = QtWidgets.QApplication(sys_argv)

    if len(sys_argv) > 1:
        filename = sys_argv[1]
    else:
        print('Enter a filename or drag and drop your .xyz or Gaussian .log file here:')
        filename = input().strip("""'" """)

    mol_renderer = MoleculeRenderer(filename)
    window = MainWindow(molecule=mol_renderer)
    window.resize(800, 600)

    sys_exit(app.exec_())


if __name__ == "__main__":
    main()
