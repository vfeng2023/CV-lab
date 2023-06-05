import stl
from stl import mesh

your_mesh = mesh.Mesh.from_file('ioscahedron.stl')
your_mesh.save('OUTPUT.stl',mode=stl.Mode.ASCII)