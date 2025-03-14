from vpython import *
import serial

current_vertices = [
    [1, -1, -1/8],
    [1,  1, -1/8],
    [-1, 1, -1/8],
    [-1,-1, -1/8],
    [1, -1,  1/8],
    [1,  1,  1/8],
    [-1, 1,  1/8],
    [-1,-1,  1/8]
]

# Save a copy of the original vertices to always use for cumulative rotation.
base_vertices = [v[:] for v in current_vertices]

edges = (
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7)
)

def quaternion_rotation(w: float, x: float, y: float, z: float, v: list) -> list:
    # Rotates vector v using the quaternion (w,x,y,z)
    new_v = [0, 0, 0]
    new_v[0] = (1 - 2 * y * y - 2 * z * z) * v[0] + (2 * x * y - 2 * w * z) * v[1] + (2 * x * z + 2 * w * y) * v[2]
    new_v[1] = (2 * x * y + 2 * w * z) * v[0] + (1 - 2 * x * x - 2 * z * z) * v[1] + (2 * y * z - 2 * w * x) * v[2]
    new_v[2] = (2 * x * z - 2 * w * y) * v[0] + (2 * y * z + 2 * w * x) * v[1] + (1 - 2 * x * x - 2 * y * y) * v[2]
    return new_v

edge_curves = []

def create_cube_edges():
    global edge_curves
    for c in edge_curves:
        c.visible = False
        del c
    edge_curves = []
    for edge in edges:
        c = curve(color=color.white, radius=0.02)
        c.append(vector(*current_vertices[edge[0]]))
        c.append(vector(*current_vertices[edge[1]]))
        edge_curves.append(c)

def update_cube_edges():
    for i, edge in enumerate(edges):
        edge_curves[i].modify(0, vector(*current_vertices[edge[0]]))
        edge_curves[i].modify(1, vector(*current_vertices[edge[1]]))

# --- New Code: Faces with Incremental Colors ---
# Define face indices and assign each face a unique color.
face_defs = [
    ([4, 5, 6, 7], color.red),      # front face (z = 1/8)
    ([0, 1, 2, 3], color.green),    # back face (z = -1/8)
    ([0, 1, 5, 4], color.blue),     # right face (x = 1)
    ([2, 3, 7, 6], color.yellow),   # left face (x = -1)
    ([1, 2, 6, 5], color.cyan),     # top face (y = 1)
    ([0, 3, 7, 4], color.magenta)   # bottom face (y = -1)
]

face_quads = []

def create_cube_faces():
    global face_quads
    for f in face_quads:
        f.visible = False
        del f
    face_quads = []
    for face_indices, face_color in face_defs:
        quad_face = quad(
            v0=vertex(pos=vector(*current_vertices[face_indices[0]]), color=face_color),
            v1=vertex(pos=vector(*current_vertices[face_indices[1]]), color=face_color),
            v2=vertex(pos=vector(*current_vertices[face_indices[2]]), color=face_color),
            v3=vertex(pos=vector(*current_vertices[face_indices[3]]), color=face_color)
        )
        face_quads.append(quad_face)

def update_cube_faces():
    for i, (face_indices, _) in enumerate(face_defs):
        face_quads[i].v0.pos = vector(*current_vertices[face_indices[0]])
        face_quads[i].v1.pos = vector(*current_vertices[face_indices[1]])
        face_quads[i].v2.pos = vector(*current_vertices[face_indices[2]])
        face_quads[i].v3.pos = vector(*current_vertices[face_indices[3]])
# --- End New Code ---

try:
    ser = serial.Serial('COM4', 115200)
    print("Connected to serial port:", ser.name)
except Exception as e:
    print("Error opening serial port:", e)
    exit(1)

create_cube_edges()
create_cube_faces()

while True:
    if ser.in_waiting:
        try:
            line = ser.readline().decode('utf-8', errors="ignore").strip()
            values = line.split(',')
            if len(values) == 4:
                qw, qx, qy, qz = map(float, values)
                print(f"Quaternion: w={qw:.4f}, x={qx:.4f}, y={qy:.4f}, z={qz:.4f}")
                new_vertices = []
                # Apply the cumulative quaternion to the original (base) vertices.
                for v in base_vertices:
                    new_v = quaternion_rotation(qw, qx, qy, qz, v)
                    new_vertices.append(new_v)
                current_vertices = new_vertices
                update_cube_edges()
                update_cube_faces()  # Update the face positions as well.
        except Exception as e:
            print("Error processing data:", e)
