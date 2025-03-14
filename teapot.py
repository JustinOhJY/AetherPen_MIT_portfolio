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

try:
    ser = serial.Serial('COM4', 115200)
    print("Connected to serial port:", ser.name)
except Exception as e:
    print("Error opening serial port:", e)
    exit(1)

create_cube_edges()

while True:
    if ser.in_waiting:
        try:
            line = ser.readline().decode('utf-8', errors="ignore").strip()
            values = line.split(',')
            if len(values) == 4:
                qw, qx, qy, qz = map(float, values)
                print(f"Quaternion: w={qw:.4f}, x={qx:.4f}, y={qy:.4f}, z={qz:.4f}")
                new_vertices = []
                for v in current_vertices:
                    new_v = quaternion_rotation(qw, qx, qy, qz, v)
                    new_vertices.append(new_v)
                current_vertices = new_vertices
                update_cube_edges()
        except Exception as e:
            print("Error processing data:", e)
