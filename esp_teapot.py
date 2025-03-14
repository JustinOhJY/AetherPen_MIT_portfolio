import sys
import math
import serial
import glfw
from OpenGL.GL import *
from OpenGL.GLU import *
from qfuncs import *

vertices = [
    [1, -1, -1/8],
    [1, 1, -1/8],
    [-1, 1, -1/8],
    [-1, -1, -1/8],
    [1, -1, 1/8],
    [1, 1, 1/8],
    [-1, 1, 1/8],
    [-1, -1, 1/8]
]

vertices1 = (
    (1, -1, -1/8),
    (1, 1, -1/8),
    (-1, 1, -1/8),
    (-1, -1, -1/8),
    (1, -1, 1/8),
    (1, 1, 1/8),
    (-1, 1, 1/8),
    (-1, -1, 1/8)
)

edges = (
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7)
)

surfaces = (
    (0, 1, 2, 3),  # front
    (3, 2, 6, 7),  # left
    (6, 5, 4, 7),  # back
    (4, 5, 1, 0),  # right
    (1, 5, 6, 2),  # top
    (4, 0, 3, 7)   # bottom
)

colors = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (0, 1, 1),
    (1, 0, 1),
    (1, 1, 0)
)


def initialize_opengl(width: int = 800, height: int = 600):
    glEnable(GL_DEPTH_TEST)
    glDepthFunc(GL_LEQUAL)
    glClearDepth(1.0)
    glClearColor(0.0, 0.0, 0.0, 1.0)
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    gluPerspective(45, width / height, 0.1, 50.0)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    glTranslatef(0.0, 0.0, -5)


def draw_cube():
    """
    Draws a cube using the global vertices, surfaces, edges, and colors.
    """
    glBegin(GL_QUADS)
    for surface in surfaces:
        for i, vertex in enumerate(surface):
            glColor3fv(colors[i % len(colors)])
            glVertex3fv(vertices[vertex])
    glEnd()

    glBegin(GL_LINES)
    glColor3fv((1, 1, 1))
    for edge in edges:
        for vertex in edge:
            glVertex3fv(vertices[vertex])
    glEnd()


def main():

    if not glfw.init():
        print("Failed to initialize GLFW")
        sys.exit(1)

    window = glfw.create_window(800, 600, "OpenGL Cube", None, None)
    if not window:
        glfw.terminate()
        print("Failed to create GLFW window")
        sys.exit(1)

    glfw.make_context_current(window)
    initialize_opengl()

    # Initialize serial connection
    try:
        ser = serial.Serial('COM4', 115200)
        print(f"Connected to serial port: {ser.name}")
    except serial.SerialException as e:
        print(f"Serial connection error: {e}")
        glfw.terminate()
        sys.exit(1)

    while not glfw.window_should_close(window):
        glfw.poll_events()

        if ser.in_waiting:
            line = ser.readline().decode('utf-8', errors="ignore").strip().split(",")
            if len(line) == 4:
                try:
                    qw, qx, qy, qz = map(float, line)
                    print(f"w={qw:.4f}, x={qx:.4f}, y={qy:.4f}, z={qz:.4f}")
                    # Apply quaternion rotation to each vertex of the cube
                    for index, v in enumerate(vertices1):
                        vertices[index] = quaternion_rotation1(qw, qx, qy, qz, v)
                except ValueError:
                    print("Error: Invalid data received; could not convert to float.")

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        draw_cube()
        glfw.swap_buffers(window)
    glfw.terminate()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Exiting")
