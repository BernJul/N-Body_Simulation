#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <tuple>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"


int n_body;
int n_iteration;

int n_omp_threads;


void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % (max_mass - min_mass) + min_mass;
        x[i] = rand() % bound_x;
        y[i] = rand() % bound_y;
        
        // For collision logic (fast paced)
        // x[i] = 2000.0f + rand() % (bound_x / 4);
        // y[i] = 2000.0f + rand() % (bound_y / 4);

        // x[i] = 3000.0f + rand() % (bound_x / 4);
        // y[i] = 3000.0f + rand() % (bound_y / 4);

        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}

std::tuple<double, double> get_accel(double m, double x_1, double x_2, double y_1, double y_2){
    // Calculate the acceleration
    double deltaX = x_2 - x_1;
    double deltaY = y_2 - y_1;
    
    double radius = (deltaX * deltaX) + (deltaY * deltaY) + err;
    radius = sqrt(radius);

    radius = (radius <= radius2) ? radius2 : radius;

    double acceleration = (gravity_const * m) / (radius * radius);

    // Collision with other bodies
    if (radius == radius2 || radius <= radius2)

        return std::make_tuple(0, 0);

    // Get acceleration of x and y direction
    double ax = deltaX/radius * acceleration;
    double ay = deltaY/radius * acceleration;

    return std::make_tuple(ax, ay);
}

void update_position(double *x, double *y, double *vx, double *vy, int n) {
    // Update position 
    for (size_t i = 0; i < n; i++) {
        x[i] += (vx[i] * dt);
        y[i] += (vy[i] * dt);

        // Border collision logic
        if (x[i] < 0) vx[i] *= -1;
        if (y[i] < 0) vy[i] *= -1;

        if (x[i] >= bound_x) vx[i] *= -1;
        if (y[i] >= bound_y) vy[i] *= -1;

    }

}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    // Calculate force and acceleration, update velocity
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) continue;
            double ax, ay;
            std::tie(ax, ay) = get_accel(m[j], x[i], x[j], y[i], y[j]);

            vx[i] = vx[i] + (dt * ax);
            vy[i] = vy[i] + (dt * ay);
        }
    }
}

void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("OpenMP", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        //Threads configuration
        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_velocity(m, x, y, vx, vy, i);
        }

        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_position(x, y, vx, vy, i);
        }

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;

        printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        l.save_frame(x, y);

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = x[i];
            yi = y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    
}


int main(int argc, char *argv[]){
    
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_omp_threads = atoi(argv[3]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation OpenMP Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

    printf("Student ID: 119010520\n");
    printf("Name: Bernaldy Jullian\n");
    printf("Assignment 2: N Body Simulation OpenMP Implementation\n");
    
    return 0;

}


