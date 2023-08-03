#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <pthread.h>
#include <tuple>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"

int n_thd; // number of threads

int n_body;
int n_iteration;
int thread_tracking;

double *m, *x, *y, *vx, *vy;

void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % (max_mass - min_mass) + min_mass;
        // x[i] = rand() % bound_x;
        // y[i] = rand() % bound_y;
        
        // collision logic
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);

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

void update_position(double *x, double *y, double *vx, double *vy, int start, int end) {
    // Update position 
    for (size_t i = start; i < end; i++) {
        x[i] += (vx[i] * dt);
        y[i] += (vy[i] * dt);

        // Border collision logic
        if (x[i] < 0) vx[i] *= -1;
        if (y[i] < 0) vy[i] *= -1;
        if (x[i] >= bound_x) vx[i] *= -1;
        if (y[i] >= bound_y) vy[i] *= -1;

    }

}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int start, int end, int n) {
    // Calculate force and acceleration, update velocity
    for (size_t i = start; i < end; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) continue;
            double ax, ay;
            std::tie(ax, ay) = get_accel(m[j], x[i], x[j], y[i], y[j]);

            vx[i] = vx[i] + (dt * ax);
            vy[i] = vy[i] + (dt * ay);
        }
    }
}
typedef struct {
    // Thread arguments
    int id;
    int start;
    int end;
} Args;

void* worker(void* args) {
    // Procedure in each threads

    Args* my_arg = (Args*) args;
    int id = my_arg->id;
    int start = my_arg->start;
    int end = my_arg->end;

    update_velocity(m, x, y, vx, vy, start, end, n_body);
    update_position(x, y, vx, vy, start, end);

    pthread_exit(NULL);
}

void master(){
    m = new double[n_body];
    x = new double[n_body];
    y = new double[n_body];
    vx = new double[n_body];
    vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("Pthread", n_body, bound_x, bound_y);

    pthread_t threads[n_thd];
    Args args[n_thd];

    int remainder = n_body % n_thd; // remaider of data

    /* Pass on arguments start_addr and end_addr */
    int temp = 0;
    for (size_t i = 0; i < n_thd; i++) {
        args[i].id = i;
        args[i].start = temp;
        int num_my_elements = (i < remainder) ? (n_body / n_thd) + 1 : (n_body / n_thd);
        temp += num_my_elements;
        args[i].end = temp - 1;
    }
    
    // Generate threads and collect from them for every frame
    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        // Assign jobs
        for (size_t i = 0; i < n_thd; i++) pthread_create(&threads[i], NULL, worker, &args[i]);
        
        // Wait for jobs to finish
        for (int i = 0; i < n_thd; i++) pthread_join(threads[i], NULL);

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


int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_thd = atoi(argv[3]);

    #ifdef GUI
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Pthread");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

    printf("Student ID: 119010520\n");
    printf("Name: Bernaldy Jullian\n");
    printf("Assignment 2: N Body Simulation Pthread Implementation\n");

	return 0;
}
