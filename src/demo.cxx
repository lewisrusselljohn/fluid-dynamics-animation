#include <stdio.h>
#include <glut.h>
#include <math.h>

// Variables and functions for fluid dynamics
const int N = 64;
const int size = (N + 2) * (N + 2);
#define IX(i, j) ((i) + (N + 2)*(j)) // I was missing a bracket around j. Caused big issues in the field arrays!
#define SWAP(x0, x) {float * tmp=x0;x0=x; x=tmp;}

// u and v are the x and y components of the velocity.
// u_prev and v_prev are the components of the user input force field
static float u[size], v[size], u_prev[size], v_prev[size];
static float dens[size], dens_prev[size];

// For drawing an animation of the density field, in place of dynamic calculation
static float dens_frames[20][size];
int record_progress = 0;
int play_progress = -1;

double visc = 0.0;
double dt = 0.099;
double diff = 0.1; // diffusion rate. diff > 0, then density will spread across the grid cells
float force = 50.0f;
int winx, winy;

int solve_depth = 20;


// 
void set_bnd(int N, int b, float* x)
{
	int i;
	for (i = 1; i <= N; i++)
	{
		x[IX(0, i)] =	  b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
		x[IX(i, 0)] =	  b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	}
	x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N+1)] = 0.5f * (x[IX(1, N+1)] + x[IX(0, N)]);
	x[IX(N+1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N+1, 1)]);
	x[IX(N+1, N+1)] = 0.5f * (x[IX(N, N+1)] + x[IX(N+1, N)]);
}

// adds source field s to field x
void add_source(int N, float* x, float* s, float dt)
{
	int i, size = (N + 2) * (N + 2);
	for (i = 0; i < size; i++) x[i] += dt * s[i];
}

void diffuse(int N, int b, float* x, float* x0, float diff, float dt)
{
	int i, j, k;
	float a = dt * diff * N * N;

	for (k = 0; k < solve_depth; k++)
	{
		for (i = 1; i <= N; i++)
		{
			for (j = 1; j <= N; j++)
			{
				x[IX(i, j)] = (x0[IX(i, j)] + a*(x[IX(i-1, j)] + x[IX(i+1, j)]
												+x[IX(i, j-1)]+x[IX(i,j+1)]))/(1+4*a);
			}
		}
		set_bnd(N, b, x);
	}
}

void advect(int N, int b, float* d, float* d0, float* u, float* v, float dt)
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt * N;
	for (i = 1; i <= N; i++)
	{
		for (j = 1; j <= N; j++)
		{
			x = i - dt0 * u[IX(i, j)]; y = j - dt0 * v[IX(i, j)];
			if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f; i0 = (int)x; i1 = i0 + 1;
			if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; j0 = (int)y; j1 = j0 + 1;
			s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
			d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
							s1*(t0*d0[IX(i1, j0)]+t1*d0[IX(i1, j1)]);

			// from demo
			//x = i - dt0 * u[IX(i, j)]; y = j - dt0 * v[IX(i, j)];
			//if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f; i0 = (int)x; i1 = i0 + 1;
			//if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f; j0 = (int)y; j1 = j0 + 1;
			//s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1;
			//d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
			//	s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
		}
	}
	set_bnd(N, b, d);
}

void dens_step(int N, float* x, float* x0, float* u, float* v, float diff,
	float dt)
{
	add_source(N, x, x0, dt);
	SWAP(x0, x); diffuse(N, 0, x, x0, diff, dt);
	SWAP(x0, x); advect(N, 0, x, x0, u, v, dt);
}

// project is the only additional function needed for the velocity step
void project(int N, float* u, float* v, float* p, float* div)
{
	int i, j, k;
	float h;

	h = 1.0f / N;
	for (i = 1; i <= N; i++)
	{
		for (j = 1; j <= N; j++)
		{
			div[IX(i, j)] = -0.5f * h * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]);

			//demo doesnt use h
			//div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
			p[IX(i, j)] = 0;
		}
	}
	set_bnd(N, 0, div); set_bnd(N, 0, p);

	// linear solve
	for (k = 0; k < solve_depth; k++)
	{
		for (i = 1; i <= N; i++)
		{
			for (j = 1; j <= N; j++)
			{
				p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] + p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
			}
		}
		set_bnd(N, 0, p);
	}

	for (i = 1; i <= N; i++)
	{
		for (j = 1; j <= N; j++)
		{
			u[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
			v[IX(i, j)] -= 0.5 * (p[IX(i, j+1)] -p[IX(i, j-1)]) / h;

			// from demo
			//u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
			//v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
		}
	}
	set_bnd(N, 1, u); set_bnd(N, 2, v);
}

void vel_step(int N, float* u, float* v, float* u0, float* v0, float visc, float dt)
{
	add_source(N, u, u0, dt); add_source(N, v, v0, dt);
	SWAP(u0, u); diffuse(N, 1, u, u0, visc, dt);
	SWAP(v0, v); diffuse(N, 2, v, v0, visc, dt);
	project(N, u, v, u0, v0);
	SWAP(u0, u); SWAP(v0, v);
	advect(N, 1, u, u0, u0, v0, dt); advect(N, 2, v, v0, u0, v0, dt);
	project(N, u, v, u0, v0);
}


// GLUT function declarations
void DisplayCallback();
void ReshapeCallback(int, int);
void TimerCallback(int);
void KeyboardCallback(unsigned char, int, int);
void KeyboardUpCallback(unsigned char, int, int);
void MouseCallback(int, int, int, int);

bool firstDisplay = true;

// User input variables
double angle = 0.0;
bool rotateRight = false;
bool rotateLeft = false;
bool go = false;
bool draw_toggle = false;
bool select_frame = false;

static int mouse_down[3];


void draw_velocity()
{
	int i, j;
	float x, y, h;

	h = 1.0f / N;

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

	for (i = 1; i <= N; i++) {
		x = (i - 0.5f) * h;
		for (j = 1; j <= N; j++) {
			y = (j - 0.5f) * h;

			glVertex2f(x, y);
			glVertex2f(x + u[IX(i, j)], y + v[IX(i, j)]);
		}
	}

	glEnd();
}

void draw_density()
{
	int i, j;
	float x, y, h, d00, d01, d10, d11;
	float f00, f01, f10, f11;
	float u00, u01, u10, u11;


	h = 1.0 / N;
	glBegin(GL_QUADS);
	for (i = 0; i <= N; i++) {
		x = (i - 0.5f) * h;
		for (j = 0; j <= N; j++) {
			y = (j - 0.5f) * h;

			d00 = dens[IX(i, j)];
			d01 = dens[IX(i, j + 1)];
			d10 = dens[IX(i + 1, j)];
			d11 = dens[IX(i + 1, j + 1)];

			f00 = v[IX(i, j)];
			f01 = v[IX(i, j + 1)];
			f10 = v[IX(i + 1, j)];
			f11 = v[IX(i + 1, j + 1)];

			u00 = v[IX(i, j)];
			u01 = v[IX(i, j + 1)];
			u10 = v[IX(i + 1, j)];
			u11 = v[IX(i + 1, j + 1)];
			
			//force 
			float f = 1.2;
			glColor4f(d00 + f00 * f, 0.3f * d00 + f00 * f, 0.011f * d00 + f00 * 0.48, 0.0f); glVertex3f(x, y, 0.0);
			glColor4f(d10 + f10 * f, 0.3f * d10 + f10 * f, 0.011f * d10 + f10 * 0.48, 0.0f); glVertex3f(x + h, y, 0.0);
			glColor4f(d11 + f11 * f, 0.3f * d11 + f11 * f, 0.011f * d11 + f11 * 0.48, 0.0f); glVertex3f(x + h, y + h, 0.0);
			glColor4f(d01 + f01 * f, 0.3f * d01 + f01 * f, 0.01f * d01 + f01 * 0.48, 0.0f); glVertex3f(x, y + h, 0.0);
		}
	}
	glEnd();
}

void draw_dens_frames()
{
	int i, j;
	float x, y, h, d00, d01, d10, d11;
	float f00, f01, f10, f11;
	float u00, u01, u10, u11;


	play_progress++;

	h = 1.0 / N;
	glBegin(GL_QUADS);
	for (i = 0; i <= N; i++) {
		x = (i - 0.5f) * h;
		for (j = 0; j <= N; j++) {
			y = (j - 0.5f) * h;

			d00 = dens_frames[play_progress][IX(i, j)];
			d01 = dens_frames[play_progress][IX(i, j + 1)];
			d10 = dens_frames[play_progress][IX(i + 1, j)];
			d11 = dens_frames[play_progress][IX(i + 1, j + 1)];

			//force 
			float f = 1.2;
			glColor4f(d00, 0.3f * d00, 0.011f * d00, 0.0f); glVertex3f(x, y, 0.0);
			glColor4f(d10, 0.3f * d10, 0.011f * d10, 0.0f); glVertex3f(x + h, y, 0.0);
			glColor4f(d11, 0.3f * d11, 0.011f * d11, 0.0f); glVertex3f(x + h, y + h, 0.0);
			glColor4f(d01, 0.3f * d01, 0.01f * d01, 0.0f); glVertex3f(x, y + h, 0.0);
		}
	}
	glEnd();

	if (play_progress == 19)
		play_progress = -1;
}

int main(int argc, char** argv)
{
	// Glut setup
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(720, 486);
	glutInitWindowPosition(200, 100);
	glutCreateWindow("Fluid Dynamics Demo");

	glutDisplayFunc(DisplayCallback);
	glutKeyboardFunc(KeyboardCallback);
	glutKeyboardUpFunc(KeyboardUpCallback);
	glutMouseFunc(MouseCallback);
	glutReshapeFunc(ReshapeCallback);

	glutMainLoop();

	return 0;
}

void DisplayCallback()
{
	if (firstDisplay)
	{
		printf("demo usage: press g to run the dynamic simulation\
			\n s to begin recording the simulation \
			\n d to toggle drawing simulation or animation");

		// Call first frame
		firstDisplay = false;
		glutTimerFunc((unsigned int)1, TimerCallback, glutGet(GLUT_ELAPSED_TIME));
	}
	else
	{
		glViewport(0, 0, winx, winy);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(0.0, 1.0, 0.0, 1.0);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//draw_velocity();
		if (!draw_toggle)
			draw_density();
		else
			draw_dens_frames();
			//draw_velocity();
		glutSwapBuffers();
	}
}

void ReshapeCallback(int pWidth, int pHeight)
{
	winx = pWidth; winy = pHeight;
	glViewport(0, 0, pWidth, pHeight);
}

void TimerCallback(int pPrevious)
{

	// Simulation steps
	// First, a makeshift "get_from_ui" section
	for (int i = 0; i < size; i++)
	{
		u_prev[i] = v_prev[i] = dens_prev[i] = 0.0f;
		//dens_prev[i] = 0.0f;
	}

	if (go)
	{

		v_prev[IX(32, 12)] = force;
		dens_prev[IX(32, 18)] = 1000.0f;
	}

	// solver steps
	vel_step(N, u, v, u_prev, v_prev, visc, dt);
	dens_step(N, dens, dens_prev, u, v, diff, dt);

	/*
		record the density into frames. 
		the last frame is set to the first to avoid a flicker.
		A simple method that probably works best with fire animations since
		the flames of a fire will naturally flick very quickly.

		Density field in 20 frames will take about 85kb, when N = 64
		add 2x that if you want to record the components of velocity field as well.
	*/
	if (go && select_frame && record_progress != 19)
	{
		for (int i = 0; i < size; i++)
		{
			dens_frames[record_progress][i] = dens[i];
		}
		record_progress++;

		if (record_progress == 19)
		{
			for (int i = 0; i < size; i++)
			{
				//dens_frames[record_progress][i] = (dens_frames[record_progress][i] + dens_frames[0][i]) / 2.0f;
				dens_frames[record_progress][i] = dens_frames[0][i];
			}
		}
	}

	// Queue a new frame to be drawn
	glutPostRedisplay();
	glutTimerFunc((unsigned int)1, TimerCallback, glutGet(GLUT_ELAPSED_TIME));
}

void KeyboardCallback(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'g':
	{
		go = true;
		return;
	}
	case 'd':
	{
		draw_toggle = !draw_toggle;
		return;
	}
	case 's':
	{
		if (!select_frame)
		{
			select_frame = true;
		}
	}
	case 'f':
	{
		if (record_progress >= 19)
		{
			printf("recording...\n");
			// open a text file and save the output the values of each frame in a readable format
			//fopen("density_frames.txt",)
		}
	}
	default:
		break;
	}
}

void KeyboardUpCallback(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'g':
		go = false;
	default:
		break;
	}
}

void MouseCallback(int pButton, int pState, int x, int y)
{

}

