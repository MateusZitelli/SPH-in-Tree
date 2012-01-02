//SPH in Tree - SPH using barnes-hut algorithm

//Copyright 2012 Mateus Zitelli <zitellimateus@gmail.com>

//This program is free software; you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation; either version 2 of the License, or
//(at your option) any later version.

//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//MA 02110-1301, USA.
#include <stdlib.h>
#include <GL/glut.h>
#include "helper.h"
#define colorMax 70E3

float angleX = 0;
float angleY = 0;
float positionZ = -10;
int paused = 0, render = 1;

void initRendering()
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void handleResize(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (float)w / (float)h, 1.0, 200.0);
}

void drawScene()
{
	int i;
	float r,g,b;
	double acel;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(0.0f, 0.0f, positionZ);

	glRotatef(angleX, 0.0f, 1.0f, 0.0f);
	glRotatef(angleY, 1.0f, 0.0f, 0.0f);
	glBegin(GL_POINTS);

	for (i = 0; i < BODIES_QUANTITY; i++) {
	        glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		glVertex3f(bodies[i].position(0) / 100.0,bodies[i].position(1) / 100.0,bodies[i].position(2) / 100.0);
		//glutSolidSphere(0.01,10, 10);
	}
	glEnd();
	glBegin(GL_LINE_LOOP);		// Draw The Cube Using quads
                glColor3f(0.0f,1.0f,0.0f);	// Color Blue
                glVertex3f( 1.0f, 3.0f,-1.0f);	// Top Right Of The Quad (Top)
                glVertex3f(-1.0f, 3.0f,-1.0f);	// Top Left Of The Quad (Top)
                glVertex3f(-1.0f, 1.0f, 1.0f);	// Bottom Left Of The Quad (Top)
                glVertex3f( 1.0f, 1.0f, 1.0f);	// Bottom Right Of The Quad (Top)
                glColor3f(1.0f,0.5f,0.0f);	// Color Orange
                glVertex3f( 1.0f,-1.0f, 1.0f);	// Top Right Of The Quad (Bottom)
                glVertex3f(-1.0f,-1.0f, 1.0f);	// Top Left Of The Quad (Bottom)
                glVertex3f(-1.0f,-1.0f,-1.0f);	// Bottom Left Of The Quad (Bottom)
                glVertex3f( 1.0f,-1.0f,-1.0f);	// Bottom Right Of The Quad (Bottom)
                glColor3f(1.0f,0.0f,0.0f);	// Color Red	
                glVertex3f( 1.0f, 3.0f, 1.0f);	// Top Right Of The Quad (Front)
                glVertex3f(-1.0f, 3.0f, 1.0f);	// Top Left Of The Quad (Front)
                glVertex3f(-1.0f,-1.0f, 1.0f);	// Bottom Left Of The Quad (Front)
                glVertex3f( 1.0f,-1.0f, 1.0f);	// Bottom Right Of The Quad (Front)
                glColor3f(1.0f,1.0f,0.0f);	// Color Yellow
                glVertex3f( 1.0f,-1.0f,-1.0f);	// Top Right Of The Quad (Back)
                glVertex3f(-1.0f,-1.0f,-1.0f);	// Top Left Of The Quad (Back)
                glVertex3f(-1.0f, 3.0f,-1.0f);	// Bottom Left Of The Quad (Back)
                glVertex3f( 1.0f, 3.0f,-1.0f);	// Bottom Right Of The Quad (Back)
                glColor3f(0.0f,0.0f,1.0f);	// Color Blue
                glVertex3f(-1.0f, 3.0f, 1.0f);	// Top Right Of The Quad (Left)
                glVertex3f(-1.0f, 3.0f,-1.0f);	// Top Left Of The Quad (Left)
                glVertex3f(-1.0f,-1.0f,-1.0f);	// Bottom Left Of The Quad (Left)
                glVertex3f(-1.0f,-1.0f, 1.0f);	// Bottom Right Of The Quad (Left)
                glColor3f(1.0f,0.0f,1.0f);	// Color Violet
                glVertex3f( 1.0f, 3.0f,-1.0f);	// Top Right Of The Quad (Right)
                glVertex3f( 1.0f, 3.0f, 1.0f);	// Top Left Of The Quad (Right)
                glVertex3f( 1.0f,-1.0f, 1.0f);	// Bottom Left Of The Quad (Right)
                glVertex3f( 1.0f,-1.0f,-1.0f);	// Bottom Right Of The Quad (Right)
	glEnd();
	glutSwapBuffers();
	glFlush();
}

void update(void)
{
	int i, j, k;

	//angleX += 0.5;
	if (angleX > 360)
		angleX -= 360;
	if (angleY > 360)
		angleY -= 360;
	if (angleX < 0)
		angleX += 360;
	if (angleY < 0)
		angleY += 360;
	roots_quantity = 0;
	resetNodes();
	divideNode(&nodes[0]);
	for (i = 0; i < node_quantity; i++) {
		setCenterOfMass(&nodes[i]);
	}
	for (i = 0; i < roots_quantity; i++) {
		for (j = 0; j < roots[i]->bodies_quantity; j++) {
			forceOverNode(roots[i], NULL, roots[i]->bodies[j], 0);
		}
	}
	fprintf(positionData, "FF\n");
	if (render)
		glutPostRedisplay();
	for (i = 0; i < BODIES_QUANTITY; i++) {
	        bodies[i].force(1) -= 0.0001;
	        if(bodies[i].position(0) >= 100) bodies[i].force(0) -= (bodies[i].position(0) - 100) / 30.0;
	        else if(bodies[i].position(0) < -100) bodies[i].force(0) -= (bodies[i].position(0) + 100) / 30.0;
	        if(bodies[i].position(1) >= 300) bodies[i].force(1) -= (bodies[i].position(1) - 300) / 30.0;
	        else if(bodies[i].position(1) < -100) bodies[i].force(1) -= (bodies[i].position(1) + 100) / 30.0;
	        if(bodies[i].position(2) >= 100) bodies[i].force(2) -= (bodies[i].position(2) - 100) / 30.0;
	        else if(bodies[i].position(2) < -100) bodies[i].force(2) -= (bodies[i].position(2) + 100) / 30.0;
		bodies[i].speed(0) += bodies[i].force(0) / bodies[i].mass;
		bodies[i].speed(1) += bodies[i].force(1) / bodies[i].mass;
		bodies[i].speed(2) += bodies[i].force(2) / bodies[i].mass;
		//bodies[i].acel = sqrt(bodies[i].force(0) * bodies[i].force(0) + bodies[i].force(1) * bodies[i].force(1) + bodies[i].force(2) * bodies[i].force(2));
		fprintf(positionData, "%i,%i,%i\n",
			(int)(bodies[i].position(0) * 10E-16),
			(int)(bodies[i].position(1) * 10E-16),
			(int)(bodies[i].position(2) * 10E-16));
		bodies[i].position(0) += bodies[i].speed(0) * 5E1;
		bodies[i].position(1) += bodies[i].speed(1) * 5E1;
		bodies[i].position(2) += bodies[i].speed(2) * 5E1;
		bodies[i].force(0) = 0;
		bodies[i].force(1) = 0;
		bodies[i].force(2) = 0;
	}
}

void handleKeypress(unsigned char key, int x, int y)
{
	switch (key) {
	case 27:		//Escape key
		exit(x);
		break;
	case GLUT_KEY_RIGHT:
		angleX += 2;
		break;
	case GLUT_KEY_LEFT:
		angleX -= 2;
		break;
	case GLUT_KEY_UP:
		angleY += 2;
		break;
	case GLUT_KEY_DOWN:
		angleY -= 2;
		break;
	case 112:
		paused = !paused;
		if (paused) {
			glutIdleFunc(NULL);
			fprintf(positionData, "P1\n");
		} else {
			glutIdleFunc(update);
			fprintf(positionData, "P0\n");
		}
		break;
	case 114:
		render = !render;
		break;
	case 43:
		positionZ -= 1;
		break;
	case 95:
		positionZ += 1;
		break;
	}
	printf("%i\n", key);
	glutPostRedisplay();
}

void handleSKeypress(int key, int x, int y)
{
	switch (key) {
	case 27:		//Escape key
		exit(x);
		break;
	case GLUT_KEY_RIGHT:
		angleX += 2;
		break;
	case GLUT_KEY_LEFT:
		angleX -= 2;
		break;
	case GLUT_KEY_UP:
		angleY += 2;
		break;
	case GLUT_KEY_DOWN:
		angleY -= 2;
		break;
	case 112:
		paused = !paused;
		if (paused) {
			glutIdleFunc(NULL);
			fprintf(positionData, "P1\n");
		} else {
			glutIdleFunc(update);
			fprintf(positionData, "P0\n");
		}
		break;
	case 114:
		render = !render;
		break;
	case 43:
		positionZ -= 1;
		break;
	case 95:
		positionZ += 1;
		break;
	}
	printf("%i\n", key);
	glutPostRedisplay();
}

int main(int argc, char **argv)
{
	int i, frame = 0;
	srand(time(0));
	init();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_DEPTH);
	glutInitWindowSize(WIDTH, HEIGHT);

	glutCreateWindow("teste");
	initRendering();

	glutDisplayFunc(drawScene);
	glutKeyboardFunc(handleKeypress);
	glutSpecialFunc(handleSKeypress);
	glutReshapeFunc(handleResize);
	glutIdleFunc(update);

	glutMainLoop();
	return 0;
}
