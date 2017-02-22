//
// Created by lurker on 2/13/17.
//

#include "view.h"

void view::_cube(point& p, scalar_t r) {
    GLfloat front  = (GLfloat) (p.data[0] - r);
    GLfloat back   = (GLfloat) (p.data[0] + r);
    GLfloat left   = (GLfloat) (p.data[1] - r);
    GLfloat right  = (GLfloat) (p.data[1] + r);
    GLfloat bottom = (GLfloat) (p.data[2] - r);
    GLfloat top    = (GLfloat) (p.data[2] + r);

    glPushMatrix();

    glBegin(GL_POLYGON);
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, bottom, front);
    glVertex3f(right, top,    front);
    glVertex3f(left,  top,    front);
    glVertex3f(left,  bottom, front);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(right, bottom, back);
    glVertex3f(right, top,    back);
    glVertex3f(left,  top,    back);
    glVertex3f(left,  bottom, back);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, top,    front);
    glVertex3f(right, top,    back);
    glVertex3f(right, bottom, back);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(1.0, 1.0, 0.0);
    glVertex3f(left, bottom, back);
    glVertex3f(left, top,    back);
    glVertex3f(left, top,    front);
    glVertex3f(left, bottom, front);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(0.0, 1.0, 1.0);
    glVertex3f(right, top, back);
    glVertex3f(right, top, front);
    glVertex3f(left,  top, front);
    glVertex3f(left,  top, back);
    glEnd();

    glBegin(GL_POLYGON);
    glColor3f(1.0, 0.0, 1.0);
    glVertex3f(right, bottom, front);
    glVertex3f(right, bottom, back);
    glVertex3f(left,  bottom, back);
    glVertex3f(left,  bottom, front);
    glEnd();

    glPopMatrix();
}

void view::_key(unsigned char k, int x, int y) {
    if (k == 'd') rotate_y += rotate_step;
    else if (k == 'a') rotate_y -= rotate_step;
    else if (k == 'w') rotate_x += rotate_step;
    else if (k == 's') rotate_x -= rotate_step;
    else if (k == 'q') zoom *= (1.0 + zoom_step);
    else if (k == 'e') zoom *= (1.0 - zoom_step);

    glutPostRedisplay();
}

void view::_display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    glRotatef((GLfloat) rotate_x, 1.0, 0.0, 0.0);
    glRotatef((GLfloat) rotate_y, 0.0, 1.0, 0.0);
    glScaled(zoom, zoom, zoom);

    /*
     * draw the tree.
     */
    if (tag == 0) {
        _tree();
    }
    else if (tag == 1){
        _grid();
    }
    else {
        _levelset();
    }

    glFlush();
    glutSwapBuffers();

}



void view::_tree() {

    int count = 0, count_o = 0, count_i=0;

    if (!cached) {
        for (int i = 0; i < treePtr->_dict.size(); ++i) {
            if (treePtr->_dict[i].onBoundary && treePtr->_dict[i].isLeaf) {
                if (treePtr->_dict[i].value == 1) count_o++;
                if (treePtr->_dict[i].value == -1) count_i++;
                if (treePtr->_dict[i].value == 0) count++;

                if (treePtr->_dict[i].value == 1 || treePtr->_dict[i].value == -1) {
                    cachePoints.push_back(treePtr->_dict[i].center - treePtr->_dict[0].center);
                    cacheRadius.push_back(treePtr->_dict[i].radius);
                }
            }
        }
        cached = true;
        std::cout << "boundary size : " << cachePoints.size() << " on " << count<< " out " << count_o << " in " << count_i<< std::endl;
    }
    else {
        for (index_t i = 0; i < cachePoints.size(); ++i) {
            if (cachePoints[i].data[0] >= 0) {
                _cube(cachePoints[i], cacheRadius[i]);
            }

        }
    }
}

void view::_grid() {

    int iX, iY, iZ;
    int len = grid->N;
    int count = 0;

    if (!cached) {
        for (iX =1; iX < len-1; ++iX) {
            for (iY = 1; iY < len-1; ++iY) {
                for (iZ = 1; iZ < len-1; ++iZ) {
                    if (grid->matrix[iX][iY][iZ].phi == 0) {
                        count++;
                        point P = {static_cast<scalar_t >(grid->matrix[iX][iY][iZ].point.x),
                                   static_cast<scalar_t >(grid->matrix[iX][iY][iZ].point.y),
                                   static_cast<scalar_t >(grid->matrix[iX][iY][iZ].point.z)};
                        scalar_t r = grid->stepSize;
                        cachePoints.push_back(P);
                        cacheRadius.push_back(r);
                    }
                }
            }
        }
        cached = true;
        std::cout << "boundary size : " << cachePoints.size() << std::endl;
    }
    else {
        for (index_t i = 0; i < cachePoints.size(); ++i) {
//            if (cachePoints[i].data[0] > -20)
            _cube(cachePoints[i], cacheRadius[i]);
        }
    }

}

void view::_levelset() {
    int iX, iY, iZ;
    int lx = ls->Nx;
    int ly = ls->Ny;
    int lz = ls->Nz;
    double dx = ls->dx;
    double sx = ls->sx;
    double sy = ls->sy;
    double sz = ls->sz;
    int count = 0;

    if (!cached) {
        for (iX =1; iX < lx-1; ++iX) {
            for (iY = 1; iY < ly-1; ++iY) {
                for (iZ = 1; iZ < lz-1; ++iZ) {
                    if (fabs(ls->get(ls->phi, iX, iY, iZ)) < ls->dx) {
                        count++;
                        point P = {
                                sx + iX * dx,
                                sy + iY * dx,
                                sz + iZ * dx
                        };
                        scalar_t r = dx / 2.0;
                        cachePoints.push_back(P);
                        cacheRadius.push_back(r);
                    }
                }
            }
        }
        cached = true;
        std::cout << "boundary size : " << cachePoints.size() << std::endl;
    }
    else {
        for (index_t i = 0; i < cachePoints.size(); ++i) {
//            if (cachePoints[i].data[0] > -20)
            _cube(cachePoints[i], cacheRadius[i]);
        }
    }


}


void view::run() {
    int argc = 1;
    char* argv[1]; argv[0] = strdup(TITLE);
    glutInit(&argc, argv);

    glutInitWindowSize(WIDTH, HEIGHT);
    glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - WIDTH)/2,
                           (glutGet(GLUT_SCREEN_HEIGHT) - HEIGHT)/2);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    glutCreateWindow(TITLE);
    glEnable(GL_DEPTH_TEST);

    glutKeyboardFunc(view::key_callback);
    glutDisplayFunc(view::display_callback);
    glutReshapeFunc(view::reshape_callback);

    glutMainLoop();

    return;
}