#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <limits>
#include <queue>

float MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1)
{
    int sign = 0;

    FaceHandle fh0 = _mesh->face_handle ( faceID0 );
    FaceHandle fh1 = _mesh->face_handle ( faceID1 );

    // On récupère les normales, calcule leur produit scalaire et on renvoie l'angle
    // non signé avec le produit scalaire
    OpenMesh::Vec3f normal0 (_mesh->normal ( fh0 ) );
    OpenMesh::Vec3f normal1 (_mesh->normal ( _mesh->face_handle ( faceID1 ) ) );
    float scalar = normal0 | normal1;

    return abs ( acos ( scalar ) );
}

void MainWindow::computeAngularDistances( MyMesh *_mesh ) {
    for ( MyMesh::FaceIter curFace = _mesh->faces_begin( ) ; curFace != _mesh->faces_end( ) ; curFace++ ) {
        for ( MyMesh::FaceFaceIter curNeigh = _mesh->ff_iter( *curFace ) ; curNeigh.is_valid() ; curNeigh++ ) {

            FaceHandle fh0 = *curFace;
            FaceHandle fh1 = *curNeigh;
            float angle = angleFF ( _mesh , fh0.idx() , fh1.idx() );
            float coef = angle > M_PI ? 0.1 : 1;
            angularDistances[std::make_pair( fh0.idx() , fh1.idx() )] = coef * ( 1 - cos( angle ) );

        }
    }
}

MyMesh::Point MainWindow::faceGravityCenter( MyMesh *_mesh , int faceID ) {
    FaceHandle fh = _mesh->face_handle( faceID );
    std::vector<VertexHandle> vertices;
    for( MyMesh::FaceVertexIter curVert = _mesh->fv_iter( fh ) ; curVert.is_valid() ; curVert++ ) {
        vertices.push_back( *curVert );
    }

    MyMesh::Point result = ( _mesh->point( vertices[0] ) + _mesh->point( vertices[1] ) + _mesh->point( vertices[2] ) ) / 3;
    return result;
}

std::vector<MyMesh::Point> MainWindow::discretizeEdge( MyMesh *_mesh , int edgeID ) {
    EdgeHandle eh = _mesh->edge_handle( edgeID );
    HalfedgeHandle heh = _mesh->halfedge_handle( eh , 0 );

    VertexHandle vh0 = _mesh->to_vertex_handle( heh );
    VertexHandle vh1 = _mesh->from_vertex_handle( heh );

    std::vector<MyMesh::Point> discretizedPoints {_mesh->point( vh0 ) , _mesh->point( vh1 )};
    std::vector<MyMesh::Point> updatedPoints ( discretizedPoints );

    while( discretizedPoints.size() < 17 ) {

        for ( int i = 0 ; i < discretizedPoints.size() ; ++i ) {

            MyMesh::Point newPoint = (discretizedPoints[i] + discretizedPoints[i + 1]) / 2;
            auto it = updatedPoints.begin();
            updatedPoints.insert( it + i + 1 , newPoint );

        }
        discretizedPoints = updatedPoints;

    }

    return discretizedPoints;
}

float MainWindow::geodesicDistance( MyMesh *_mesh , int edgeID ) {
    EdgeHandle eh = _mesh->edge_handle( edgeID );
    HalfedgeHandle heh0 = _mesh->halfedge_handle( eh , 0 );
    HalfedgeHandle heh1 = _mesh->halfedge_handle( eh , 1 );

    FaceHandle fh0 = _mesh->face_handle ( heh0 );
    FaceHandle fh1 = _mesh->face_handle ( heh1 );

    MyMesh::Point center0 = faceGravityCenter( _mesh , fh0.idx() );
    MyMesh::Point center1 = faceGravityCenter( _mesh , fh1.idx() );

     std::vector<MyMesh::Point> edgePoints = discretizeEdge( _mesh , edgeID );

     std::sort ( edgePoints.begin ( ) , edgePoints.end ( ) , [=](MyMesh::Point p1, MyMesh::Point p2) {
         VectorT<float, 3> vector0 = p1 - center0;
         VectorT<float, 3> vector1 = center1 - p1;
         float length0 = abs( vector0.norm() ) + abs( vector1.norm() );


         vector0 = p2 - center0;
         vector1 = center1 - p2;
         float length1 = abs( vector0.norm() ) + abs( vector1.norm() );

         return length0 < length1;
     } );

     MyMesh::Point pointForShortestPath = edgePoints.front();
     VectorT<float, 3> vector0 = pointForShortestPath - center0;
     VectorT<float, 3> vector1 = center1 - pointForShortestPath;
     return ( abs( vector0.norm() ) + abs( vector1.norm() ) );
}

void MainWindow::computeGeodesicDistances( MyMesh *_mesh ) {
    for ( MyMesh::FaceIter curFace = _mesh->faces_begin( ) ; curFace != _mesh->faces_end( ) ; curFace++ ) {
        for ( MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter( *curFace ) ; curEdge.is_valid() ; curEdge++ ) {

            FaceHandle fh0 = *curFace;
            EdgeHandle eh = *curEdge;
            HalfedgeHandle heh0 = _mesh->halfedge_handle( eh , 0 );
            HalfedgeHandle heh1 = _mesh->halfedge_handle( eh , 1 );

            FaceHandle fh1 = _mesh->face_handle( heh0 );
            if( fh0.idx() == fh1.idx() ) fh1 = _mesh->face_handle( heh1 );

            float distance = geodesicDistance( _mesh , eh.idx() );
            geodesicDistances[std::make_pair(  fh0.idx() , fh1.idx() )] = distance;

        }
    }
}

float MainWindow::avgAngularDistance() {
    float result = 0.0f;
    for (std::map<std::pair<int , int> , float>::iterator it = angularDistances.begin() ; it != angularDistances.end() ; it++ ) {
        result += it->second;
    }
    return result / angularDistances.size();
}

float MainWindow::avgGeodesicDistances() {
    float result = 0.0f;
    for (std::map<std::pair<int , int> , float>::iterator it = geodesicDistances.begin() ; it != geodesicDistances.end() ; it++ ) {
        result += it->second;
    }
    return result / geodesicDistances.size();
}

void MainWindow::segmentationSimple(MyMesh* _mesh, int k) {
    for(MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++) {
        if(curFace->is_valid()) {
            _mesh->set_color(*curFace, MyMesh::Color(255,0,0));
        }
    }
    displayMesh(_mesh);
}

/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_segmentation_clicked() {
    int k = ui->spinBox_segmentation->value();
    segmentationSimple(&mesh, k);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

/* **** fin de la partie boutons et IHM **** */



/* **** fonctions supplémentaires **** */

// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
    MyMesh::ConstFaceVertexIter fvIt;
    int i = 0;
    for (; fIt!=fEnd; ++fIt)
    {
        fvIt = _mesh->cfv_iter(*fIt);

        triCols[3*i+0] = _mesh->color(*fIt)[0];
        triCols[3*i+1] = _mesh->color(*fIt)[1];
        triCols[3*i+2] = _mesh->color(*fIt)[2];

        triVerts[3*i+0] = _mesh->point(*fvIt)[0];
        triVerts[3*i+1] = _mesh->point(*fvIt)[1];
        triVerts[3*i+2] = _mesh->point(*fvIt)[2];

        triIndiceArray[i] = i;

        i++;
        ++fvIt;

        triCols[3*i+0] = _mesh->color(*fIt)[0];
        triCols[3*i+1] = _mesh->color(*fIt)[1];
        triCols[3*i+2] = _mesh->color(*fIt)[2];

        triVerts[3*i+0] = _mesh->point(*fvIt)[0];
        triVerts[3*i+1] = _mesh->point(*fvIt)[1];
        triVerts[3*i+2] = _mesh->point(*fvIt)[2];

        triIndiceArray[i] = i;

        i++;
        ++fvIt;

        triCols[3*i+0] = _mesh->color(*fIt)[0];
        triCols[3*i+1] = _mesh->color(*fIt)[1];
        triCols[3*i+2] = _mesh->color(*fIt)[2];

        triVerts[3*i+0] = _mesh->point(*fvIt)[0];
        triVerts[3*i+1] = _mesh->point(*fvIt)[1];
        triVerts[3*i+2] = _mesh->point(*fvIt)[2];

        triIndiceArray[i] = i;

        i++;
    }

    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

