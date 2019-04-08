#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <limits>
#include <queue>

void displayVertex(MyMesh* _mesh, int id, MyMesh::Color c){
    _mesh->set_color(_mesh->vertex_handle(id), c);
    _mesh->data(_mesh->vertex_handle(id)).thickness = 10;
}

void displayNeighborhoodVertex(MyMesh* _mesh, int id, MyMesh::Color c){
    displayVertex(_mesh, id, c);

    for (MyMesh::VertexEdgeIter curEdge = _mesh->ve_iter(_mesh->vertex_handle(id)); curEdge.is_valid(); curEdge++) {
        _mesh->set_color((* curEdge), c);
        _mesh->data((* curEdge)).thickness = 5;
    }

}

void displayEdge(MyMesh* _mesh, int id, MyMesh::Color c){
    _mesh->set_color(_mesh->edge_handle(id), c);
    _mesh->data(_mesh->edge_handle(id)).thickness = 5;

    displayVertex(_mesh,_mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(id),0)).idx(), c);
    displayVertex(_mesh,_mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(id),1)).idx(), c);
}

void displayNeighborhoodEdge(MyMesh* _mesh, int id, MyMesh::Color c){
    displayEdge(_mesh, id, MyMesh::Color(0, 255, 0));

    FaceHandle face1 = _mesh->face_handle(_mesh->halfedge_handle(_mesh->edge_handle(id),0));
    FaceHandle face2 = _mesh->face_handle(_mesh->halfedge_handle(_mesh->edge_handle(id),1));

    if(face1.is_valid()){
        _mesh->set_color(face1, c);
    }
    if(face2.is_valid()){
        _mesh->set_color(face2, c);
    }

}

void displayFace(MyMesh* _mesh, int id, MyMesh::Color c){
    _mesh->set_color(_mesh->face_handle(id), MyMesh::Color(50, 50, 255));
    _mesh->set_color(_mesh->face_handle(id), MyMesh::Color(50, 50, 255));
    for (MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter(_mesh->face_handle(id)); curEdge.is_valid(); curEdge++) {
        displayEdge(_mesh, (* curEdge).idx(), MyMesh::Color(0, 0, 255));
    }
}

void displayNeighborhoodFace(MyMesh* _mesh, int id, MyMesh::Color c){


    for (MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter(_mesh->face_handle(id)); curEdge.is_valid(); curEdge++) {
        FaceHandle face1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),1));
        FaceHandle face2 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),0));
        if(face1.is_valid()){
            _mesh->set_color(face1, c);
        }
        if(face2.is_valid()){
            _mesh->set_color(face2, c);
        }
    }
    displayFace(_mesh, id, MyMesh::Color(50, 50, 255));
}

void MainWindow::showSelections(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);


    int f = ui->labelFace->text().toInt();
    if(f != -1 && f < _mesh->n_faces()){
        displayFace(_mesh, f, MyMesh::Color(50, 50, 255));

    }

    int e = ui->labelEdge->text().toInt();
    if(e != -1 && e < _mesh->n_edges()){
        displayEdge(_mesh, e, MyMesh::Color(0, 255, 0));
    }

    int v = ui->labelVertex->text().toInt();
    if(v != -1 && v < _mesh->n_vertices()){
        displayVertex(_mesh, v, MyMesh::Color(255, 0, 0));
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}


void MainWindow::showSelectionsNeighborhood(MyMesh* _mesh)
{
    // on réinitialise les couleurs de tout le maillage
    resetAllColorsAndThickness(_mesh);

    int f = ui->labelFace->text().toInt();
    if(f != -1 && f < _mesh->n_faces()){
        displayNeighborhoodFace(_mesh, f, MyMesh::Color(100, 100, 255));

    }

    int e = ui->labelEdge->text().toInt();
    if(e != -1 && e < _mesh->n_edges()){
        displayNeighborhoodEdge(_mesh, e, MyMesh::Color(150, 255, 150));
    }

    int v = ui->labelVertex->text().toInt();
    if(v != -1 && v < _mesh->n_vertices()){
        displayNeighborhoodVertex(_mesh, v, MyMesh::Color(255, 0, 0));
    }
    // on affiche le nouveau maillage
    displayMesh(_mesh);
}



void MainWindow::showBorder(MyMesh* _mesh)
{
    // on réinitialise l'affichage
    resetAllColorsAndThickness(_mesh);

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        if(!_mesh->face_handle(_mesh->halfedge_handle((*curEdge), 0)).is_valid() || !_mesh->face_handle(_mesh->halfedge_handle((*curEdge), 1)).is_valid()){
            _mesh->set_color((*curEdge), MyMesh::Color(255, 0, 255));
            _mesh->set_color((*curEdge), MyMesh::Color(255, 0, 255));
            _mesh->data((*curEdge)).thickness = 5;
            _mesh->data((*curEdge)).thickness = 5;
        }
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

auto compInt = [](const QPair<int, float> &a, const QPair<int, float> &b) {return a.second > b.second; };

// calcul le plus court chemin en taille à l'aide de l'algorithme de dijkstra
QVector<int> dijkstraSize(MyMesh* _mesh, int v1, int v2) {

  int NbNodes = _mesh->n_vertices();
  int NbEdges = _mesh->n_edges();

  QVector<QVector<QPair<int, float>>> G(NbNodes);

  for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++) {
      MyMesh::Point x = _mesh->point(_mesh->to_vertex_handle(_mesh->halfedge_handle((*curEdge), 0)));
      MyMesh::Point y = _mesh->point(_mesh->from_vertex_handle(_mesh->halfedge_handle((*curEdge), 0)));

      float sizeXY = (x - y).norm();
      G[_mesh->from_vertex_handle(_mesh->halfedge_handle(*curEdge, 0)).idx()].push_back(QPair<int, float>(_mesh->to_vertex_handle(_mesh->halfedge_handle(*curEdge, 0)).idx(), sizeXY));
      G[_mesh->from_vertex_handle(_mesh->halfedge_handle(*curEdge, 1)).idx()].push_back(QPair<int, float>(_mesh->to_vertex_handle(_mesh->halfedge_handle(*curEdge, 1)).idx(), sizeXY));
  }


  int StartNode = v1;
  qDebug() << StartNode;

  QVector<float> Distances(NbNodes, FLT_MAX);

  Distances[StartNode] = 0;

  QVector<int> Parents(NbNodes, -1);

  priority_queue<QPair<int, float>, QVector<QPair<int, float>>, decltype(compInt)> Q(compInt);
  Q.push(QPair<int, float>(StartNode, 0));

  while (!Q.empty()) {
    int v = Q.top().first;
    float w = Q.top().second;
    Q.pop();

    if (w <= Distances[v]) {
      for (const auto& i : G[v]) {
        int v2 = i.first;
        float w2 = i.second;

        if (Distances[v] + w2 < Distances[v2]) {
          Distances[v2] = Distances[v] + w2;
          Parents[v2] = v;
          Q.push(QPair<int, float>(v2, Distances[v2]));
        }
      }
    }
  }

  QVector<int> chemin;

  if(Parents[v2] == -1){
       qDebug() << "Erreur : Chemin impossible entre le vertex" << v1 << "et le vertex" << v2 << "car composante non connexe";
  } else {


    chemin.push_back(v2);
    qDebug() << "Vertex depart :" << v2;
    for (int p = Parents[v2]; p != -1; p = Parents[p]){
      qDebug() << " <- " << p;
      chemin.push_back(p);
    }

    qDebug() << "Chemin depuis le vertex" << StartNode << "au vertex" << v2 << "a une taille de" << Distances[v2] << endl;
  }

    return chemin;
}

auto compFloat = [](const QPair<int, int> &a, const QPair<int, int> &b) {return a.second > b.second; };

// calcul le plus court chemin en nombre d'arretes à l'aide de l'algorithme de dijkstra
QVector<int> dijkstraNbr(MyMesh* _mesh, int v1, int v2) {

  int NbNodes = _mesh->n_vertices();
  int NbEdges = _mesh->n_edges();

  QVector<QVector<QPair<int, int>>> G(NbNodes);

  for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++) {
      G[_mesh->from_vertex_handle(_mesh->halfedge_handle(*curEdge, 0)).idx()].push_back(QPair<int, int>(_mesh->to_vertex_handle(_mesh->halfedge_handle(*curEdge, 0)).idx(), 1));
      G[_mesh->from_vertex_handle(_mesh->halfedge_handle(*curEdge, 1)).idx()].push_back(QPair<int, int>(_mesh->to_vertex_handle(_mesh->halfedge_handle(*curEdge, 1)).idx(), 1));
  }


  int StartNode = v1;
  qDebug() << StartNode;

  QVector<int> Distances(NbNodes, INT_MAX);

  Distances[StartNode] = 0;

  QVector<int> Parents(NbNodes, -1);

  priority_queue<QPair<int, int>, QVector<QPair<int, int>>, decltype(compFloat)> Q(compFloat);
  Q.push(QPair<int, int>(StartNode, 0));

  while (!Q.empty()) {
    int v = Q.top().first;
    int w = Q.top().second;
    Q.pop();

    if (w <= Distances[v]) {
      for (const auto& i : G[v]) {
        int v2 = i.first;
        int w2 = i.second;

        if (Distances[v] + w2 < Distances[v2]) {
          Distances[v2] = Distances[v] + w2;
          Parents[v2] = v;
          Q.push(QPair<int, int>(v2, Distances[v2]));
        }
      }
    }
  }

  QVector<int> chemin;

  if(Parents[v2] == -1){
       qDebug() << "Erreur : Chemin impossible entre le vertex" << v1 << "et le vertex" << v2 << "car composante non connexe";
  } else {


    chemin.push_back(v2);
    qDebug() << "Vertex depart :" << v2;
    for (int p = Parents[v2]; p != -1; p = Parents[p]){
      qDebug() << " <- " << p;
      chemin.push_back(p);
    }

    qDebug() << "Chemin depuis le vertex" << StartNode << "au vertex" << v2 << "passe par" << Distances[v2] << "edges" << endl;
  }

    return chemin;
}

void MainWindow::showPath(MyMesh* _mesh, int v1, int v2, bool bySize)
{
    // on réinitialise l'affichage
    resetAllColorsAndThickness(_mesh);

    // point de départ et point d'arrivée en vert et en gros
    _mesh->set_color(_mesh->vertex_handle(v1), MyMesh::Color(0, 255, 0));
    _mesh->set_color(_mesh->vertex_handle(v2), MyMesh::Color(0, 255, 0));
    _mesh->data(_mesh->vertex_handle(v1)).thickness = 12;
    _mesh->data(_mesh->vertex_handle(v2)).thickness = 12;
    QVector<int> chemin;
    if(bySize){
        chemin = dijkstraSize(_mesh, v1, v2);
    } else {
        chemin = dijkstraNbr(_mesh, v1, v2);
    }

    // pour chaque sommet du chemin
    for (int i = 0; i < chemin.size()-1; i++) {
        // on regarde toutes les arretes de ce sommet
        for (MyMesh::VertexEdgeIter curEdge = _mesh->ve_iter(_mesh->vertex_handle(chemin[i])); curEdge.is_valid(); curEdge++) {
            // si l'arret relie le sommet ainsi que le suivant dans le chemin on l'affiche
            if((_mesh->to_vertex_handle(_mesh->halfedge_handle(*curEdge, 1)).idx() == chemin[i+1]) or (_mesh->to_vertex_handle(_mesh->halfedge_handle(*curEdge, 0)).idx() == chemin[i+1])){

                _mesh->set_color((* curEdge), MyMesh::Color(255, 0, 0));
                _mesh->data((* curEdge)).thickness = 5;
            }
        }
       displayVertex(_mesh, chemin[i], MyMesh::Color(0, 255, 0));
    }

    // on affiche le nouveau maillage
    displayMesh(_mesh);
}

/* **** fin de la partie à compléter **** */


/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_bordure_clicked()
{
    showBorder(&mesh);
}

void MainWindow::on_pushButton_voisinage_clicked()
{
    // changement de mode entre avec et sans voisinage
    if(modevoisinage)
    {
        ui->pushButton_voisinage->setText("Mode normal");
        modevoisinage = false;
    }
    else
    {
        ui->pushButton_voisinage->setText("Mode voisinage");
        modevoisinage = true;
    }

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}


void MainWindow::on_pushButton_vertexMoins_clicked()
{
    // mise à jour de l'interface
    if(vertexSelection - 1 >= -1){
        vertexSelection = vertexSelection - 1;
        ui->labelVertex->setText(QString::number(vertexSelection));
    }

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_vertexPlus_clicked()
{
    // mise à jour de l'interface
    if(vertexSelection + 1 < mesh.n_vertices()){
        vertexSelection = vertexSelection + 1;
        ui->labelVertex->setText(QString::number(vertexSelection));
    }

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_edgeMoins_clicked()
{
    // mise à jour de l'interface
    if(edgeSelection - 1 >= -1){
        edgeSelection = edgeSelection - 1;
        ui->labelEdge->setText(QString::number(edgeSelection));
    }

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_edgePlus_clicked()
{
    // mise à jour de l'interface
    if(edgeSelection + 1 < mesh.n_edges()){
        edgeSelection = edgeSelection + 1;
        ui->labelEdge->setText(QString::number(edgeSelection));
    }

    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_faceMoins_clicked()
{
    // mise à jour de l'interface
    if(faceSelection - 1 >= -1){
        faceSelection = faceSelection - 1;
        ui->labelFace->setText(QString::number(faceSelection));
    }
    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_facePlus_clicked()
{
    // mise à jour de l'interface
    if(faceSelection + 1 < mesh.n_faces()){
        faceSelection = faceSelection + 1;
        ui->labelFace->setText(QString::number(faceSelection));
    }
    // on montre la nouvelle selection
    if(!modevoisinage)
        showSelections(&mesh);
    else
        showSelectionsNeighborhood(&mesh);
}

void MainWindow::on_pushButton_afficherCheminNbr_clicked()
{
    // on récupère les sommets de départ et d'arrivée
    int indexV1 = ui->spinBox_v1_chemin->value();
    int indexV2 = ui->spinBox_v2_chemin->value();

    showPath(&mesh, indexV1, indexV2, false);
}

void MainWindow::on_pushButton_afficherCheminSize_clicked()
{
    // on récupère les sommets de départ et d'arrivée
    int indexV1 = ui->spinBox_v1_chemin->value();
    int indexV2 = ui->spinBox_v2_chemin->value();

    showPath(&mesh, indexV1, indexV2, true);
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

    ui->spinBox_v1_chemin->setValue(0);
    ui->spinBox_v1_chemin->setRange(0, mesh.n_vertices()-1);
    ui->spinBox_v2_chemin->setValue(0);
    ui->spinBox_v2_chemin->setRange(0, mesh.n_vertices()-1);
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
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

