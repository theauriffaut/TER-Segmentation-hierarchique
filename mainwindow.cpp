#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <limits>
#include <queue>
#include <time.h>

#define myqDebug() qDebug() << fixed << qSetRealNumberPrecision(8)

auto compWeight = [](const QPair<int, double> &a, const QPair<int, double> &b) {return a.second > b.second; };

// calcul le plus court chemin en taille à l'aide de l'algorithme de dijkstra
double MainWindow::dijkstraDual(int v1, int v2) {

  int NbNodes = dual.getNbVertices();
  //int NbEdges = dual.getNbEdges();

  int StartNode = v1;

  QMap<int, double> Distances;

  for(int i = 0; i < NbNodes; i++) {
      Distances.insert(dual[i].getId(), DBL_MAX);
  }

  Distances[StartNode] = 0;

  QMap<int, int> Parents;
  for(int i = 0; i < NbNodes; i++) {
      Parents.insert(dual[i].getId(), -1);
  }

  priority_queue<QPair<int, double>, QVector<QPair<int, double>>, decltype(compWeight)> Q(compWeight);
  Q.push(QPair<int, double>(StartNode, 0));

  QMap<int, std::map<int, double>> G;
  for(int i = 0; i < NbNodes; i++) {
      G.insert(dual[i].getId(), dual[i].m_adjacencyList);

  }

  while (!Q.empty()) {
    int v = Q.top().first;
    double w = Q.top().second;
    Q.pop();

    if (w <= Distances[v]) {
      for (const auto& i : G[v]) {
        int v2 = i.first;
        double w2 = i.second;

        if (Distances[v] + w2 < Distances[v2]) {
          Distances[v2] = Distances[v] + w2;
          Parents[v2] = v;
          Q.push(QPair<int, double>(v2, Distances[v2]));
        }
      }
    }
  }

  QVector<int> chemin;

  if(Parents[v2] == -1){
      if(v2 != StartNode) {
        return DBL_MAX;
         //qDebug() << "Erreur : Chemin impossible entre le vertex" << v1 << "et le vertex" << v2 << "car composante non connexe";
      } else {
        return 0;
      }

  } else {


    /*chemin.push_back(v2);
    //qDebug() << "Vertex depart :" << v2;
    for (int p = Parents[v2]; p != -1; p = Parents[p]){
<<<<<<< Updated upstream
      qDebug() << " <- " << p << Distances[p] << Parents[p];
      if(p == 0 && Parents[p] == 0){
          return Distances[v2];
      }
=======
      qDebug() << " <- " << p;
>>>>>>> Stashed changes
      chemin.push_back(p);
    }*/

    //qDebug() << "Chemin depuis le vertex" << StartNode << "au vertex" << v2 << "a une taille de" << Distances[v2] << endl;
  }

    return Distances[v2];
}

void MainWindow::dijkstraByREP(MyMesh* _mesh, int IdFaceREP, int k) {
    for(int i = 0; i < patches[k].size(); i++){
        _mesh->property(dist, _mesh->face_handle(patches[k][i])).push_back(dijkstraDual(IdFaceREP, patches[k][i]));

        nbDijkstraDone++;
        ui->progressDijkstra->setValue(nbDijkstraDone);
    }
}

double MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1)
{
    FaceHandle fh0 = _mesh->face_handle ( faceID0 );
    FaceHandle fh1 = _mesh->face_handle ( faceID1 );

    // On récupère les normales, calcule leur produit scalaire et on renvoie l'angle
    // non signé avec le produit scalaire
    OpenMesh::Vec3f normal0 (_mesh->normal ( fh0 ) );
    OpenMesh::Vec3f normal1 (_mesh->normal ( fh1 ) );
    double scalar = normal0 | normal1;

    return acos ( scalar );
}

void MainWindow::computeAngularDistances( MyMesh *_mesh ) {
    angularDistances.clear();
    for ( MyMesh::FaceIter curFace = _mesh->faces_begin( ) ; curFace != _mesh->faces_end( ) ; curFace++ ) {
        for ( MyMesh::FaceFaceIter curNeigh = _mesh->ff_iter( *curFace ) ; curNeigh.is_valid() ; curNeigh++ ) {

            FaceHandle fh0 = *curFace;
            FaceHandle fh1 = *curNeigh;
            double angle = angleFF( _mesh , fh0.idx() , fh1.idx() );
            double coef = angle > M_PI ? 0.1 : 1.0;
            double result = coef * ( 1.0 - cos( angle ) );
            if( isnan( result ) ) result = 0.0;
            angularDistances[std::make_pair( fh0.idx() , fh1.idx() )] = result;

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

double MainWindow::geodesicDistance( MyMesh *_mesh , int edgeID ) {
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
         double length0 = abs( vector0.norm() ) + abs( vector1.norm() );


         vector0 = p2 - center0;
         vector1 = center1 - p2;
         double length1 = abs( vector0.norm() ) + abs( vector1.norm() );

         return length0 < length1;
     } );

     for( int i = 0 ; i != edgePoints.size() ; ++i ) {
         MyMesh::Point p1 = edgePoints[i];
         VectorT<float, 3> vector0 = p1 - center0;
         VectorT<float, 3> vector1 = center1 - p1;
         double length0 = abs( vector0.norm() ) + abs( vector1.norm() );
     }

     MyMesh::Point pointForShortestPath = edgePoints.front();
     VectorT<float, 3> vector0 = pointForShortestPath - center0;
     VectorT<float, 3> vector1 = center1 - pointForShortestPath;
     return ( abs( vector0.norm() ) + abs( vector1.norm() ) );
}

void MainWindow::computeGeodesicDistances( MyMesh *_mesh ) {
    geodesicDistances.clear();
    for ( MyMesh::FaceIter curFace = _mesh->faces_begin( ) ; curFace != _mesh->faces_end( ) ; curFace++ ) {
        for ( MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter( *curFace ) ; curEdge.is_valid() ; curEdge++ ) {

            FaceHandle fh0 = *curFace;
            EdgeHandle eh = *curEdge;
            HalfedgeHandle heh0 = _mesh->halfedge_handle( eh , 0 );
            HalfedgeHandle heh1 = _mesh->halfedge_handle( eh , 1 );

            FaceHandle fh1 = _mesh->face_handle( heh0 );
            if( fh0.idx() == fh1.idx() ) fh1 = _mesh->face_handle( heh1 );

            //Si la face est en bordure, on passe a l'arête voisine suivante si il n'y pas de face opposée liée
            if ( fh1.idx ( ) > _mesh->n_faces ( ) )
                continue;

            double distance = geodesicDistance( _mesh , eh.idx() );
            geodesicDistances[std::make_pair(  fh0.idx() , fh1.idx() )] = distance;

        }
    }
}

void MainWindow::computeWeight( MyMesh *_mesh , double coefGeod , int patch) {
    int progress = 0;
    ui->progressDual->setRange(0, patches[patch].size());
    ui->progressDual->setFormat("%p%");
    ui->progressDual->setValue(progress);

    double avgGeodesic = avgGeodesicDistances();
    double avgAngular = avgAngularDistances();

    for ( MyMesh::FaceIter curFace = _mesh->faces_begin( ) ; curFace != _mesh->faces_end( ) ; curFace++ ) {
        if(_mesh->property(patchId, *curFace) == patch){
            for ( MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter( *curFace ) ; curEdge.is_valid() ; curEdge++ ) {

                FaceHandle fh0 = *curFace;
                EdgeHandle eh = *curEdge;
                HalfedgeHandle heh0 = _mesh->halfedge_handle( eh , 0 );
                HalfedgeHandle heh1 = _mesh->halfedge_handle( eh , 1 );

                FaceHandle fh1 = _mesh->face_handle( heh0 );
                if( fh0.idx() == fh1.idx() ) fh1 = _mesh->face_handle( heh1 );

                //Si la face est en bordure, on passe a l'arête voisine suivante si il n'y pas de face opposée liée
                if ( fh1.idx () > _mesh->n_faces () || _mesh->property(patchId, fh1) != patch){
                    continue;
                }

                std::pair<int , int> key = std::make_pair( fh0.idx() , fh1.idx() );
                double geodesicDistance = geodesicDistances[key];
                double angularDistance = angularDistances[key];

                double weight = ( coefGeod * ( geodesicDistance / avgGeodesic ) ) +
                                ( ( 1 - coefGeod ) * ( angularDistance / avgAngular ) );

                //qDebug() << weight;
                 //Pour le moment
                 //double weight = geodesicDistance / avgGeodesicDistances();

                 dual.addEdge( fh0.idx() , fh1.idx() , weight );
            }
            progress++;
            ui->progressDual->setValue(progress);
        }
    }
}


double MainWindow::avgAngularDistances() {
    double result = 0.0;
    for (std::map<std::pair<int , int> , double>::iterator it = angularDistances.begin() ; it != angularDistances.end() ; it++ ) {
        result += it->second;
    }
    return result / angularDistances.size();
}

double MainWindow::avgGeodesicDistances() {
    double result = 0.0;
    for (std::map<std::pair<int , int> , double>::iterator it = geodesicDistances.begin() ; it != geodesicDistances.end() ; it++ ) {
        result += it->second;
    }
    return result / geodesicDistances.size();
}

void MainWindow::displayGeodesicDistances(){
    for (std::map<std::pair<int , int> , double>::iterator it = geodesicDistances.begin() ; it != geodesicDistances.end() ; it++ ) {
        std::pair<int, int> key = it->first;
        float value = it->second;

        myqDebug() << "La distance géodésique entre la face " << std::get<0>(key)
                  << " et la face " << std::get<1>(key) << " est de " << value << ".";
    }
}

void MainWindow::displayAngularDistances() {
    for (std::map<std::pair<int , int> , double>::iterator it = angularDistances.begin() ; it != angularDistances.end() ; it++ ) {
        std::pair<int, int> key = it->first;
        float value = it->second;

        myqDebug() << "La distance angulaire entre la face " << std::get<0>(key)
                  << " et la face " << std::get<1>(key) << " est de " << value << ".";
    }
}

void MainWindow::computeDirectDistances(MyMesh *_mesh, int patch) {
    int progress = 0;
    ui->progressChoix->setRange(0, patches[patch].size());
    ui->progressChoix->setFormat("%p%");
    ui->progressChoix->setValue(progress);
     for(auto curFace : patches[patch]) {
         for(auto curFace2 : patches[patch]) {
             FaceHandle fh0 = _mesh->face_handle(curFace);
             FaceHandle fh1 = _mesh->face_handle(curFace2);
             if ( fh0.idx( ) == fh1.idx( ) ||
                  directDistances.find( std::make_pair( fh0.idx() , fh1.idx() ) ) != directDistances.end() ||
                  directDistances.find( std::make_pair( fh1.idx() , fh0.idx() ) ) != directDistances.end() ) continue;

             MyMesh::Point center0 = faceGravityCenter( _mesh , fh0.idx( ) );
             MyMesh::Point center1 = faceGravityCenter( _mesh , fh1.idx() );

             VectorT<float, 3> vector = center1 - center0;
             double length = abs( vector.norm() );

             directDistances[std::make_pair(  fh0.idx() , fh1.idx() )] = length;
         }
         progress++;
         ui->progressChoix->setValue(progress);
     }
}

void MainWindow::displayDirectDistances() {
    for (std::map<std::pair<int , int> , double>::iterator it = directDistances.begin() ; it != directDistances.end() ; it++ ) {
        std::pair<int, int> key = it->first;
        float value = it->second;

        myqDebug() << "La distance directe entre la face " << std::get<0>(key)
                  << " et la face " << std::get<1>(key) << " est de " << value << ".";
    }
}

void MainWindow::computeProbabilities(MyMesh *_mesh, QVector<int> IdReps, int k) {
    int progress = 0;
    ui->progressProba->setRange(0, patches[k].size());
    ui->progressProba->setFormat("%p%");
    ui->progressProba->setValue(progress);
    for(int i = 0; i < patches[k].size(); i++){

        QVector<double> distRep = _mesh->property(dist, _mesh->face_handle(patches[k][i]));
        double ProbB = distRep[0] / (distRep[0] + distRep[1]);

        _mesh->property(PB, _mesh->face_handle(patches[k][i])).push_back(1-ProbB);
        _mesh->property(PB, _mesh->face_handle(patches[k][i])).push_back(ProbB);

        nbDijkstraDone++;
        ui->progressDijkstra->setValue(nbDijkstraDone);
        progress++;
        ui->progressProba->setValue(progress);
    }
}

void MainWindow::segmentationSimple(MyMesh* _mesh, int k) {



    colors = { MyMesh::Color(102,0,255), MyMesh::Color(254,231,240), MyMesh::Color(212,115,212), MyMesh::Color(255,0,255), MyMesh::Color(121,248,248), MyMesh::Color(223,109,20),
               MyMesh::Color(115,8,0), MyMesh::Color(1,215,88), MyMesh::Color(240,195,0), MyMesh::Color(255,9,33), MyMesh::Color(231,62,1), MyMesh::Color(4,139,154), MyMesh::Color(135,233,144),
               MyMesh::Color(63,34,4), MyMesh::Color(49,140,231) };

    patches = QVector<QVector<int>>(k);
    currentId = 0;

    ui->progressTotal->setValue(0);
    ui->progressColor->setValue(0);
    ui->progressDual->setValue(0);
    ui->progressDijkstra->setValue(0);
    ui->progressChoix->setValue(0);
    ui->progressProba->setValue(0);

    int nbStepsDone = 0;

    displayMesh(_mesh);
    _mesh->add_property(dist);
    _mesh->add_property(patchId);    
    _mesh->add_property(PB);

    for ( MyMesh::FaceIter curFace = _mesh->faces_begin( ) ; curFace != _mesh->faces_end( ) ; curFace++ ) {
        _mesh->property(patchId, *curFace) = 0;
        patches[0].push_back(curFace->idx());
    }

    computeAngularDistances( _mesh );
    computeGeodesicDistances( _mesh );

    ui->progressTotal->setRange(0, (k-1)*11);
    ui->progressTotal->setFormat("%p%");
    ui->progressTotal->setValue(nbStepsDone);

    for(int i = 1; i < k; i++){

        int chosenPatch = 0;
        for(int i = 0; i < patches.size(); i++){
            if(patches[chosenPatch].size() < patches[i].size()){
                chosenPatch = i;
            }
        }
        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);

        dual.clear();
        computeWeight( _mesh , coefGeod , chosenPatch);
        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);

        directDistances.clear();
        computeDirectDistances(_mesh, chosenPatch);


        std::pair<int,int> reps;
        double max = 0;
        for(auto it : directDistances){
            if(it.second > max){
                max = it.second;
                reps = it.first;
            }
        }

        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);

        //dual.displayGraph();

        QVector<int> REPs = {reps.first, reps.second};

        ui->progressDijkstra->setRange(0, patches[chosenPatch].size());
        ui->progressDijkstra->setFormat("%p%");
        ui->progressDijkstra->setValue(nbDijkstraDone);

        for(int i = 0; i < patches[chosenPatch].size(); i++){
            _mesh->property(dist, _mesh->face_handle(patches[chosenPatch][i])).clear();
        }
        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);

        for(int i = 0; i < 2; i++) {
            dijkstraByREP(_mesh, REPs[i], chosenPatch);
            nbStepsDone++;
            ui->progressTotal->setValue(nbStepsDone);
            nbDijkstraDone = 0;
        }

        for(int i = 0; i < patches[chosenPatch].size(); i++){
            _mesh->property(PB, _mesh->face_handle(patches[chosenPatch][i])).clear();
        }

        computeProbabilities(_mesh, REPs, chosenPatch);
        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);

        /*PA = 1 - PB
        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);*/

        int nbFacesColored = 0;
        ui->progressColor->setRange(0, _mesh->n_faces());
        ui->progressColor->setFormat("%p%");
        ui->progressColor->setValue(nbFacesColored);


        currentId++;

        double distMAX = dijkstraDual(REPs[0], REPs[1]);

        ambiguousFaces.clear();

        QVector<int>::iterator it = patches[chosenPatch].begin();
        while(it != patches[chosenPatch].end()){
            bool alreadyColored = false;

            if(_mesh->property(PB, _mesh->face_handle(*it))[0] >= minProba) {
                alreadyColored = true;
            }
            if(!alreadyColored){
                if(_mesh->property(PB, _mesh->face_handle(*it))[1] >= minProba){
                    qDebug() << *it << "Appartient au patch " << currentId;
                    _mesh->property(patchId, _mesh->face_handle(*it)) = currentId;
                    patches[currentId].push_back(*it);
                    patches[chosenPatch].removeOne(*it);
                    alreadyColored = true;
                } else {
                    qDebug() << *it << "Appartient à personne ";
                    _mesh->property(patchId, _mesh->face_handle(*it)) = -1;
                    ambiguousFaces.push_back(*it);
                    patches[chosenPatch].removeOne(*it);
                }
            } else {
                if(_mesh->property(PB, _mesh->face_handle(*it))[1] >= minProba){
                    qDebug() << *it << "Face ambigue ";
                    _mesh->property(patchId, _mesh->face_handle(*it)) = -1;
                    ambiguousFaces.push_back(*it);
                    patches[chosenPatch].removeOne(*it);
                } else {
                    qDebug() << *it << "Appartient au patch " << chosenPatch;
                    it++;
                }
            }
        }
        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);

        int biggest = chosenPatch;
        if(patches[chosenPatch].size() < patches[currentId].size()){
            biggest = currentId;
        }
        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);

        for(int i = 0; i < ambiguousFaces.size(); i++){
            _mesh->property(patchId, _mesh->face_handle(ambiguousFaces[i])) = biggest;
            patches[biggest].push_back(ambiguousFaces[i]);
        }

        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);
        int nb =0;
        for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
        {
            qDebug() << curFace->idx() << _mesh->property(patchId, *curFace);
            if(_mesh->property(patchId, *curFace) == -1){
                _mesh->set_color(*curFace, MyMesh::Color(255, 0, 0));
                nb++;
            } else if(_mesh->property(patchId, *curFace) == -2){
                _mesh->set_color(*curFace, MyMesh::Color(255, 255, 255));
                nb++;
            } else {
                _mesh->set_color(*curFace, colors[_mesh->property(patchId, *curFace)]);
            }
            nbFacesColored++;
            ui->progressColor->setValue(nbFacesColored);
        }
        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);
        qDebug() << nb;
        for(int i = 0; i < 2; i++) {
            //_mesh->set_color(_mesh->face_handle(REPs[i]), MyMesh::Color(0,255,0));
        }

        qDebug() << patches[chosenPatch].size();
        qDebug() << patches[currentId].size();
        qDebug() << nb+patches[chosenPatch].size()+patches[currentId].size();
        displayMesh(_mesh);
    }

    //double dist = dijkstraDual(3, 135);

    /*Graph g;
    g.addVertex(0);
    g.addVertex(1);
    g.addVertex(2);
    g.addVertex(10);
    g.addVertex(3);
    g.addVertex(0);
    g.addEdge(0, 1, 20.3);
    g.addEdge(1, 3, 10.0);
    g.addEdge(2, 10, 14.0);
    g.addEdge(200, 100, 0.1);*/


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

    mesh.update_normals();

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

