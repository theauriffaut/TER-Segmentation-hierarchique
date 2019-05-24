#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QThread>
#include <limits>
#include <queue>
#include <time.h>
#include <fstream>

#define myqDebug() qDebug() << fixed << qSetRealNumberPrecision(8)

auto compWeight = [](const QPair<int, double> &a, const QPair<int, double> &b) {return a.second > b.second; };

///Calcul le plus court chemin en taille à l'aide de l'algorithme de dijkstra
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
         qDebug() << "Erreur : Chemin impossible entre le vertex" << v1 << "et le vertex" << v2 << "car composante non connexe";
      } else {

        return 0;
      }

  }

    return Distances[v2];
}

///Calcul de l'angle entre deux faces
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

///Calcul de la distance angulaire pour toutes les faces adjacentes
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

///Renvoi le centre de gravité d'une face
MyMesh::Point MainWindow::faceGravityCenter( MyMesh *_mesh , int faceID ) {
    FaceHandle fh = _mesh->face_handle( faceID );
    std::vector<VertexHandle> vertices;
    for( MyMesh::FaceVertexIter curVert = _mesh->fv_iter( fh ) ; curVert.is_valid() ; curVert++ ) {
        vertices.push_back( *curVert );
    }

    MyMesh::Point result = ( _mesh->point( vertices[0] ) + _mesh->point( vertices[1] ) + _mesh->point( vertices[2] ) ) / 3;
    return result;
}

///Discrétise une arête
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

///Calcul de la distance géodésique entre deux faces à partir d'une arête
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

///Calcul de la distance géodésique pour toutes les arêtes
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

///Calcul le dual du mesh et ajoute un poid à chaque arête de ce dual
void MainWindow::computeWeight( MyMesh *_mesh , double coefGeod) {
    int progress = 0;
    ui->progressDual->setRange(0, _mesh->n_faces());
    ui->progressDual->setFormat("%p%");
    ui->progressDual->setValue(progress);

    double avgGeodesic = avgGeodesicDistances();
    double avgAngular = avgAngularDistances();

    for ( MyMesh::FaceIter curFace = _mesh->faces_begin( ) ; curFace != _mesh->faces_end( ) ; curFace++ ) {
            for ( MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter( *curFace ) ; curEdge.is_valid() ; curEdge++ ) {

                FaceHandle fh0 = *curFace;
                EdgeHandle eh = *curEdge;
                HalfedgeHandle heh0 = _mesh->halfedge_handle( eh , 0 );
                HalfedgeHandle heh1 = _mesh->halfedge_handle( eh , 1 );

                FaceHandle fh1 = _mesh->face_handle( heh0 );
                if( fh0.idx() == fh1.idx() )
                    fh1 = _mesh->face_handle( heh1 );

                //Si la face est en bordure, on passe a l'arête voisine suivante si il n'y pas de face opposée liée
                if ( fh1.idx () > _mesh->n_faces ()){
                    continue;
                }

                std::pair<int , int> key = std::make_pair( fh0.idx() , fh1.idx() );
                double geodesicDistance = geodesicDistances[key];
                double angularDistance = angularDistances[key];

                double weight = ( coefGeod * ( geodesicDistance / avgGeodesic ) ) +
                                ( ( 1 - coefGeod ) * ( angularDistance / avgAngular ) );


                 dual.addEdge( fh0.idx() , fh1.idx() , weight );
            }
            progress++;
            ui->progressDual->setValue(progress);

    }
}

///Moyenne de la distance angulaire du mesh
double MainWindow::avgAngularDistances() {
    double result = 0.0;
    for (std::map<std::pair<int , int> , double>::iterator it = angularDistances.begin() ; it != angularDistances.end() ; it++ ) {
        result += it->second;
    }
    return result / angularDistances.size();
}

///Moyenne de la distance géodésique du mesh
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

///Calcul de la distance géodésique entre toutes les faces du mesh (très gourmand)
void MainWindow::computeDirectDistances(MyMesh *_mesh) {
    int progress = 0;
    ui->progressComputeDistance->setRange(0, _mesh->n_faces());
    ui->progressComputeDistance->setFormat("%p%");
    ui->progressComputeDistance->setValue(progress);
     for(MyMesh::FaceIter curFace = _mesh->faces_begin( ) ; curFace != _mesh->faces_end( ) ; curFace++) {
         for(MyMesh::FaceIter curFace2 = _mesh->faces_begin( ) ; curFace2 != _mesh->faces_end( ) ; curFace2++) {
             FaceHandle fh0 = *curFace;
             FaceHandle fh1 = *curFace2;
             if ( directDistances.find( std::make_pair( fh0.idx() , fh1.idx() ) ) != directDistances.end() ||
                  directDistances.find( std::make_pair( fh1.idx() , fh0.idx() ) ) != directDistances.end() ) continue;

             double length = dijkstraDual(fh0.idx(), fh1.idx());

             directDistances[std::make_pair(  fh0.idx() , fh1.idx() )] = length;
         }
         progress++;
         ui->progressComputeDistance->setValue(progress);
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

///Calcul des probabilités d'appartenir à une face représentante pour chaque face d'un patch
void MainWindow::computeProbabilities(MyMesh *_mesh, QVector<int> IdReps, int k) {
    int progress = 0;
    ui->progressProba->setRange(0, patches[k].size());
    ui->progressProba->setFormat("%p%");
    ui->progressProba->setValue(progress);
    for(int i = 0; i < patches[k].size(); i++){
        double ProbB = 0;
        if(directDistances.find( std::make_pair( patches[k][i] , IdReps[0]) ) != directDistances.end()){
            if(directDistances.find( std::make_pair( patches[k][i] , IdReps[1]) ) != directDistances.end()){
                ProbB = directDistances.find( std::make_pair( patches[k][i] , IdReps[0]))->second /
                        (directDistances.find( std::make_pair( patches[k][i] , IdReps[0]))->second + directDistances.find( std::make_pair( patches[k][i] , IdReps[1]))->second);
            } else if(directDistances.find( std::make_pair(IdReps[1],  patches[k][i]) ) != directDistances.end()){
                ProbB = directDistances.find( std::make_pair( patches[k][i] , IdReps[0]))->second /
                        (directDistances.find( std::make_pair( patches[k][i] , IdReps[0]))->second + directDistances.find( std::make_pair( IdReps[1], patches[k][i]))->second);
            }
        } else if(directDistances.find( std::make_pair(IdReps[0], patches[k][i] ) ) != directDistances.end()){
            if(directDistances.find( std::make_pair( patches[k][i] , IdReps[1]) ) != directDistances.end()){
                ProbB = directDistances.find( std::make_pair( IdReps[0], patches[k][i]))->second /
                        (directDistances.find( std::make_pair( IdReps[0], patches[k][i]))->second + directDistances.find( std::make_pair( patches[k][i] , IdReps[1]))->second);
            } else if(directDistances.find( std::make_pair(IdReps[1],  patches[k][i]) ) != directDistances.end()){
                ProbB = directDistances.find( std::make_pair( IdReps[0], patches[k][i]))->second /
                        (directDistances.find( std::make_pair( IdReps[0], patches[k][i]))->second + directDistances.find( std::make_pair( IdReps[1], patches[k][i]))->second);
            }
        }

        _mesh->property(PB, _mesh->face_handle(patches[k][i])).push_back(1-ProbB);
        _mesh->property(PB, _mesh->face_handle(patches[k][i])).push_back(ProbB);

        progress++;
        ui->progressProba->setValue(progress);
    }
}

vector<string> split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

///Rempli le vector directDistance avec les données stocké dans un fichier txt
void MainWindow::on_pushButtonLoadDistance_clicked(){
    fileNameDistance = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Fichier de distance (*.txt)"));

    QStringList nameAll = fileNameDistance.split("/");

    QString name = nameAll.back();
    ui->labelLoadDistance->setText(name);
    ifstream fichier(fileNameDistance.toStdString(), ios::in);

    ui->label->setText("Processing ...");
    ui->label->setStyleSheet("color: #FF9900;");

    if(fichier) {
        directDistances.clear();
        std::string ligne;
        while(getline(fichier, ligne)) {
            std::vector<std::string> values = split(ligne, " ");
            directDistances[std::make_pair(std::stoi(values[0]), std::stoi(values[1]))] = std::stod(values[2]);
        }
        fichier.close();
        ui->label->setText("Finished !");
        ui->label->setStyleSheet("color: #00FF00;");

    } else {
        cerr << "Impossible d'ouvrir le fichier !" << endl;
        ui->label->setText("ERROR !!");
        ui->label->setStyleSheet("color: #FF0000;");
    }
}

///Crée une fichier txt rempli des valeurs de distance géodésique
void MainWindow::on_pushButtonComputeDistance_clicked(){

    dual.clear();

    computeAngularDistances( &mesh );
    computeGeodesicDistances( &mesh );

    /* a mon moi du futur, il faut regarder ces resultats.*/
    computeWeight(&mesh , coefGeod);
    dual.displayGraph();


    directDistances.clear();
    computeDirectDistances(&mesh);
    QStringList name = fileName.split(".obj");
    string finalName = name.at(0).toLocal8Bit().data();
    finalName.append("_" + to_string(coefGeod) + ".txt");
    ofstream fichier(finalName, ios::out | ios::trunc);

    if(fichier) {
        std::map<std::pair<int, int>,double>::iterator it;
        for(it = directDistances.begin(); it != directDistances.end(); it++){
            fichier << it->first.first << " " << it->first.second << " " << it->second << endl;
        }

        fichier.close();
    } else {
        cerr << "Impossible d'ouvrir le fichier !" << endl;
    }

}

///Permet la répartition du patch par propagation de maniere récurente
void MainWindow::faceCourse(int id, int chosenPatch, int currentId, bool firstCourse, MyMesh* _mesh){
    if(_mesh->face_handle(id).is_valid()){
        if(_mesh->property(marked,_mesh->face_handle(id)) == false){
            if(firstCourse){
                if(_mesh->property(patchId, _mesh->face_handle(id)) == -1){

                    if(_mesh->property(PB, _mesh->face_handle(id))[0] >= minProba)
                    {
                        _mesh->property(patchId, _mesh->face_handle(id)) = chosenPatch;
                        _mesh->property(marked, _mesh->face_handle(id)) = true;

                        patches[chosenPatch].push_back(id);
                        ambiguousFaces.removeOne(id);

                        for (MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter(_mesh->face_handle(id)); curEdge.is_valid(); curEdge++) {
                            int fh1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),0)).idx();
                            if(id == fh1){
                                fh1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),1)).idx();
                            }
                            faceCourse(fh1, chosenPatch, currentId, firstCourse, _mesh);
                        }
                    }
                }
            } else {
                if(_mesh->property(patchId, _mesh->face_handle(id)) == -1){

                    _mesh->property(patchId, _mesh->face_handle(id)) = currentId;
                    _mesh->property(marked, _mesh->face_handle(id)) = true;

                    patches[currentId].push_back(id);
                    ambiguousFaces.removeOne(id);

                    for (MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter(_mesh->face_handle(id)); curEdge.is_valid(); curEdge++) {
                        int fh1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),0)).idx();
                        if(id == fh1){
                            fh1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),1)).idx();
                        }
                        faceCourse(fh1, chosenPatch, currentId, firstCourse, _mesh);
                    }
                }
            }
        } else {
            if(!firstCourse){
                if(_mesh->property(patchId, _mesh->face_handle(id)) == chosenPatch){
                    if(_mesh->property(PB, _mesh->face_handle(id))[1] >= minProba){

                        _mesh->property(patchId, _mesh->face_handle(id)) = -1;
                        ambiguousFaces.push_back(id);
                        patches[chosenPatch].removeOne(id);

                        for (MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter(_mesh->face_handle(id)); curEdge.is_valid(); curEdge++) {
                            int fh1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),0)).idx();
                            if(id == fh1){
                                fh1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),1)).idx();
                            }
                            faceCourse(fh1, chosenPatch, currentId, firstCourse, _mesh);
                        }
                    }
                }
            }
        }
    }
}

///Regarde si les faces représentante ont déjà été tiré
bool MainWindow::findREPs(QVector<QPair<int,int>> vec, int repA, int repB){
    if(!vec.empty()){
        for(int i = 0; i < vec.size(); i++){
            if(vec[i] == QPair<int,int>(repA, repB)){
                return true;
            }
        }
    }
    return false;
}

///Renvoie la première occurence de face ambigüe trouvé
int MainWindow::findOneAmbiguous(){
    if(!ambiguousFaces.empty()){
        return ambiguousFaces[0];
    }
    return -1;
}

///Répartition des faces ambigüe
void MainWindow::ambiguousCourse(int id, MyMesh* _mesh){
    if(!_mesh->property(ambiguousMarked, _mesh->face_handle(id)) && _mesh->property(patchId, _mesh->face_handle(id)) == -1){
        _mesh->property(ambiguousMarked, _mesh->face_handle(id)) = true;
        waitingAmbiguous.push_back(id);
        for (MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter(_mesh->face_handle(id)); curEdge.is_valid(); curEdge++) {
            int fh1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),0)).idx();
            if(ambiguousFaces[0] == fh1){
                fh1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),1)).idx();
            }
            if(_mesh->face_handle(fh1).is_valid()){
                ambiguousCourse(fh1, _mesh);
            }

        }
    } else if(!_mesh->property(ambiguousMarked, _mesh->face_handle(id)) && _mesh->property(patchId, _mesh->face_handle(id)) != -1){
        _mesh->property(ambiguousMarked, _mesh->face_handle(id)) = true;
        if(repartAmbiguous.find(_mesh->property(patchId, _mesh->face_handle(id))) != repartAmbiguous.end()){
            repartAmbiguous.at(_mesh->property(patchId, _mesh->face_handle(id)))++;
        } else {
            repartAmbiguous.insert(make_pair(_mesh->property(patchId, _mesh->face_handle(id)), 0));
        }
    }
}

///Cherche le patch avec le poid le plus important (à améliorer)
double MainWindow::patchCourse(int patch){
    double poids = 0;

    for(int i = 0; i < patches[patch].size(); i++){
        std::map<int , double>::iterator it;
        int pos = dual.getVertexById(patches[patch][i]);
        if(pos == -1){
            qDebug() << "ERREUR pos";
        }
        std::map<int , double> map = dual[pos].m_adjacencyList;
        for(it = map.begin(); it != map.end(); it++){
            if(directDistances.find( std::make_pair( patches[patch][i] , it->first)) != directDistances.end()){
                poids += directDistances.find( std::make_pair( patches[patch][i] , it->first))->second;
                nbEdges++;
            } else if(directDistances.find( std::make_pair(it->first, patches[patch][i])) != directDistances.end()){
                poids += directDistances.find( std::make_pair(it->first, patches[patch][i]))->second;
                nbEdges++;
            }
        }

    }
    return poids;
}

///Fonction principale qui permet la segmentation hiérarchique
void MainWindow::segmentationSimple(MyMesh* _mesh, int k) {

    colors = { MyMesh::Color(102,0,255), MyMesh::Color(254,231,240), MyMesh::Color(212,115,212), MyMesh::Color(255,0,255), MyMesh::Color(121,248,248), MyMesh::Color(223,109,20),
               MyMesh::Color(115,8,0), MyMesh::Color(1,215,88), MyMesh::Color(240,195,0), MyMesh::Color(255,9,33), MyMesh::Color(231,62,1), MyMesh::Color(4,139,154), MyMesh::Color(135,233,144),
               MyMesh::Color(63,34,4), MyMesh::Color(49,140,231) };

    patches = QVector<QVector<int>>(k);
    currentId = 0;


    dual.clear();

    computeAngularDistances( &mesh );
    computeGeodesicDistances( &mesh );

    /* a mon moi du futur, il faut regarder ces resultats.*/
    computeWeight(&mesh , coefGeod);

    ui->progressTotal->setValue(0);
    ui->progressColor->setValue(0);
    ui->progressDual->setValue(0);
    ui->progressProba->setValue(0);
    ui->progressComputeDistance->setValue(0);
    ui->labelChoix->setText("Waiting ...");
    ui->labelChoix->setStyleSheet("color: #909090;");

    ui->label->setText("Waiting ...");
    ui->label->setStyleSheet("color: #909090;");

    int nbStepsDone = 0;

    displayMesh(_mesh);
    _mesh->add_property(patchId);    
    _mesh->add_property(marked);
    _mesh->add_property(ambiguousMarked);
    _mesh->add_property(PB);

    for ( MyMesh::FaceIter curFace = _mesh->faces_begin( ) ; curFace != _mesh->faces_end( ) ; curFace++ ) {
        _mesh->property(patchId, *curFace) = 0;
        patches[0].push_back(curFace->idx());
    }


    ui->progressTotal->setRange(0, (k-1)*5);
    ui->progressTotal->setFormat("%p%");
    ui->progressTotal->setValue(nbStepsDone);



    for(int i = 1; i < k; i++){

        ui->labelChoix->setText("Waiting ...");
        ui->labelChoix->setStyleSheet("color: #909090;");

        int chosenPatch = 0;
        int nbFaces = (_mesh->n_faces() * minSizePatch) / 100;

        for(int j = 0; j < patches.size(); j++){
            if(patches[chosenPatch].size() < patches[j].size()){
                chosenPatch = j;
            }
        }




        /*
        patchSizes.clear();
        for(int j = 0; j < patches.size(); j++){
            patchSizes.push_back(0);
        }
        for(int j = 0; j < patches.size(); j++){
            nbEdges = 0;
            if(!patches[j].empty()){
                patchSizes[j] = patchCourse(j);
            }
        }
        bool finish = true;
        for(int j = 0; j < patches.size(); j++){
            if(patches[j].size() >= nbFaces){
                finish = false;
            }
        }
        if(finish){
            qDebug() << "Programme fini avec" << i << "patches.";
            return;
        }
        for(int j = 0; j < patchSizes.size(); j++){
            if(patchSizes[chosenPatch] < patchSizes[j] && patches[j].size() >= nbFaces){
                chosenPatch = j;
            }
        }*/

        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);


        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);

        std::pair<int,int> reps;
        double max = 0;

        for(auto it : directDistances){
            if(_mesh->property(patchId, _mesh->face_handle(it.first.first)) == chosenPatch && _mesh->property(patchId, _mesh->face_handle(it.first.second)) == chosenPatch){
                if(it.second > max){
                    max = it.second;
                    reps = it.first;
                }
            }
        }


        int nbFacesColored = 0;

        ui->progressColor->setRange(0, _mesh->n_faces());
        ui->progressColor->setFormat("%p%");
        ui->progressColor->setValue(nbFacesColored);

        QVector<int> REPs;
 //-------------------------------------------------------------------------------------------

        ui->labelChoix->setText("Calculating ...");
        ui->labelChoix->setStyleSheet("color: #FF9900;");

        currentId++;

        bool complete = false;
        QVector<QPair<int,int>> oldREPs;
        while(!complete){
            ui->progressDual->setRange(0, patches[chosenPatch].size());

            REPs = {reps.first, reps.second};

            for(int i = 0; i < patches[currentId].size(); i++){
                _mesh->property(patchId, _mesh->face_handle(patches[currentId][i])) = chosenPatch;
                patches[chosenPatch].push_back(patches[currentId][i]);
            }
            patches[currentId].clear();

            for(int i = 0; i < patches[chosenPatch].size(); i++){
                _mesh->property(PB, _mesh->face_handle(patches[chosenPatch][i])).clear();
            }

            computeProbabilities(_mesh, REPs, chosenPatch);

            double distMAX;

            if(directDistances.find( std::make_pair( REPs[0], REPs[1]) ) != directDistances.end()){
                distMAX = directDistances.find( std::make_pair( REPs[0], REPs[1]) )->second;
            } else {
                distMAX = directDistances.find( std::make_pair( REPs[1], REPs[0]) )->second;
            }


            ambiguousFaces.clear();

            for(int i = 0; i < patches[currentId].size(); i++){
                _mesh->property(patchId, _mesh->face_handle(patches[currentId][i])) = -1;
                _mesh->property(marked, _mesh->face_handle(patches[currentId][i])) = false;
                ambiguousFaces.push_back(patches[currentId][i]);
            }
            patches[currentId].clear();

            for(int i = 0; i < patches[chosenPatch].size(); i++){
                _mesh->property(patchId, _mesh->face_handle(patches[chosenPatch][i])) = -1;
                _mesh->property(marked, _mesh->face_handle(patches[chosenPatch][i])) = false;
                ambiguousFaces.push_back(patches[chosenPatch][i]);
            }
            patches[chosenPatch].clear();

            faceCourse(REPs[0], chosenPatch, currentId, true, _mesh);
            displayMesh(_mesh);

            faceCourse(REPs[1], chosenPatch, currentId, false, _mesh);

            for(int i = 0; i < ambiguousFaces.size(); i++){
                _mesh->property(ambiguousMarked, _mesh->face_handle(ambiguousFaces[i])) = false;
            }

            int ambiguous = findOneAmbiguous();
            while(ambiguous != -1){
                qDebug() << ambiguous;
                repartAmbiguous.clear();
                waitingAmbiguous.clear();
                ambiguousCourse(ambiguous, _mesh);
                std::map<int,int>::iterator it;
                std::map<int,int>::iterator itMax = repartAmbiguous.begin();
                for (it = repartAmbiguous.begin(); it != repartAmbiguous.end(); it++) {
                    if(itMax->second < it->second){
                        itMax = it;
                    }
                }
                for(int i = 0; i < waitingAmbiguous.size(); i++){
                    _mesh->property(patchId, _mesh->face_handle(waitingAmbiguous[i])) = itMax->first;
                    patches[itMax->first].push_back(waitingAmbiguous[i]);
                    ambiguousFaces.removeOne(waitingAmbiguous[i]);
                }
                ambiguous = findOneAmbiguous();
            }

            std::pair<int, int> newReps;

            int idMinA = reps.first;
            double min = DBL_MAX;
            for(int f = 0; f < patches[chosenPatch].size(); f++){
                double minSum = 0;
                for(int fi = 0; fi < patches[chosenPatch].size(); fi++){
                    double distance = 0;
                    std::map<std::pair<int,int>,double>::iterator it = directDistances.find( std::make_pair( patches[chosenPatch][f] , patches[chosenPatch][fi] ) );
                    if(it != directDistances.end()){
                        distance = it->second;
                    } else {
                        std::map<std::pair<int,int>,double>::iterator it2 = directDistances.find( std::make_pair( patches[chosenPatch][fi] , patches[chosenPatch][f] ) );
                        if(it2 != directDistances.end()){
                            distance = it2->second;
                        }
                    }
                    minSum += _mesh->property(PB, _mesh->face_handle(patches[chosenPatch][fi]))[0] * distance;
                }
                if(min > minSum){
                    min = minSum;
                    idMinA = patches[chosenPatch][f];
                }
            }

            newReps.first = idMinA;

            int idMinB = reps.second;
            min = DBL_MAX;
            for(int f = 0; f < patches[currentId].size(); f++){
                double minSum = 0;
                for(int fi = 0; fi < patches[currentId].size(); fi++){
                    double distance = 0;
                    std::map<std::pair<int,int>,double>::iterator it = directDistances.find( std::make_pair( patches[currentId][f] , patches[currentId][fi] ) );
                    if(it != directDistances.end()){
                        distance = it->second;
                    } else {
                        std::map<std::pair<int,int>,double>::iterator it2 = directDistances.find( std::make_pair( patches[currentId][fi] , patches[currentId][f] ) );
                        if(it2 != directDistances.end()){
                            distance = it2->second;
                        }
                    }
                    minSum += _mesh->property(PB, _mesh->face_handle(patches[currentId][fi]))[1] * distance;
                }
                if(min > minSum){
                    min = minSum;
                    idMinB = patches[currentId][f];
                }
            }


            newReps.second = idMinB;

//            _mesh->set_color(_mesh->face_handle(newReps.first), MyMesh::Color(0,255,0));
//            _mesh->set_color(_mesh->face_handle(newReps.second), MyMesh::Color(0,255,0));
//            _mesh->set_color(_mesh->face_handle(reps.first), colors[chosenPatch]);
//            _mesh->set_color(_mesh->face_handle(reps.second), colors[currentId]);

            qDebug() << "Nouvelles Faces :" << newReps << "Ancienne Face :" << reps;
            if(!findREPs(oldREPs, newReps.first, newReps.second)){
                complete = false;
                oldREPs.push_back(QPair<int,int>(newReps.first, newReps.second));
                reps = newReps;
            } else if(newReps.first == newReps.second){
                qDebug() << "Erreur: faces représentantes égales.";
                break;
            } else {
                complete = true;
                reps = newReps;
                REPs = {newReps.first, newReps.second};
            }



            displayMesh(_mesh);
        }

        ui->labelChoix->setText("Finished !");
        ui->labelChoix->setStyleSheet("color: #00FF00;");
        qDebug() << "Nb k =" << i;
        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);
 //-------------------------------------------------------------------------------------------


        int nb =0;
        for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
        {
            //qDebug() << curFace->idx() << _mesh->property(patchId, *curFace);
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
        //qDebug() << nb;


//        qDebug() << patches[chosenPatch].size();
//        qDebug() << patches[currentId].size();
//        qDebug() << nb+patches[chosenPatch].size()+patches[currentId].size();

/*
        bool modif = false;
        for(int y = 0; y < patches.size(); y++){
            if(patches[y].size() <= nbFaces && patches[y].size() != 0){
                QVector<int> neighb;
                for(int x = 0; x < patches.size(); x++){
                    neighb.push_back(0);
                }
                for(int j = 0; j < patches[y].size(); j++){
                    for (MyMesh::FaceEdgeIter curEdge = _mesh->fe_iter(_mesh->face_handle(patches[y][j])); curEdge.is_valid(); curEdge++) {
                        FaceHandle face1 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),1));
                        FaceHandle face2 = _mesh->face_handle(_mesh->halfedge_handle((* curEdge),0));
                        if(face1.is_valid()){
                            neighb[_mesh->property(patchId, face1)]++;
                        }
                        if(face2.is_valid()){
                            neighb[_mesh->property(patchId, face2)]++;
                        }
                    }
                }

                int newPatch = 0;
                for(int j = 0; j < neighb.size(); j++){
                    if(j != currentId && j != chosenPatch){
                        if(neighb[j] > neighb[newPatch]){
                            newPatch = j;
                        }
                    }
                }

                for(int j = 0; j < patches[y].size(); j++){
                    _mesh->property(patchId, _mesh->face_handle(patches[y][j])) = newPatch;
                    patches[newPatch].push_back(patches[y][j]);
                }

                patches[currentId].clear();
                currentId--;
                i--;
                modif = true;
            }
        }
        if(modif){
            compteur++;
        } else {
            compteur = 0;
        }

        nbStepsDone++;
        ui->progressTotal->setValue(nbStepsDone);

        if(compteur == 5){
            break;
        }*/
        displayMesh(_mesh);
    }
    displayMesh(_mesh);
}

/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_segmentation_clicked() {
    int k = ui->spinBox_segmentation->value();
    coefGeod = ui->doubleSpinBox_coef->value();
    minProba = ui->doubleSpinBox_probMin->value();
    minSizePatch = ui->spinBox_pourcentage->value();

    segmentationSimple(&mesh, k);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);

    ui->doubleSpinBox_coef->setValue(coefGeod);
    ui->doubleSpinBox_probMin->setValue(minProba);
    ui->spinBox_pourcentage->setValue(minSizePatch);
    ui->spinBox_segmentation->setValue(nbPatch);
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

